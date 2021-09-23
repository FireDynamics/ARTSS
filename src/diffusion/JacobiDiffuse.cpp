/// \file       JacobiDiffuse.cpp
/// \brief      Solves diffusion equation with Jacobian method
/// \details    Solves Diffusion equation \f$ \partial_t \phi_2 = \nu \nabla^2 \phi_2 \f$ via calculated iterations of Jacobi step (dependent on residual/ maximal 1000)
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>

#include "JacobiDiffuse.h"
#include "../utility/Parameters.h"
#include "../boundary/BoundaryController.h"
#include "../Domain.h"
#include "../utility/Utility.h"

JacobiDiffuse::JacobiDiffuse() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    auto params = Parameters::getInstance();

    m_dt = params->get_real("physical_parameters/dt");
    m_dsign = 1.;
    m_w = params->get_real("solver/diffusion/w");

    m_max_iter = static_cast<size_t>(params->get_int("solver/diffusion/max_iter"));
    m_tol_res = params->get_real("solver/diffusion/tol_res");
}

// ============================ Diffuse =====================================
// *****************************************************************************
/// \brief  solves diffusion equation \f$ \partial_t \phi_2 = \nu \ nabla^2 \phi_2 \f$
///     via calculated iterations of Jacobi step (dependent on residual/ maximal iterations)
/// \param  out     output pointer
/// \param  in      input pointer
/// \param  b       source pointer
/// \param  D       diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void JacobiDiffuse::diffuse(Field &out, const Field &in, const Field &b, const real D, bool sync) {
    auto domain = Domain::getInstance();
    auto boundary = BoundaryController::getInstance();
    auto bsize_i = boundary->get_size_inner_list();
    auto bsize_b = boundary->get_size_boundary_list();

    size_t *d_inner_list = boundary->get_inner_list_level_joined();
    size_t *d_bList = boundary->get_boundary_list_level_joined();

#pragma acc data present(out, in, b)
    {
        const real dx = domain->get_dx();  // due to unnecessary parameter passing of *this
        const real dy = domain->get_dy();
        const real dz = domain->get_dz();

        const real reciprocal_dx = 1. / dx;  // due to unnecessary parameter passing of *this
        const real reciprocal_dy = 1. / dy;
        const real reciprocal_dz = 1. / dz;

        const real alpha_x = D * m_dt * reciprocal_dx * reciprocal_dx;  // due to better pgi handling of scalars (instead of arrays)
        const real alpha_y = D * m_dt * reciprocal_dy * reciprocal_dy;
        const real alpha_z = D * m_dt * reciprocal_dz * reciprocal_dz;

        const real reciprocal_beta = (1. + 2. * (alpha_x + alpha_y + alpha_z));

        const real dsign = m_dsign;
        const real w = m_w;

        size_t it = 0;
        const size_t max_it = m_max_iter;
        const real tol_res = m_tol_res;
        real sum;
        real res = 1.;

        while (res > tol_res && it < max_it) {
            JacobiStep(out, in, b, alpha_x, alpha_y, alpha_z, reciprocal_beta, dsign, w, sync);
            boundary->apply_boundary(out, sync);

            sum = 0.;

#pragma acc parallel loop independent present(out, in, d_inner_list[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_inner_list[j];
                res = reciprocal_beta * (out[i] - in[i]);  // = reciprocal_beta*(beta*(b - sum_j!=i(alpha*in))) - reciprocal_beta*in = b - sum(alpha*in) = b - A*x(k)
                sum += res * res;
            }
            // info: in nvvp profile 8byte size copy from to device to/from pageable due to sum!

#pragma acc wait
            res = sqrt(sum);
            it++;

            // swap (no pointer swap due to uncontrolled behavior in TimeIntegration Update)
#pragma acc parallel loop independent present(out, in, d_inner_list[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_inner_list[j];
                in[i] = out[i];
            }

#pragma acc parallel loop independent present(out, in, d_bList[:bsize_b]) async
            for (size_t j = 0; j < bsize_b; ++j) {
                const size_t i = d_bList[j];
                in[i] = out[i];
            }
        }

        if (it % 2 != 0)  // swap necessary when odd number of iterations
        {
#pragma acc parallel loop independent present(out, in, d_inner_list[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_inner_list[j];
                out[i] = in[i];
            }

#pragma acc parallel loop independent present(out, in, d_bList[:bsize_b]) async
            for (size_t j = 0; j < bsize_b; ++j) {
                const size_t i = d_bList[j];
                out[i] = in[i];
            }
        }

        if (sync) {
#pragma acc wait
        }

#ifndef BENCHMARKING
        m_logger->info("Number of iterations: {}", it);
        m_logger->info("Jacobi ||res|| = {:0.5e}", res);
#endif
    }
}

// ======================= Turbulent version ================================
// ************************************************************************
/// \brief  solves turbulent diffusion equation \f$ \partial_t \phi_2 = \nu \ nabla^2 \phi_2 \f$
///     via calculated iterations of Jacobi step (dependent on residual/ maximal iterations)
/// \param  out     output pointer
/// \param  in      input pointer
/// \param  b       source pointer
/// \param  D     diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  EV      turbulent diffusion coefficient (eddy viscosity)
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void JacobiDiffuse::diffuse(
        Field &out, const Field &in, Field const &b,
        real const D, Field const &EV, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    auto boundary = BoundaryController::getInstance();

    size_t *d_inner_list = boundary->get_inner_list_level_joined();
    size_t *d_bList = boundary->get_boundary_list_level_joined();

    auto bsize_i = boundary->get_size_inner_list();
    auto bsize_b = boundary->get_size_boundary_list();

#pragma acc data present(out, in, b, EV)
    {
        const real dx = domain->get_dx();
        const real dy = domain->get_dy();
        const real dz = domain->get_dz();

        const real reciprocal_dx = 1. / dx;
        const real reciprocal_dy = 1. / dy;
        const real reciprocal_dz = 1. / dz;

        real dt = m_dt;

        real alpha_x, alpha_y, alpha_z, reciprocal_beta;  // calculated in JacobiStep!

        const real dsign = m_dsign;
        const real w = m_w;

        size_t it = 0;
        const size_t max_it = m_max_iter;
        const real tol_res = m_tol_res;
        real sum;
        real res = 1.;

        while (res > tol_res && it < max_it) {
            JacobiStep(out, in, b, dsign, w, D, EV, dt, sync);
            boundary->apply_boundary(out, sync);

            sum = 0.;

#pragma acc parallel loop independent present(out, b, EV, d_inner_list[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_inner_list[j];
                alpha_x = (D + EV[i]) * dt * reciprocal_dx * reciprocal_dx;
                alpha_y = (D + EV[i]) * dt * reciprocal_dy * reciprocal_dy;
                alpha_z = (D + EV[i]) * dt * reciprocal_dz * reciprocal_dz;
                reciprocal_beta = (1. + 2. * (alpha_x + alpha_y + alpha_z));

                res = reciprocal_beta * (out[i] - in[i]);
                sum += res * res;
            }
            // info: in nvvp profile 8byte size copy from to device to/from pageable due to sum!
#pragma acc wait
            res = sqrt(sum);
            it++;

// swap (no pointer swap due to uncontrolled behavior in TimeIntegration Update)
#pragma acc parallel loop independent present(out, in, d_inner_list[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_inner_list[j];
                in[i] = out[i];
            }
#pragma acc parallel loop independent present(out, in, d_bList[:bsize_b]) async
            for (size_t j = 0; j < bsize_b; ++j) {
                const size_t i = d_bList[j];
                in[i] = out[i];
            }

        } //end while

        if (it % 2 != 0)// swap necessary when odd number of iterations
        {
#pragma acc parallel loop independent present(out, in, d_inner_list[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_inner_list[j];
                out[i] = in[i];
            }
#pragma acc parallel loop independent present(out, in, d_bList[:bsize_b]) async
            for (size_t j = 0; j < bsize_b; ++j) {
                const size_t i = d_bList[j];
                out[i] = in[i];
            }
        }

        if (sync) {
#pragma acc wait
        }
#ifndef BENCHMARKING
        m_logger->info("Number of iterations: {}", it);
        m_logger->info("Jacobi ||res|| = {:0.5e}", res);
#endif
    }
}

// ======================= Jacobian stencil ================================
// ************************************************************************
/// \brief  performs one (weighted) Jacobi step \f$ x_i^{(n+1)}:=\frac1{a_{ii}}\left(b_i-\sum_{j\not=i} a_{ij}\cdot x_j^{(n)}\right), \, i=0,\dots,N_x-2 \f$
/// \param  out   output pointer
/// \param  in    input pointer
/// \param  b     source pointer
/// \param  alpha 3-dimensional array;
///           \f$ (1/dx^2, 1/dy^2, 1/dz^2)\f$ for pressure,
///           \f$ (\nu\cdot dt\cdot 1/dx^2, \nu\cdot dt\cdot 1/dy^2, \nu\cdot dt\cdot 1/dz^2)\f$  for velocity,
///           \f$ (\kappa\cdot dt\cdot 1/dx^2, \kappa\cdot dt\cdot 1/dy^2, \kappa\cdot dt\cdot 1/dz^2)\f$ for temperature
/// \param  beta  3-dimensional array;
///           \f$ 1./(2\cdot(\alpha_0 + \alpha_1 + \alpha_2) + 0.5)\f$ for pressure
///           \f$ 1/(2\cdot(\alpha_0 + \alpha_1 + \alpha_2) + 1)\f$ for velocity
/// \param  dsign sign (\a -1. for pressure, \a 1. else)
/// \param  w   weight (1. - diffusion, 2./3. - multigrid)
/// \param  sync  synchronous kernel launching (true, default: false)
// ***************************************************************************************
void JacobiDiffuse::JacobiStep(
        Field &out, Field const &in, Field const &b,
        real const alpha_x, real const alpha_y, real const alpha_z,
        real const reciprocal_beta, real const dsign, real const w, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(out.get_level()); //due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny(out.get_level());

    auto boundary = BoundaryController::getInstance();

    size_t *d_inner_list = boundary->get_inner_list_level_joined();
    auto bsize_i = boundary->get_size_inner_list();

    size_t neighbour_i = 1;
    size_t neighbour_j = Nx;
    size_t neighbour_k = Nx * Ny;
#pragma acc parallel loop independent present(out, in, b, d_inner_list[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_inner_list[j];
        real out_h = (dsign * b[i]
                + alpha_x * (in[i + neighbour_i] + in[i - neighbour_i])
                + alpha_y * (in[i + neighbour_j] + in[i - neighbour_j])
                + alpha_z * (in[i + neighbour_k] + in[i - neighbour_k])) / reciprocal_beta;
        out[i] = (1 - w) * in[i] + w * out_h;
    }

    if (sync) {
#pragma acc wait
    }
}

// =============== Multigrid version for Jacobian stencil ===============
// ************************************************************************
/// \brief  performs one (weighted) Jacobi step at multigrid level, \f$ x_i^{(n+1)}:=\frac1{a_{ii}}\left(b_i-\sum_{j\not=i} a_{ij}\cdot x_j^{(n)}\right), \, i=0,\dots,N_x-2 \f$
/// \param  level multigrid level
/// \param  out   output pointer
/// \param  in    input pointer
/// \param  b     source pointer
/// \param  alpha 3-dimensional array;
///           \f$ (1/dx^2, 1/dy^2, 1/dz^2)\f$ for pressure,
///           \f$ (\nu\cdot dt\cdot 1/dx^2, \nu\cdot dt\cdot 1/dy^2, \nu\cdot dt\cdot 1/dz^2)\f$  for velocity,
///           \f$ (\kappa\cdot dt\cdot 1/dx^2, \kappa\cdot dt\cdot 1/dy^2, \kappa\cdot dt\cdot 1/dz^2)\f$ for temperature
/// \param  beta  3-dimensional array;
///           \f$ 1./(2\cdot(\alpha_0 + \alpha_1 + \alpha_2) + 0.5)\f$ for pressure
///           \f$ 1/(2\cdot(\alpha_0 + \alpha_1 + \alpha_2) + 1)\f$ for velocity
/// \param  dsign sign (\a -1. for pressure, \a 1. else)
/// \param  w   weight (1. - diffusion, 2./3. - multigrid)
/// \param  sync  synchronous kernel launching (true, default: false)
// ***************************************************************************************
void JacobiDiffuse::JacobiStep(
        size_t level, Field &out, Field const &in, Field const &b,
        real const alpha_x, real const alpha_y, real const alpha_z,
        real const beta, real const dsign, real const w, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(level); //due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny(level);

    auto boundary = BoundaryController::getInstance();

    size_t *d_inner_list = boundary->get_inner_list_level_joined();
    size_t start_i = boundary->get_inner_list_level_joined_start(level);
    size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;

    size_t neighbour_i = 1;
    size_t neighbour_j = Nx;
    size_t neighbour_k = Nx * Ny;
#pragma acc parallel loop independent present(out, in, b, d_inner_list[start_i:(end_i-start_i)]) async
    for (size_t j = start_i; j < end_i; ++j) {
        const size_t i = d_inner_list[j];
        real out_h = beta * (dsign * b[i]
                + alpha_x * (in[i + neighbour_i] + in[i - neighbour_i])
                + alpha_y * (in[i + neighbour_j] + in[i - neighbour_j])
                + alpha_z * (in[i + neighbour_k] + in[i - neighbour_k]));
        out[i] = (1 - w) * in[i] + w * out_h;
    }

    if (sync) {
#pragma acc wait
    }
}

// =============== Turbulent version for Jacobian stencil ===============
// ************************************************************************
/// \brief  performs one (weighted) Jacobi step for turbulent diffusion, \f$ x_i^{(n+1)}:=\frac1{a_{ii}}\left(b_i-\sum_{j\not=i} a_{ij}\cdot x_j^{(n)}\right), \, i=0,\dots,N_x-2 \f$
/// \param  out   output pointer
/// \param  in    input pointer
/// \param  b     source pointer
/// \param  dsign sign (\a -1. for pressure, \a 1. else)
/// \param  w   weight (1. - diffusion, 2./3. - multigrid)
/// \param  D   diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  EV    turbulent diffusion coefficient (eddy viscosity)
/// \param  dt    time step
/// \param  sync  synchronous kernel launching (true, default: false)
// ***************************************************************************************
void JacobiDiffuse::JacobiStep(
        Field &out, Field const &in, Field const &b,
        real const dsign, real const w, real const D,
        Field const &EV, real const dt, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx();
    const size_t Ny = domain->get_Ny();

    const real dx = domain->get_dx();
    const real dy = domain->get_dy();
    const real dz = domain->get_dz();

    const real reciprocal_dx = 1. / dx;
    const real reciprocal_dy = 1. / dy;
    const real reciprocal_dz = 1. / dz;

    auto boundary = BoundaryController::getInstance();

    size_t *d_inner_list = boundary->get_inner_list_level_joined();
    auto bsize_i = boundary->get_size_inner_list();

    size_t neighbour_i = 1;
    size_t neighbour_j = Nx;
    size_t neighbour_k = Nx * Ny;
#pragma acc parallel loop independent present(out, in, b, EV, d_inner_list[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_inner_list[j];

        auto aX = (D + EV[i]) * dt * reciprocal_dx * reciprocal_dx;
        auto aY = (D + EV[i]) * dt * reciprocal_dy * reciprocal_dy;
        auto aZ = (D + EV[i]) * dt * reciprocal_dz * reciprocal_dz;

        auto rb = (1. + 2. * (aX + aY + aZ));
        auto bb = 1. / rb;

        auto d_aX = aX * (in[i + neighbour_i] + in[i - neighbour_i]) ;
        auto d_aY = aY * (in[i + neighbour_j] + in[i - neighbour_j]) ;
        auto d_aZ = aZ * (in[i + neighbour_k] + in[i - neighbour_k]);

        real out_h = bb * (dsign * b[i] + d_aX + d_aY + d_aZ);
        out[i] = (1 - w) * in[i] + w * out_h;
    }

    if (sync) {
#pragma acc wait
    }
}
