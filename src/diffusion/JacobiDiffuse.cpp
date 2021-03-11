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
/// \param  D     diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void JacobiDiffuse::diffuse(Field &out, Field &in, Field const &b,
        real const D, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    FieldType type = out.get_type();

    auto d_out = out.data;

    auto boundary = BoundaryController::getInstance();

    auto bsize_i = boundary->getSize_innerList();
    auto bsize_b = boundary->getSize_boundaryList();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t *d_bList = boundary->get_boundaryList_level_joined();

#pragma acc data present(out, in, b)
    {
        const real dx = domain->get_dx(out.get_level());  // due to unnecessary parameter passing of *this
        const real dy = domain->get_dy(out.get_level());
        const real dz = domain->get_dz(out.get_level());

        const real rdx = 1. / dx;  // due to unnecessary parameter passing of *this
        const real rdy = 1. / dy;
        const real rdz = 1. / dz;

        const real alphaX = D * m_dt * rdx * rdx;  // due to better pgi handling of scalars (instead of arrays)
        const real alphaY = D * m_dt * rdy * rdy;
        const real alphaZ = D * m_dt * rdz * rdz;

        const real rbeta = (1. + 2. * (alphaX + alphaY + alphaZ));
        const real beta = 1. / rbeta;

        const real dsign = m_dsign;
        const real w = m_w;

        size_t it = 0;
        const size_t max_it = m_max_iter;
        const real tol_res = m_tol_res;
        real sum;
        real res = 1.;

        while (res > tol_res && it < max_it) {
            JacobiStep(out, in, b, alphaX, alphaY, alphaZ, rbeta, dsign, w, sync);
            boundary->applyBoundary(d_out, type, sync);

            sum = 0.;

#pragma acc parallel loop independent present(out, in, d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
                res = rbeta * (out[i] - in[i]);  // = rbeta*(beta*(b - sum_j!=i(alpha*in))) - rbeta*in = b - sum(alpha*in) = b - A*x(k)
                sum += res * res;
            }

// info: in nvvp profile 8byte size copy from to device to/from pageable due to sum!

#pragma acc wait
            res = sqrt(sum);
            it++;

// swap (no pointer swap due to uncontrolled behavior in TimeIntegration Update)
#pragma acc parallel loop independent present(out, in, d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
                in[i] = out[i];
            }

#pragma acc parallel loop independent present(out, in, d_bList[:bsize_b]) async
            for (size_t j = 0; j < bsize_b; ++j) {
                const size_t i = d_bList[j];
                in[i] = out[i];
            }
        }

        if (it % 2 != 0) // swap necessary when odd number of iterations
        {
#pragma acc parallel loop independent present(out, in, d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
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
    }  // end data region
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
void JacobiDiffuse::diffuse(Field &out, Field &in, Field const &b,
        real const D, Field const &EV, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    auto bsize = domain->get_size(out.get_level());
    FieldType type = out.get_type();

    auto d_out = out.data;

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t *d_bList = boundary->get_boundaryList_level_joined();

    auto bsize_i = boundary->getSize_innerList();
    auto bsize_b = boundary->getSize_boundaryList();

#pragma acc data present(d_out[:bsize], d_in[:bsize], d_b[:bsize], d_EV[:bsize])
    {
        const real dx = domain->get_dx(out.get_level()); //due to unnecessary parameter passing of *this
        const real dy = domain->get_dy(out.get_level());
        const real dz = domain->get_dz(out.get_level());

        const real rdx = 1. / dx; //due to unnecessary parameter passing of *this
        const real rdy = 1. / dy;
        const real rdz = 1. / dz;

        real dt = m_dt;

        real alphaX, alphaY, alphaZ, rbeta; //calculated in JacobiStep!

        const real dsign = m_dsign;
        const real w = m_w;

        size_t it = 0;
        const size_t max_it = m_max_iter;
        const real tol_res = m_tol_res;
        real sum;
        real res = 1.;

        while (res > tol_res && it < max_it) {
            JacobiStep(out, in, b, dsign, w, D, EV, dt, sync);
            boundary->applyBoundary(d_out, type, sync);

            sum = 0.;

#pragma acc parallel loop independent present(out, b, EV, d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
                alphaX = (D + EV[i]) * dt * rdx * rdx;
                alphaY = (D + EV[i]) * dt * rdy * rdy;
                alphaZ = (D + EV[i]) * dt * rdz * rdz;
                rbeta = (1. + 2. * (alphaX + alphaY + alphaZ));

                res = rbeta * (out[i] - in[i]);
                sum += res * res;
            }

// info: in nvvp profile 8byte size copy from to device to/from pageable due to sum!

#pragma acc wait
            res = sqrt(sum);
            it++;

// swap (no pointer swap due to uncontrolled behavior in TimeIntegration Update)
#pragma acc parallel loop independent present(out, in, d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
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
#pragma acc parallel loop independent present(out, in, d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
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
void JacobiDiffuse::JacobiStep(Field &out, Field const &in, Field const &b,
        real const alphaX, real const alphaY, real const alphaZ,
        real const rbeta, real const dsign, real const w, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(out.get_level()); //due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny(out.get_level());

    auto bsize = domain->get_size(out.get_level());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

#pragma acc parallel loop independent present(out, in, b, d_iList[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
        real out_h = (dsign * b[i] + alphaX * (in[i + 1] + in[i - 1])
                     + alphaY * (in[i + Nx] + in[i - Nx])
                     + alphaZ * (in[i + Nx * Ny] + in[i - Nx * Ny])) / rbeta;
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
void JacobiDiffuse::JacobiStep(size_t level, Field &out, Field const &in, Field const &b,
        real const alphaX, real const alphaY, real const alphaZ,
        real const beta, real const dsign, real const w, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(level); //due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny(level);

    auto bsize = domain->get_size(out.get_level());
    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t start_i = boundary->get_innerList_level_joined_start(level);
    size_t end_i = boundary->get_innerList_level_joined_end(level) + 1;

#pragma acc parallel loop independent present(out, in, b, d_iList[start_i:(end_i-start_i)]) async
    for (size_t j = start_i; j < end_i; ++j) {
        const size_t i = d_iList[j];
        real out_h = beta * (dsign * b[i] + alphaX * (in[i + 1] + in[i - 1])
                     + alphaY * (in[i + Nx] + in[i - Nx])
                     + alphaZ * (in[i + Nx * Ny] + in[i - Nx * Ny]));
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
void JacobiDiffuse::JacobiStep(Field &out, Field const &in, Field const &b,
        real const dsign, real const w, real const D,
        Field const &EV, real const dt, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(out.get_level()); //due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny(out.get_level());

    const real dx = domain->get_dx(out.get_level()); //due to unnecessary parameter passing of *this
    const real dy = domain->get_dy(out.get_level());
    const real dz = domain->get_dz(out.get_level());

    const real rdx = 1. / dx; //due to unnecessary parameter passing of *this
    const real rdy = 1. / dy;
    const real rdz = 1. / dz;

    real aX, aY, aZ, bb, rb; //multipliers calculated

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

#pragma acc parallel loop independent present(out, in, b, EV, d_iList[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
        aX = (D + EV[i]) * dt * rdx * rdx;
        aY = (D + EV[i]) * dt * rdy * rdy;
        aZ = (D + EV[i]) * dt * rdz * rdz;

        rb = (1. + 2. * (aX + aY + aZ));
        bb = 1. / rb;

        real out_h = bb * (dsign * b[i] + aX * (in[i + 1] + in[i - 1])
                     + aY * (in[i + Nx] + in[i - Nx])
                     + aZ * (in[i + Nx * Ny] + in[i - Nx * Ny]));
        out[i] = (1 - w) * in[i] + w * out_h;
    }

    if (sync) {
#pragma acc wait
    }
}
