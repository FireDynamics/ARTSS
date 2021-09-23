/// \file         ColoredGaussSeidelDiffuse.cpp
/// \brief        solves Diffusion equation \f$ \partial_t \phi_2 = \nu \nabla^2 \phi_2 \f$
///               via calculated iterations of CGS step (dependent on residual/ maximal 1000)
/// \date         June 21, 2016
/// \author       My Linh Wuerzburger
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.


#include <cmath>

#include "ColoredGaussSeidelDiffuse.h"
#include "../boundary/BoundaryController.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../utility/Utility.h"


// ========================== Constructor =================================
ColoredGaussSeidelDiffuse::ColoredGaussSeidelDiffuse() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    auto params = Parameters::getInstance();

    m_dt = params->get_real("physical_parameters/dt");
    m_dsign = 1.;
    m_w = params->get_real("solver/diffusion/w");

    if (params->get("solver/diffusion/type") == "ColoredGaussSeidel") {
        m_max_iter = static_cast<size_t>(params->get_int("solver/diffusion/max_iter"));
        m_tol_res = params->get_real("solver/diffusion/tol_res");
    } else {
        m_max_iter = 10000;
        m_tol_res = 1e-16;
    }
}

// ========================== Diffuse =================================
// ************************************************************************
/// \brief  solves Diffusion equation \f$ \partial_t \phi_2 = \nu \ nabla^2 \phi_2 \f$
///     via calculated iterations of CGS step (dependent on residual/ maximal iterations)
/// \param  out     output pointer
/// \param  in      input pointer
/// \param  b       source pointer
/// \param  D     diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::diffuse(Field &out, const Field &, Field const &b, const real D, bool sync) {
    auto domain = Domain::getInstance();
    auto boundary = BoundaryController::getInstance();

    auto bsize_i = boundary->get_size_inner_list();
    auto bsize_b = boundary->get_size_boundary_list();

    size_t* d_inner_list = boundary->get_inner_list_level_joined();
    size_t* d_boundary_list = boundary->get_boundary_list_level_joined();

//#pragma acc data present(d_out[:bsize], d_b[:bsize])
{
    const size_t Nx = domain->get_Nx();  // due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny();

    const real dx = domain->get_dx();  // due to unnecessary parameter passing of *this
    const real dy = domain->get_dy();
    const real dz = domain->get_dz();

    const real reciprocal_dx = 1. / dx;
    const real reciprocal_dy = 1. / dy;
    const real reciprocal_dz = 1. / dz;

    const real alpha_x = D * m_dt * reciprocal_dx * reciprocal_dx;  // due to better pgi handling of scalars (instead of arrays)
    const real alpha_y = D * m_dt * reciprocal_dy * reciprocal_dy;
    const real alpha_z = D * m_dt * reciprocal_dz * reciprocal_dz;

    const real rbeta    = (1. + 2. * (alpha_x + alpha_y + alpha_z));
    const real beta     = 1. / rbeta;

    const real dsign    = m_dsign;
    const real w        = m_w;

    size_t it           = 0;
    const size_t max_it = m_max_iter;
    const real tol_res  = m_tol_res;

    real sum;
    real res = 1.;

    const size_t neighbour_i = 1;
    const size_t neighbour_j = Nx;
    const size_t neighbour_k = Nx * Ny;
    while (res > tol_res && it < max_it) {
        colored_gauss_seidel_step(out, b, alpha_x, alpha_y, alpha_z, beta, dsign, w, sync);
        boundary->apply_boundary(out, sync);

        sum = 0;

//#pragma acc parallel loop independent present(d_out[:bsize], d_b[:bsize], d_inner_list[:bsize_i]) async
        for (size_t j = 0; j < bsize_i; ++j){
        const size_t index = d_inner_list[j];
            res = (- alpha_x * (out[index + neighbour_i] + out[index - neighbour_i])
                   - alpha_y * (out[index + neighbour_j] + out[index - neighbour_j])
                   - alpha_z * (out[index + neighbour_k] + out[index - neighbour_k])
                   + rbeta * out[index] - b[index]);
            // TODO: find a way to exclude obstacle indices! Jacobi now uses res=D*(out-in)
            sum += res * res;
        }

//#pragma acc wait
        res = sqrt(sum);
        it++;
    } //end while

    if (sync) {
//#pragma acc wait
    }

#ifndef BENCHMARKING
    m_logger->info("Number of iterations: {}", it);
    m_logger->info("Colored Gauss-Seidel ||res|| = {:.5e}", res);
#endif
} //end data region
}

// ===================== Turbulent version ============================
// ************************************************************************
/// \brief  solves Diffusion equation \f$ \partial_t \phi_2 = \nu_{eff} \ nabla^2 \phi_2 \f$
///         via calculated iterations of CGS step (dependent on residual/ maximal iterations)
/// \param  out         output pointer
/// \param  in          input pointer
/// \param  b           source pointer
/// \param  D           diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  EV          eddy viscosity (nu_turb)
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// ************************************************************************
void ColoredGaussSeidelDiffuse::diffuse(Field &out, const Field &, const Field &b,
        const real D, const Field &EV, bool sync) {
    auto domain = Domain::getInstance();
    auto boundary = BoundaryController::getInstance();

    auto bsize_i = boundary->get_size_inner_list();
    auto bsize_b = boundary->get_size_boundary_list();

    size_t* d_inner_list = boundary->get_inner_list_level_joined();
    size_t* d_boundary_list = boundary->get_boundary_list_level_joined();

//#pragma acc data present(d_out[:bsize], d_b[:bsize], d_EV[:bsize])
{
    const size_t Nx = domain->get_Nx();
    const size_t Ny = domain->get_Ny();

    const real dx = domain->get_dx();
    const real dy = domain->get_dy();
    const real dz = domain->get_dz();

    const real reciprocal_dx = 1. / dx;
    const real reciprocal_dy = 1. / dy;
    const real reciprocal_dz = 1. / dz;

    real dt = m_dt;

    real alpha_x, alpha_y, alpha_z, rbeta;  // calculated in colored_gauss_seidel_step!

    const real dsign = m_dsign;
    const real w = m_w;

    size_t it = 0;
    const size_t max_it = m_max_iter;
    const real tol_res = m_tol_res;
    real sum;
    real res = 1.;

    const size_t neighbour_i = 1;
    const size_t neighbour_j = Nx;
    const size_t neighbour_k = Nx * Ny;
    while (res > tol_res && it < max_it) {
        colored_gauss_seidel_step(out, b, dsign, w, D, EV, dt, sync);
        boundary->apply_boundary(out, sync);

        sum = 0;

//#pragma acc parallel loop independent present(d_out[:bsize], d_b[:bsize], d_EV[:bsize], d_inner_list[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j){
        const size_t i = d_inner_list[j];
        alpha_x = (D + EV[i]) * dt * reciprocal_dx * reciprocal_dx;
        alpha_y = (D + EV[i]) * dt * reciprocal_dy * reciprocal_dy;
        alpha_z = (D + EV[i]) * dt * reciprocal_dz * reciprocal_dz;
        rbeta = (1. + 2. * (alpha_x + alpha_y + alpha_z));

            res = (- alpha_x * (out[i + neighbour_i] + out[i - neighbour_i])
                   - alpha_y * (out[i + neighbour_j] + out[i - neighbour_j])
                   - alpha_z * (out[i + neighbour_k] + out[i - neighbour_k])
                   + rbeta * out[i] - b[i]);
            sum += res * res;
        }

//#pragma acc wait
        res = sqrt(sum);
        it++;

    } //end while

    if (sync) {
//#pragma acc wait
    }

#ifndef BENCHMARKING
    m_logger->info("Number of iterations: {}", it);
    m_logger->info("Colored Gauss-Seidel ||res|| = {.5e}", res);
#endif
    }
}

//========================== Iteration step =================================
// *****************************************************************************
/// \brief  applies single CGS step on red-black grid
/// \param  out      output pointer
/// \param  b        source pointer
/// \param  alpha    2-dimensional array;
///                      \f$ (reciprocal_dx^2, reciprocal_dy^2)\f$ for pressure,
///                      \f$ (\nu\cdot dt\cdot reciprocal_dx^2, \nu\cdot dt\cdot reciprocal_dy^2)\f$  for velocity,
///                      \f$ (\kappa\cdot dt\cdot reciprocal_dx^2, \kappa\cdot dt\cdot reciprocal_dy^2)\f$ for temperature
/// \param  beta     \f$ 1./(2\cdot(\alpha_0 + \alpha_1) + 0.5)\f$ for pressure
///                  \f$ 1/(2\cdot(\alpha_0 + \alpha_1) + 1)\f$ for velocity
/// \param  dsign    sign (\a -1. for pressure, \a 1. else)
/// \param  w        weight (1. - diffusion, 2./3. - multigrid)
/// \param  sync     synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(
        Field &out, const Field &b,
        const real alpha_x, const real alpha_y, const real alpha_z,
        const real beta, const real dsign, const real w, bool) {

    auto domain = Domain::getInstance();
    // local parameters for GPU
    const size_t nx = domain->get_Nx();
    const size_t ny = domain->get_Ny();
    const size_t nz = domain->get_Nz();

    auto d_out = out.data;
    auto d_b = b.data;

//#pragma acc kernels present(d_out[:bsize], d_b[:bsize]) async
{
    // TODO: exclude obstacles!
    // red
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
} //end data region

//#pragma acc wait

//#pragma acc kernels present(d_out[:bsize], d_b[:bsize]) async
{
    //black
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
}

#pragma acc wait
}

// ============== Turbulent version of iteration step =================
// *****************************************************************************
/// \brief  applies single CGS step on red-black grid
/// \param  out      output pointer
/// \param  b        source pointer
/// \param  dsign    sign (\a -1. for pressure, \a 1. else)
/// \param  w        weight (1. - diffusion, 2./3. - multigrid)
/// \param  D        diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  EV       eddy viscosity (nu_turb)
/// \param  dt       time step
/// \param  sync     synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(
        Field &out, const Field &b,
        const real dsign, const real w, const real D,
        const Field &EV,
        const real dt, bool) {

    auto domain = Domain::getInstance();
    // local parameters for GPU
    const size_t Nx = domain->get_Nx(out.get_level());
    const size_t Ny = domain->get_Ny(out.get_level());
    const size_t Nz = domain->get_Nz(out.get_level());

    const real dx = domain->get_dx(out.get_level());  // due to unnecessary parameter passing of *this
    const real dy = domain->get_dy(out.get_level());
    const real dz = domain->get_dz(out.get_level());

    const real reciprocal_dx = 1. / dx;  // due to unnecessary parameter passing of *this
    const real reciprocal_dy = 1. / dy;
    const real reciprocal_dz = 1. / dz;

    real aX, aY, aZ, bb, rb;  // multipliers calculated

    auto d_out  = out.data;
    auto d_b    = b.data;
    auto d_EV   = EV.data;

// TODO: exclude obstacles!
//#pragma acc kernels present(d_out[:bsize], d_b[:bsize], d_EV[:bsize]) async
{
    // red
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < Nz - 1; k += 2) {
        for (size_t j = 1; j < Ny - 1; j += 2) {
            for (size_t i = 1; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dx * reciprocal_dx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dy * reciprocal_dy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dz * reciprocal_dz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
        for (size_t j = 2; j < Ny - 1; j += 2) {
            for (size_t i = 2; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dx * reciprocal_dx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dy * reciprocal_dy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dz * reciprocal_dz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < Nz - 1; k += 2) {
        for (size_t j = 1; j < Ny - 1; j += 2) {
            for (size_t i = 2; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dx * reciprocal_dx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dy * reciprocal_dy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dz * reciprocal_dz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
        for (size_t j = 2; j < Ny - 1; j += 2) {
            for (size_t i = 1; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dx * reciprocal_dx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dy * reciprocal_dy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dz * reciprocal_dz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
}
#pragma acc wait
//#pragma acc kernels present(d_out[:bsize], d_b[:bsize], d_EV[:bsize]) async
{
    // black
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < Nz - 1; k += 2) {
        for (size_t j = 1; j < Ny - 1; j += 2) {
            for (size_t i = 2; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dx * reciprocal_dx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dy * reciprocal_dy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dz * reciprocal_dz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < Nz - 1; k += 2) {
        for (size_t j = 2; j < Ny - 1; j += 2) {
            for (size_t i = 1; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dx * reciprocal_dx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dy * reciprocal_dy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dz * reciprocal_dz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < Nz - 1; k += 2) {
        for (size_t j = 1; j < Ny - 1; j += 2) {
            for (size_t i = 1; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dx * reciprocal_dx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dy * reciprocal_dy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dz * reciprocal_dz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < Nz - 1; k += 2) {
        for (size_t j = 2; j < Ny - 1; j += 2) {
            for (size_t i = 2; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dx * reciprocal_dx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dy * reciprocal_dy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * reciprocal_dz * reciprocal_dz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
} // end data region
//#pragma acc wait
}

// ========================= CGS stencil ==============================
// ************************************************************************
/// \brief  single CGS step
/// \param  i        index in x-direction
/// \param  j        index in y-direction
/// \param  k        index in z-direction
/// \param  out      output pointer
/// \param  b        source pointer
/// \param  alpha    2-dimensional array;
///                      \f$ (reciprocal_dx^2, reciprocal_dy^2)\f$ for pressure,
///                      \f$ (\nu\cdot dt\cdot reciprocal_dx^2, \nu\cdot dt\cdot reciprocal_dy^2)\f$  for velocity,
///                      \f$ (\kappa\cdot dt\cdot reciprocal_dx^2, \kappa\cdot dt\cdot reciprocal_dy^2)\f$ for temperature
/// \param  beta     \f$ 1./(2\cdot(\alpha_0 + \alpha_1) + 0.5)\f$ for pressure
///                  \f$ 1/(2\cdot(\alpha_0 + \alpha_1) + 1)\f$ for velocity
/// \param  dsign    sign (\a -1. for pressure, \a 1. else)
/// \param  w        weight (1. - diffusion, 2./3. - multigrid)
/// \param  Nx       number of cells in x-direction of computational domain
/// \param  Ny       number of cells in y-direction
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::colored_gauss_seidel_stencil(
        const size_t i, const size_t j, const size_t k,
        real *out, const real *b,
        const real alpha_x, const real alpha_y, const real alpha_z,
        const real dsign, const real beta, const real w,
        const size_t Nx, const size_t Ny) {
    real d_out_x    = *(out + IX(i + 1, j, k, Nx, Ny)); // per value (not access) necessary due to performance issues
    real d_out_x2   = *(out + (IX(i - 1, j, k, Nx, Ny)));
    real d_out_y    = *(out + (IX(i, j + 1, k, Nx, Ny)));
    real d_out_y2   = *(out + IX(i, j - 1, k, Nx, Ny));
    real d_out_z    = *(out + (IX(i, j, k + 1, Nx, Ny)));
    real d_out_z2   = *(out + (IX(i, j, k - 1, Nx, Ny)));
    real d_b        = *(b + IX(i, j, k, Nx, Ny));
    real r_out      = *(out + IX(i, j, k, Nx, Ny));

    real out_h      = beta * (dsign * d_b\
                         + alpha_x * (d_out_x + d_out_x2)\
                         + alpha_y * (d_out_y + d_out_y2)\
                         + alpha_z * (d_out_z + d_out_z2));

    *(out + IX(i, j, k, Nx, Ny)) = (1 - w) * r_out + w * out_h;
};
