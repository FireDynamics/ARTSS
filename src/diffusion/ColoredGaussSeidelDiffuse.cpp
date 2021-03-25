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
void ColoredGaussSeidelDiffuse::diffuse(Field &out, Field &, Field const &b, const real D, bool sync) {
    auto domain = Domain::getInstance();
    // local parameters for GPU
    FieldType type = out.getType();

    auto d_out = out.data;
    auto d_b = b.data;

    auto boundary = BoundaryController::getInstance();

    auto bsize_i = boundary->getSize_innerList();
    size_t* d_iList = boundary->get_innerList_level_joined();

{
    const size_t Nx = domain->get_Nx(out.getLevel());  // due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny(out.getLevel());

    const real dx = domain->get_dx(out.getLevel());  // due to unnecessary parameter passing of *this
    const real dy = domain->get_dy(out.getLevel());
    const real dz = domain->get_dz(out.getLevel());

    const real rdx = 1. / dx;
    const real rdy = 1. / dy;
    const real rdz = 1. / dz;

    const real alpha_x = D * m_dt * rdx * rdx;  // due to better pgi handling of scalars (instead of arrays)
    const real alpha_y = D * m_dt * rdy * rdy;
    const real alpha_z = D * m_dt * rdz * rdz;

    const real rbeta    = (1. + 2. * (alpha_x + alpha_y + alpha_z));
    const real beta     = 1. / rbeta;

    const real dsign    = m_dsign;
    const real w        = m_w;

    size_t it           = 0;
    const size_t max_it = m_max_iter;
    const real tol_res  = m_tol_res;

    real sum;
    real res = 1.;

    while (res > tol_res && it < max_it) {
        colored_gauss_seidel_step(out, b, alpha_x, alpha_y, alpha_z, beta, dsign, w, sync);
        boundary->applyBoundary(d_out, type, sync);

        sum = 0;

        for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
            res = (- alpha_x * (d_out[i + 1] + d_out[i - 1])
                   - alpha_y * (d_out[i + Nx] + d_out[i - Nx])
                   - alpha_z * (d_out[i + Nx * Ny] + d_out[i - Nx * Ny])
                   + rbeta * d_out[i] - d_b[i]);
            // TODO: find a way to exclude obstacle indices! Jacobi now uses res=D*(out-in)
            sum += res * res;
        }

        res = sqrt(sum);
        it++;
    }

#ifndef BENCHMARKING
    m_logger->info("Number of iterations: {}", it);
    m_logger->info("Colored Gauss-Seidel ||res|| = {:.5e}", res);
#endif
}
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
void ColoredGaussSeidelDiffuse::diffuse(Field &out, Field &, Field const &b,
        real const D, Field const &EV, bool sync) {
    auto domain = Domain::getInstance();
    // local parameters for GPU
    FieldType type = out.getType();

    auto d_out  = out.data;
    auto d_b    = b.data;
    auto d_EV   = EV.data;

    auto boundary = BoundaryController::getInstance();

    auto bsize_i = boundary->getSize_innerList();
    size_t* d_iList = boundary->get_innerList_level_joined();

{
    const size_t Nx = domain->get_Nx(out.getLevel());
    const size_t Ny = domain->get_Ny(out.getLevel());

    const real dx = domain->get_dx(out.getLevel());
    const real dy = domain->get_dy(out.getLevel());
    const real dz = domain->get_dz(out.getLevel());

    const real rdx = 1. / dx;
    const real rdy = 1. / dy;
    const real rdz = 1. / dz;

    real dt = m_dt;

    real alpha_x, alpha_y, alpha_z, rbeta;  // calculated in colored_gauss_seidel_step!

    const real dsign = m_dsign;
    const real w = m_w;

    size_t it = 0;
    const size_t max_it = m_max_iter;
    const real tol_res = m_tol_res;
    real sum;
    real res = 1.;

    while (res > tol_res && it < max_it) {
        colored_gauss_seidel_step(out, b, dsign, w, D, EV, dt, sync);
        boundary->applyBoundary(d_out, type, sync);

        sum = 0;

        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            alpha_x = (D + d_EV[i]) * dt * rdx * rdx;
            alpha_y = (D + d_EV[i]) * dt * rdy * rdy;
            alpha_z = (D + d_EV[i]) * dt * rdz * rdz;
            rbeta = (1. + 2. * (alpha_x + alpha_y + alpha_z));

            res = (- alpha_x * (d_out[i + 1      ] + d_out[i - 1      ])
                   - alpha_y * (d_out[i + Nx     ] + d_out[i - Nx     ])
                   - alpha_z * (d_out[i + Nx * Ny] + d_out[i - Nx * Ny])
                   + rbeta * d_out[i] - d_b[i]);
            sum += res * res;
        }

        res = sqrt(sum);
        it++;
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
///                      \f$ (rdx^2, rdy^2)\f$ for pressure,
///                      \f$ (\nu\cdot dt\cdot rdx^2, \nu\cdot dt\cdot rdy^2)\f$  for velocity,
///                      \f$ (\kappa\cdot dt\cdot rdx^2, \kappa\cdot dt\cdot rdy^2)\f$ for temperature
/// \param  beta     \f$ 1./(2\cdot(\alpha_0 + \alpha_1) + 0.5)\f$ for pressure
///                  \f$ 1/(2\cdot(\alpha_0 + \alpha_1) + 1)\f$ for velocity
/// \param  dsign    sign (\a -1. for pressure, \a 1. else)
/// \param  w        weight (1. - diffusion, 2./3. - multigrid)
/// \param  sync     synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(
        Field &out, Field const &b,
        real const alpha_x, real const alpha_y, real const alpha_z,
        real const beta, real const dsign, real const w, bool) {
    auto domain = Domain::getInstance();
    // local parameters for GPU
    const size_t nx = domain->get_Nx(out.getLevel());
    const size_t ny = domain->get_Ny(out.getLevel());
    const size_t nz = domain->get_Nz(out.getLevel());

    auto d_out = out.data;
    auto d_b = b.data;

{
    // TODO: exclude obstacles!
    // red
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
}

{
    //black
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, nx, ny);
            }
        }
    }
}
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
        Field &out, Field const &b,
        real const dsign, real const w, real const D,
        Field const &EV, real const dt, bool) {
    auto domain = Domain::getInstance();
    // local parameters for GPU
    const size_t Nx = domain->get_Nx(out.getLevel());
    const size_t Ny = domain->get_Ny(out.getLevel());
    const size_t Nz = domain->get_Nz(out.getLevel());

    const real dx = domain->get_dx(out.getLevel());  // due to unnecessary parameter passing of *this
    const real dy = domain->get_dy(out.getLevel());
    const real dz = domain->get_dz(out.getLevel());

    const real rdx = 1. / dx;  // due to unnecessary parameter passing of *this
    const real rdy = 1. / dy;
    const real rdz = 1. / dz;

    real aX, aY, aZ, bb, rb;  // multipliers calculated

    auto d_out  = out.data;
    auto d_b    = b.data;
    auto d_EV   = EV.data;

// TODO: exclude obstacles!
{
    // red
    for (size_t k = 1; k < Nz - 1; k += 2) {
        for (size_t j = 1; j < Ny - 1; j += 2) {
            for (size_t i = 1; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
        for (size_t j = 2; j < Ny - 1; j += 2) {
            for (size_t i = 2; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
    for (size_t k = 2; k < Nz - 1; k += 2) {
        for (size_t j = 1; j < Ny - 1; j += 2) {
            for (size_t i = 2; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
        for (size_t j = 2; j < Ny - 1; j += 2) {
            for (size_t i = 1; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
}
{
    // black
    for (size_t k = 1; k < Nz - 1; k += 2) {
        for (size_t j = 1; j < Ny - 1; j += 2) {
            for (size_t i = 2; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
    for (size_t k = 1; k < Nz - 1; k += 2) {
        for (size_t j = 2; j < Ny - 1; j += 2) {
            for (size_t i = 1; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
    for (size_t k = 2; k < Nz - 1; k += 2) {
        for (size_t j = 1; j < Ny - 1; j += 2) {
            for (size_t i = 1; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
    for (size_t k = 2; k < Nz - 1; k += 2) {
        for (size_t j = 2; j < Ny - 1; j += 2) {
            for (size_t i = 2; i < Nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, Nx, Ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                colored_gauss_seidel_stencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
            }
        }
    }
}
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
///                      \f$ (rdx^2, rdy^2)\f$ for pressure,
///                      \f$ (\nu\cdot dt\cdot rdx^2, \nu\cdot dt\cdot rdy^2)\f$  for velocity,
///                      \f$ (\kappa\cdot dt\cdot rdx^2, \kappa\cdot dt\cdot rdy^2)\f$ for temperature
/// \param  beta     \f$ 1./(2\cdot(\alpha_0 + \alpha_1) + 0.5)\f$ for pressure
///                  \f$ 1/(2\cdot(\alpha_0 + \alpha_1) + 1)\f$ for velocity
/// \param  dsign    sign (\a -1. for pressure, \a 1. else)
/// \param  w        weight (1. - diffusion, 2./3. - multigrid)
/// \param  Nx       number of cells in x-direction of computational domain
/// \param  Ny       number of cells in y-direction
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::colored_gauss_seidel_stencil(
        size_t i, size_t j, size_t k,
        real *out, real *b,
        real const alpha_x, real const alpha_y, real const alpha_z,
        real const dsign, real const beta, real const w,
        size_t const Nx, size_t const Ny) {
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
