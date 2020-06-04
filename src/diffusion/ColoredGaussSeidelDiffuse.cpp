/// \file         ColoredGaussSeidelDiffuse.cpp
/// \brief        solves Diffusion equation \f$ \partial_t \phi_2 = \nu \nabla^2 \phi_2 \f$
///               via calculated iterations of CGS step (dependent on residual/ maximal 1000)
/// \date         June 21, 2016
/// \author       lgewuerz
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>
#ifndef PROFILING
#include <spdlog/spdlog.h>
#endif

#include "ColoredGaussSeidelDiffuse.h"
#include "../boundary/BoundaryController.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../utility/Utility.h"

//=============================== Constructor ======================================
ColoredGaussSeidelDiffuse::ColoredGaussSeidelDiffuse() {
#ifndef PROFILING
    m_logger = Utility::createLogger(typeid(this).name());
#endif
    auto params = Parameters::getInstance();

    m_dt = params->getReal("physical_parameters/dt");
    m_dsign = 1.;
    m_w = params->getReal("solver/diffusion/w");

    if (params->get("solver/diffusion/type") == "ColoredGaussSeidel") {
        m_max_iter = static_cast<size_t>(params->getInt("solver/diffusion/max_iter"));
        m_tol_res = params->getReal("solver/diffusion/tol_res");
    }
    else {
        m_max_iter = 10000;
        m_tol_res = 1e-16;
    }
}

//==================================== Diffuse ===========================================
// ***************************************************************************************
/// \brief  solves Diffusion equation \f$ \partial_t \phi_2 = \nu \ nabla^2 \phi_2 \f$
/// 		via calculated iterations of CGS step (dependent on residual/ maximal iterations)
/// \param  out			output pointer
/// \param	in			input pointer
/// \param	b 			source pointer
/// \param	D			diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  sync		synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::diffuse(Field *out, Field *in, const Field *b, const real D, bool sync) {

    auto domain = Domain::getInstance();
    // local parameters for GPU
    auto bsize = domain->GetSize(out->GetLevel());
    FieldType type = out->GetType();

    auto d_out = out->data;
    auto d_b = b->data;

    auto boundary = BoundaryController::getInstance();

	auto bsize_i = boundary->getSize_innerList();
	auto bsize_b = boundary->getSize_boundaryList();

    size_t* d_iList = boundary->get_innerList_level_joined();
    size_t* d_bList = boundary->get_boundaryList_level_joined();

//#pragma acc data present(d_out[:bsize], d_b[:bsize])
{
    const size_t Nx = domain->GetNx(out->GetLevel()); //due to unnecessary parameter passing of *this
    const size_t Ny = domain->GetNy(out->GetLevel());
    const size_t Nz = domain->GetNz(out->GetLevel());

    const real dx = domain->Getdx(out->GetLevel()); //due to unnecessary parameter passing of *this
    const real dy = domain->Getdy(out->GetLevel());
    const real dz = domain->Getdz(out->GetLevel());

    const real rdx = 1. / dx;
    const real rdy = 1. / dy;
    const real rdz = 1. / dz;

    const real alphaX = D * m_dt * rdx * rdx; //due to better pgi handling of scalars (instead of arrays)
    const real alphaY = D * m_dt * rdy * rdy;
    const real alphaZ = D * m_dt * rdz * rdz;

    const real rbeta    = (1. + 2. * (alphaX + alphaY + alphaZ));
    const real beta     = 1. / rbeta;

    const real dsign    = m_dsign;
    const real w        = m_w;

    size_t it           = 0;
    const size_t max_it = m_max_iter;
    const real tol_res  = m_tol_res;

    real sum;
    real res = 1.;

    while (res > tol_res && it < max_it) {
        ColoredGaussSeidelStep(out, b, alphaX, alphaY, alphaZ, beta, dsign, w, sync);
        boundary->applyBoundary(d_out, type, sync);

        sum = 0;

//#pragma acc parallel loop independent present(d_out[:bsize], d_b[:bsize], d_iList[:bsize_i]) async
        for (size_t j = 0; j < bsize_i; ++j){
        const size_t i = d_iList[j];
            res = (	- alphaX * (d_out[i + 1] + d_out[i - 1])\
                    - alphaY * (d_out[i + Nx] + d_out[i - Nx])\
                    - alphaZ * (d_out[i + Nx * Ny] + d_out[i - Nx * Ny])\
                    + rbeta * d_out[i] - d_b[i]);
            //TODO: find a way to exclude obstacle indices! Jacobi now uses res=D*(out-in)
            sum += res * res;
        }

//#pragma acc wait
        res = sqrt(sum);
        it++;
    } //end while

    if (sync)
{
//#pragma acc wait
    }

#ifndef PROFILING
    m_logger->info("Number of iterations: {}", it);
    m_logger->info("Colored Gauss-Seidel ||res|| = {:.5e}", res);
#endif
} //end data region
}

//=============================== Turbulent version ======================================
// ***************************************************************************************
/// \brief  solves Diffusion equation \f$ \partial_t \phi_2 = \nu_{eff} \ nabla^2 \phi_2 \f$
///         via calculated iterations of CGS step (dependent on residual/ maximal iterations)
/// \param  out         output pointer
/// \param  in          input pointer
/// \param  b           source pointer
/// \param  D           diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  EV          eddy viscosity (nu_turb)
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::diffuse(Field *out, Field *in, const Field *b, const real D, const Field *EV, bool sync) {

    auto domain = Domain::getInstance();
    // local parameters for GPU
    auto bsize = domain->GetSize(out->GetLevel());
    FieldType type = out->GetType();

    auto d_out  = out->data;
    auto d_b    = b->data;
    auto d_EV   = EV->data;

    auto boundary = BoundaryController::getInstance();

    auto bsize_i = boundary->getSize_innerList();
    auto bsize_b = boundary->getSize_boundaryList();

    size_t* d_iList = boundary->get_innerList_level_joined();
    size_t* d_bList = boundary->get_boundaryList_level_joined();

//#pragma acc data present(d_out[:bsize], d_b[:bsize], d_EV[:bsize])
{
    const size_t Nx = domain->GetNx(out->GetLevel()); //due to unnecessary parameter passing of *this
    const size_t Ny = domain->GetNy(out->GetLevel());
    const size_t Nz = domain->GetNz(out->GetLevel());

    const real dx = domain->Getdx(out->GetLevel()); //due to unnecessary parameter passing of *this
    const real dy = domain->Getdy(out->GetLevel());
    const real dz = domain->Getdz(out->GetLevel());

    const real rdx = 1. / dx;
    const real rdy = 1. / dy;
    const real rdz = 1. / dz;

    real dt = m_dt;

    real alphaX, alphaY, alphaZ, rbeta; //calculated in ColoredGaussSeidelStep!

    const real dsign = m_dsign;
    const real w = m_w;

    size_t it = 0;
    const size_t max_it = m_max_iter;
    const real tol_res = m_tol_res;
    real sum;
    real res = 1.;

    while (res > tol_res && it < max_it) {
        ColoredGaussSeidelStep(out, b, dsign, w, D, EV, dt, sync);
        boundary->applyBoundary(d_out, type, sync);

        sum = 0;

//#pragma acc parallel loop independent present(d_out[:bsize], d_b[:bsize], d_EV[:bsize], d_iList[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j){
        const size_t i = d_iList[j];
        alphaX = (D + d_EV[i]) * dt * rdx * rdx;
        alphaY = (D + d_EV[i]) * dt * rdy * rdy;
        alphaZ = (D + d_EV[i]) * dt * rdz * rdz;
        rbeta = (1. + 2. * (alphaX + alphaY + alphaZ));

        res = ( - alphaX * (d_out[i + 1      ] + d_out[i - 1      ])\
                - alphaY * (d_out[i + Nx     ] + d_out[i - Nx     ])\
                - alphaZ * (d_out[i + Nx * Ny] + d_out[i - Nx * Ny])\
                + rbeta * d_out[i] - d_b[i]);
        sum += res * res;
    }

//#pragma acc wait
        res = sqrt(sum);
        it++;

    } //end while

    if (sync)
{
//#pragma acc wait
    }

#ifndef PROFILING
    m_logger->info("Number of iterations: {}", it);
    m_logger->info("Colored Gauss-Seidel ||res|| = {.5e}", res);
#endif
} //end data region
}

//=============================== Iteration step ======================================
// ***************************************************************************************
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
void ColoredGaussSeidelDiffuse::ColoredGaussSeidelStep(Field *out, const Field *b, const real alphaX, const real alphaY, const real alphaZ, const real beta, const real dsign, const real w, bool sync) {

    auto domain = Domain::getInstance();
    // local parameters for GPU
    const size_t nx = domain->GetNx(out->GetLevel());
    const size_t ny = domain->GetNy(out->GetLevel());
    const size_t nz = domain->GetNz(out->GetLevel());

    const size_t bsize = domain->GetSize(out->GetLevel());

    auto d_out = out->data;
    auto d_b = b->data;

//#pragma acc kernels present(d_out[:bsize], d_b[:bsize]) async
{
//TODO: exclude obstacles!
    //red
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, alphaX, alphaY, alphaZ, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, alphaX, alphaY, alphaZ, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, alphaX, alphaY, alphaZ, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, alphaX, alphaY, alphaZ, dsign, beta, w, nx, ny);
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
                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, alphaX, alphaY, alphaZ, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, alphaX, alphaY, alphaZ, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, alphaX, alphaY, alphaZ, dsign, beta, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, alphaX, alphaY, alphaZ, dsign, beta, w, nx, ny);
            }
        }
    }
} //end data region

//#pragma acc wait
}

//======================== Turbulent version of iteration step ===========================
// ***************************************************************************************
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
void ColoredGaussSeidelDiffuse::ColoredGaussSeidelStep(Field *out, const Field *b, const real dsign, const real w, const real D, const Field *EV, const real dt, bool sync) {

    auto domain = Domain::getInstance();
    // local parameters for GPU
    const size_t nx = domain->GetNx(out->GetLevel());
    const size_t ny = domain->GetNy(out->GetLevel());
    const size_t nz = domain->GetNz(out->GetLevel());

    const real dx = domain->Getdx(out->GetLevel()); //due to unnecessary parameter passing of *this
    const real dy = domain->Getdy(out->GetLevel());
    const real dz = domain->Getdz(out->GetLevel());

    const real rdx = 1. / dx; //due to unnecessary parameter passing of *this
    const real rdy = 1. / dy;
    const real rdz = 1. / dz;

    size_t bsize = domain->GetSize(out->GetLevel());

    real aX, aY, aZ, bb, rb; //multipliers calculated

    auto d_out  = out->data;
    auto d_b    = b->data;
    auto d_EV   = EV->data;

//TODO: exclude obstacles!
//#pragma acc kernels present(d_out[:bsize], d_b[:bsize], d_EV[:bsize]) async
{

    //red
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, nx, ny);
            }
        }
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, nx, ny);
            }
        }
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, nx, ny);
            }
        }
    }
} //end data region
//#pragma acc wait
//#pragma acc kernels present(d_out[:bsize], d_b[:bsize], d_EV[:bsize]) async
{
    //black
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 1; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 1; j < ny - 1; j += 2) {
            for (size_t i = 1; i < nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, nx, ny);
            }
        }
    }
//#pragma acc loop independent collapse(3)
    for (size_t k = 2; k < nz - 1; k += 2) {
        for (size_t j = 2; j < ny - 1; j += 2) {
            for (size_t i = 2; i < nx - 1; i += 2) {
                aX = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdx * rdx;
                aY = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdy * rdy;
                aZ = (D + d_EV[IX(i, j, k, nx, ny)]) * dt * rdz * rdz;

                rb = (1. + 2. * (aX + aY + aZ));
                bb = 1. / rb;

                ColoredGaussSeidelStencil(i, j, k, d_out, d_b, aX, aY, aZ, dsign, bb, w, nx, ny);
            }
        }
    }
} //end data region
//#pragma acc wait
}

//=================================== CGS stencil ========================================
// ***************************************************************************************
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
/// \param  nx       number of cells in x-direction of computational domain
/// \param  ny       number of cells in y-direction
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::ColoredGaussSeidelStencil(size_t i, size_t j, size_t k, real *out, real *b, const real alphaX, const real alphaY, const real alphaZ, const real dsign, const real beta, const real w, const size_t nx, const size_t ny) {

    real d_out_x    = *(out + IX(i + 1, j, k, nx, ny)); // per value (not access) necessary due to performance issues
    real d_out_x2   = *(out + (IX(i - 1, j, k, nx, ny)));
    real d_out_y    = *(out + (IX(i, j + 1, k, nx, ny)));
    real d_out_y2   = *(out + IX(i, j - 1, k, nx, ny));
    real d_out_z    = *(out + (IX(i, j, k + 1, nx, ny)));
    real d_out_z2   = *(out + (IX(i, j, k - 1, nx, ny)));
    real d_b        = *(b + IX(i, j, k, nx, ny));
    real r_out      = *(out + IX(i, j, k, nx, ny));

    real out_h      = beta * (dsign * d_b\
                         + alphaX * (d_out_x + d_out_x2)\
                         + alphaY * (d_out_y + d_out_y2)\
                         + alphaZ * (d_out_z + d_out_z2));

    *(out + IX(i, j, k, nx, ny)) = (1 - w) * r_out + w * out_h;
}
