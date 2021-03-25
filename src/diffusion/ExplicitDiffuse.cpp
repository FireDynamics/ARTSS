/// \file       ExplicitDiffuse.cpp
/// \brief      Solves diffusion equation with an explicit method
/// \date       December 12, 2019
/// \author     Max Boehler
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>

#ifdef _OPENACC
#include <accelmath.h>
#endif

#include "ExplicitDiffuse.h"
#include "../utility/Parameters.h"
#include "../boundary/BoundaryController.h"
#include "../Domain.h"

ExplicitDiffuse::ExplicitDiffuse() {
    auto params = Parameters::getInstance();

    m_dt = params->get_real("physical_parameters/dt");
}

//====================================== Diffuse ===============================================
// ***************************************************************************************
/// \brief  solves diffusion equation \f$ \partial_t \phi_2 = \nu \ nabla^2 \phi_2 \f$
/// \param  out     output pointer
/// \param  in      input pointer
/// \param  b       source pointer
/// \param  D       diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ExplicitDiffuse::diffuse(Field &out, Field &in, Field const &, real const D, bool sync) {
    // local variables and parameters for GPU
    FieldType type = out.get_type();

    auto d_out = out.data;

    auto boundary = BoundaryController::getInstance();

    ExplicitStep(out, in, D, sync);
    boundary->applyBoundary(d_out, type, sync);
}

//====================================== Turbulent Diffuse ===============================================
// ***************************************************************************************
/// \brief  solves diffusion equation \f$ \partial_t \phi_2 = \nu \ nabla^2 \phi_2 \f$
/// \param  out     output pointer
/// \param  in      input pointer
/// \param  b       source pointer
/// \param  D       diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ExplicitDiffuse::diffuse(
        Field &out, Field &in,
        Field const &, real const D, Field const &EV, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    auto bsize = domain->get_size(out.get_level());
    FieldType type = out.get_type();

    auto d_out = out.data;

    auto boundary = BoundaryController::getInstance();

    ExplicitStep(out, in, D, EV, sync);
    boundary->applyBoundary(d_out, type, sync);
}


void ExplicitDiffuse::ExplicitStep(Field &out, Field const &in, real const D, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(out.get_level()); //due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny(out.get_level());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

    real rdx = D / (domain->get_dx(out.get_level()) * domain->get_dx(out.get_level()));
    real rdy = D / (domain->get_dy(out.get_level()) * domain->get_dy(out.get_level()));
    real rdz = D / (domain->get_dz(out.get_level()) * domain->get_dz(out.get_level()));

#pragma acc parallel loop independent present(out, in, d_iList[:bsize_i]) async
    for (size_t ii = 0; ii < bsize_i; ++ii) {
        const size_t idx = d_iList[ii];

        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);

        out[idx] = in[idx] +
                     m_dt * ((in[IX(i + 1, j, k, Nx, Ny)] - 2 * in[idx] + in[IX(i - 1, j, k, Nx, Ny)]) * rdx
                           + (in[IX(i, j + 1, k, Nx, Ny)] - 2 * in[idx] + in[IX(i, j - 1, k, Nx, Ny)]) * rdy
                           + (in[IX(i, j, k + 1, Nx, Ny)] - 2 * in[idx] + in[IX(i, j, k - 1, Nx, Ny)]) * rdz
                     );
    }

    if (sync) {
#pragma acc wait
    }
}


// Turbulent Diffuse
void ExplicitDiffuse::ExplicitStep(Field &out, const Field &in, real const D, Field const &EV, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(out.get_level());  // due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny(out.get_level());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

#pragma acc parallel loop independent present(out, in, d_iList[:bsize_i], ev) async
    for (size_t ii = 0; ii < bsize_i; ++ii) {
        const size_t idx = d_iList[ii];

        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);

        real nu_x1 = 0.5 * (EV[IX(i - 1, j, k, Nx, Ny)] + EV[idx]) + D;
        real nu_x2 = 0.5 * (EV[IX(i + 1, j, k, Nx, Ny)] + EV[idx]) + D;
        real nu_y1 = 0.5 * (EV[IX(i, j - 1, k, Nx, Ny)] + EV[idx]) + D;
        real nu_y2 = 0.5 * (EV[IX(i, j + 1, k, Nx, Ny)] + EV[idx]) + D;
        real nu_z1 = 0.5 * (EV[IX(i, j, k - 1, Nx, Ny)] + EV[idx]) + D;
        real nu_z2 = 0.5 * (EV[IX(i, j, k + 1, Nx, Ny)] + EV[idx]) + D;

        auto di0 = (in[idx] - in[IX(i - 1, j, k, Nx, Ny)]);  // u_i - u_{i-1}
        auto di1 = (in[IX(i + 1, j, k, Nx, Ny)] - in[idx]);  // u_{i+1} - u_i
        auto dj0 = (in[idx] - in[IX(i, j - 1, k, Nx, Ny)]);
        auto dj1 = (in[IX(i, j + 1, k, Nx, Ny)] - in[idx]);
        auto dk0 = (in[idx] - in[IX(i, j, k - 1, Nx, Ny)]);
        auto dk1 = (in[IX(i, j, k + 1, Nx, Ny)] - in[idx]);

        out[idx] = in[idx] +
                     m_dt * ((nu_x1 * di0 - nu_x2 * di1)   // dx // * 4 / dx;  // u_{i-0.5} - u_{i+0.5}
                           + (nu_y1 * dj0 - nu_y2 * dj1)   // dy
                           + (nu_z1 * dk0 - nu_z2 * dk1)); // dz
    }

    if (sync) {
#pragma acc wait
    }
}
