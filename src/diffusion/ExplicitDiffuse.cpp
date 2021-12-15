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
#include "../boundary/BoundaryController.h"
#include "../DomainData.h"

//====================================== Diffuse ===============================================
// ***************************************************************************************
/// \brief  solves diffusion equation \f$ \partial_t \phi_2 = \nu \ nabla^2 \phi_2 \f$
/// \param  out     output pointer
/// \param  in      input pointer
/// \param  b       source pointer
/// \param  D       diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ExplicitDiffuse::diffuse(Field &out, const Field &in, Field const &, real const D, bool sync) {
    ExplicitStep(out, in, D, sync);
    BoundaryController::getInstance()->apply_boundary(out, sync);
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
        Field &out, Field const &in,
        Field const &, real const D, Field const &EV, bool sync) {
    ExplicitStep(out, in, D, EV, sync);
    BoundaryController::getInstance()->apply_boundary(out, sync);
}


void ExplicitDiffuse::ExplicitStep(Field &out, Field const &in, real const D, bool sync) {
    auto domain = DomainData::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx();  // due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny();

    const real dt = m_settings.get_real("physical_parameters/dt");

    auto boundary = BoundaryController::getInstance();

    size_t *d_inner_list = boundary->get_domain_inner_list_level_joined();
    auto bsize_i = boundary->get_size_domain_inner_list_level_joined(0);

    real reciprocal_dx = D / (domain->get_dx() * domain->get_dx());
    real reciprocal_dy = D / (domain->get_dy() * domain->get_dy());
    real reciprocal_dz = D / (domain->get_dz() * domain->get_dz());

    const size_t neighbour_i = 1;
    const size_t neighbour_j = Nx;
    const size_t neighbour_k = Nx * Ny;
#pragma acc parallel loop independent present(out, in, d_inner_list[:bsize_i]) async
    for (size_t ii = 0; ii < bsize_i; ++ii) {
        const size_t idx = d_inner_list[ii];
        out[idx] = in[idx] +
                     dt * ((in[idx + neighbour_i] - 2 * in[idx] + in[idx - neighbour_i]) * reciprocal_dx
                         + (in[idx + neighbour_j] - 2 * in[idx] + in[idx - neighbour_j]) * reciprocal_dy
                         + (in[idx + neighbour_k] - 2 * in[idx] + in[idx - neighbour_k]) * reciprocal_dz
                     );
    }

    if (sync) {
#pragma acc wait
    }
}

// Turbulent Diffuse
void ExplicitDiffuse::ExplicitStep(Field &out, const Field &in, real const D, Field const &EV, bool sync) {
    auto domain = DomainData::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx();  // due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny();

    const real dt = m_settings.get_real("physical_parameters/dt");

    auto boundary = BoundaryController::getInstance();

    size_t *d_inner_list = boundary->get_domain_inner_list_level_joined();
    auto bsize_i = boundary->get_size_domain_inner_list_level_joined(0);

    const size_t neighbour_i = 1;
    const size_t neighbour_j = Nx;
    const size_t neighbour_k = Nx * Ny;
#pragma acc parallel loop independent present(out, in, d_inner_list[:bsize_i], EV) async
    for (size_t ii = 0; ii < bsize_i; ++ii) {
        const size_t idx = d_inner_list[ii];
        real nu_x1 = 0.5 * (EV[idx + neighbour_i] + EV[idx]) + D;
        real nu_x2 = 0.5 * (EV[idx - neighbour_i] + EV[idx]) + D;
        real nu_y1 = 0.5 * (EV[idx + neighbour_j] + EV[idx]) + D;
        real nu_y2 = 0.5 * (EV[idx - neighbour_j] + EV[idx]) + D;
        real nu_z1 = 0.5 * (EV[idx + neighbour_k] + EV[idx]) + D;
        real nu_z2 = 0.5 * (EV[idx - neighbour_k] + EV[idx]) + D;

        real di0 = (in[idx] - in[idx - neighbour_i]);  // u_i - u_{i-1}
        real di1 = (in[idx + neighbour_i] - in[idx]);  // u_{i+1} - u_i
        real dj0 = (in[idx] - in[idx - neighbour_j]);
        real dj1 = (in[idx + neighbour_j] - in[idx]);
        real dk0 = (in[idx] - in[idx - neighbour_k]);
        real dk1 = (in[idx + neighbour_k] - in[idx]);

        out[idx] = in[idx] +
                     dt * ((nu_x1 * di0 - nu_x2 * di1)    // dx // * 4 / dx;  // u_{i-0.5} - u_{i+0.5}
                         + (nu_y1 * dj0 - nu_y2 * dj1)    // dy
                         + (nu_z1 * dk0 - nu_z2 * dk1));  // dz
    }

    if (sync) {
#pragma acc wait
    }
}
