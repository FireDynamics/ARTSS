/// \file         ConstSmagorinsky.cpp
/// \brief        calculates eddy viscosity based on Constant Smagorinsky-Lilly LES model
/// \date         Aug 18, 2016
/// \author       Suryanarayana Maddu
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>

#ifdef _OPENACC
#include <accelmath.h>
#endif

#include "ConstSmagorinsky.h"
#include "../DomainData.h"
#include "../boundary/DomainController.h"

ConstSmagorinsky::ConstSmagorinsky(Settings::Settings const &settings) {
    // Cs value of 0.1 is found to yield the best results for wide range of flows
    // reference from Ansys Fluent Subgrid Scale models
    m_Cs = settings.get_real("solver/turbulence/Cs");
}

//============================ Calculate turbulent viscosity =======================================
// *************************************************************************************************
/// \brief  calculates turbulent viscosity
/// \param  ev            output pointer
/// \param  in_u          input pointer of x-velocity
/// \param  in_v          input pointer of y-velocity
/// \param  in_w          input pointer of z-velocity
/// \param  sync          synchronization boolean (true=sync (default), false=async)
// *************************************************************************************************
void ConstSmagorinsky::calc_turbulent_viscosity(
            Field &ev,
            Field const &in_u, Field const &in_v, Field const &in_w,
            bool sync) {
    auto domain_data = DomainData::getInstance();
#pragma acc data present(ev, in_u, in_v, in_w)
    {
        const size_t Nx = domain_data->get_Nx();
        const size_t Ny = domain_data->get_Ny();

        const real dx = domain_data->get_dx();
        const real dy = domain_data->get_dy();
        const real dz = domain_data->get_dz();

        const real reciprocal_dx = 1. / dx;
        const real reciprocal_dy = 1. / dy;
        const real reciprocal_dz = 1. / dz;

        const real delta_s = cbrt(dx * dy * dz);  // implicit filter
        real S_bar;
        real S11, S22, S33, S12, S13, S23;
        real Cs = m_Cs;

        auto domain_controller = DomainController::getInstance();
        size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
        auto size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

        size_t neighbour_i = 1;
        size_t neighbour_j = Nx;
        size_t neighbour_k = Nx * Ny;
#pragma acc parallel loop independent present(ev, in_u, in_v, in_w, domain_inner_list[:size_domain_inner_list]) async
        for (size_t j = 0; j < size_domain_inner_list; ++j) {
            const size_t i = domain_inner_list[j];
            S11 = (in_u[i + neighbour_i] - in_u[i - neighbour_i]) * 0.5 * reciprocal_dx;
            S22 = (in_v[i + neighbour_j] - in_v[i - neighbour_j]) * 0.5 * reciprocal_dy;
            S33 = (in_w[i + neighbour_k] - in_w[i - neighbour_k]) * 0.5 * reciprocal_dz;
            S12 = 0.5 * ((in_u[i + neighbour_j] - in_u[i - neighbour_j]) * 0.5 * reciprocal_dy
                       + (in_v[i + neighbour_i] - in_v[i - neighbour_i]) * 0.5 * reciprocal_dx);
            S13 = 0.5 * ((in_u[i + neighbour_k] - in_u[i - neighbour_k]) * 0.5 * reciprocal_dz
                       + (in_w[i + neighbour_i] - in_w[i - neighbour_i]) * 0.5 * reciprocal_dx);
            S23 = 0.5 * ((in_v[i + neighbour_k] - in_v[i - neighbour_k]) * 0.5 * reciprocal_dz
                       + (in_w[i + neighbour_j] - in_w[i - neighbour_j]) * 0.5 * reciprocal_dy);

            S_bar = sqrt(2. * (S11 * S11 + S22
                            *  S22 + S33 * S33
                       + 2. * (S12 * S12)
                       + 2. * (S13 * S13)
                       + 2. * (S23 * S23)));

            ev[i] = Cs * Cs * delta_s * delta_s * S_bar;
        }

        if (sync) {
#pragma acc wait
        }
        domain_controller->apply_boundary(ev, sync);
    }
}

//============================ Explicit filtering ==============================
// *****************************************************************************
/// \brief  explicitly filters variables
/// \param  out           output pointer
/// \param  in            input pointer
/// \param  sync          synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void ConstSmagorinsky::explicit_filtering(Field &, Field const &, bool) {
}
