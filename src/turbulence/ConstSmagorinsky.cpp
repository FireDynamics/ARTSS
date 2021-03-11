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
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"

ConstSmagorinsky::ConstSmagorinsky() {
    auto params = Parameters::getInstance();

    m_nu = params->get_real("physical_parameters/nu");
    m_dt = params->get_real("physical_parameters/dt");

    m_Cs = 0.1; //Cs value of 0.1 is found to yield the best results for wide range of flows
    // reference from Ansys Fluent Subgrid Scale models
    m_Cs = params->get_real("solver/turbulence/Cs");
}

//============================ Calculate turbulent viscosity =============================
// ***************************************************************************************
/// \brief  calculates turbulent viscosity
/// \param  ev            output pointer
/// \param  in_u          input pointer of x-velocity
/// \param  in_v          input pointer of y-velocity
/// \param  in_w          input pointer of z-velocity
/// \param  sync          synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ConstSmagorinsky::CalcTurbViscosity(Field &ev,
        Field const &in_u, Field const &in_v, Field const &in_w, bool sync) {
    auto domain = Domain::getInstance();

#pragma acc data present(ev, u, v, w)
    {
        const size_t Nx = domain->get_Nx(in_u.get_level());
        const size_t Ny = domain->get_Ny(in_v.get_level());

        const real dx = domain->get_dx(in_u.get_level());
        const real dy = domain->get_dy(in_v.get_level());
        const real dz = domain->get_dz(in_w.get_level());

        const real rdx = 1. / dx;
        const real rdy = 1. / dy;
        const real rdz = 1. / dz;

        const real delta_s = cbrt(dx * dy * dz);  // implicit filter
        real S_bar;
        real S11, S22, S33, S12, S13, S23;
        real Cs = m_Cs;

        auto boundary = BoundaryController::getInstance();
        size_t *d_iList = boundary->get_innerList_level_joined();
        auto bsize_i = boundary->getSize_innerList();

#pragma acc parallel loop independent present(ev, u, v, w, d_iList[:bsize_i]) async
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            S11 = (in_u[i + 1] - in_u[i - 1]) * 0.5 * rdx;
            S22 = (in_v[i + Nx] - in_v[i - Nx]) * 0.5 * rdy;
            S33 = (in_w[i + Nx * Ny] - in_w[i - Nx * Ny]) * 0.5 * rdz;
            S12 = 0.5 * ((in_u[i + Nx] - in_u[i - Nx]) * 0.5 * rdy \
 + (in_v[i + 1] - in_v[i - 1]) * 0.5 * rdx);
            S13 = 0.5 * ((in_u[i + Nx * Ny] - in_u[i - Nx * Ny]) * 0.5 * rdz  \
 + (in_w[i + 1] - in_w[i - 1]) * 0.5 * rdx);
            S23 = 0.5 * ((in_v[i + Nx * Ny] - in_v[i - Nx * Ny]) * 0.5 * rdz  \
 + (in_w[i + Nx] - in_w[i - Nx]) * 0.5 * rdy);

            S_bar = sqrt(2. * (S11 * S11 + S22 * S22 + S33 * S33 + 2. * (S12 * S12) + 2. * (S13 * S13) + 2. * (S23 * S23)));

            ev[i] = Cs * Cs * delta_s * delta_s * S_bar;
        }

        if (sync) {
#pragma acc wait
        }
    }
}

//============================ Explicit filtering =============================
// ***************************************************************************************
/// \brief  explicitly filters variables
/// \param  out           output pointer
/// \param  in            input pointer
/// \param  sync          synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ConstSmagorinsky::ExplicitFiltering(Field &, Field const &, bool) {
}
