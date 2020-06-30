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
void ConstSmagorinsky::CalcTurbViscosity(Field *ev, Field *in_u, Field *in_v, Field *in_w, bool sync) {

    auto domain = Domain::getInstance();
    // local parameters for GPU
    auto bsize = domain->get_size(in_u->GetLevel());

    auto d_u = in_u->data;
    auto d_v = in_v->data;
    auto d_w = in_w->data;
    auto d_ev = ev->data;

#pragma acc data present(d_ev[:bsize], d_u[:bsize], d_v[:bsize], d_w[:bsize])
    {
        const size_t Nx = domain->get_Nx(in_u->GetLevel());
        const size_t Ny = domain->get_Ny(in_v->GetLevel());

        const real dx = domain->get_dx(in_u->GetLevel());
        const real dy = domain->get_dy(in_v->GetLevel());
        const real dz = domain->get_dz(in_w->GetLevel());

        const real rdx = 1. / dx;
        const real rdy = 1. / dy;
        const real rdz = 1. / dz;

        const real delta_s = cbrt(dx * dy * dz); // implicit filter
        real S_bar;
        real S11, S22, S33, S12, S13, S23;
        real Cs = m_Cs;

        auto boundary = BoundaryController::getInstance();
        size_t *d_iList = boundary->get_innerList_level_joined();
        auto bsize_i = boundary->getSize_innerList();

#pragma acc parallel loop independent present(d_ev[:bsize], d_u[:bsize], d_v[:bsize], d_w[:bsize], d_iList[:bsize_i]) async
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            S11 = (d_u[i + 1] - d_u[i - 1]) * 0.5 * rdx;
            S22 = (d_v[i + Nx] - d_v[i - Nx]) * 0.5 * rdy;
            S33 = (d_w[i + Nx * Ny] - d_w[i - Nx * Ny]) * 0.5 * rdz;
            S12 = 0.5 * ((d_u[i + Nx] - d_u[i - Nx]) * 0.5 * rdy \
 + (d_v[i + 1] - d_v[i - 1]) * 0.5 * rdx);
            S13 = 0.5 * ((d_u[i + Nx * Ny] - d_u[i - Nx * Ny]) * 0.5 * rdz  \
 + (d_w[i + 1] - d_w[i - 1]) * 0.5 * rdx);
            S23 = 0.5 * ((d_v[i + Nx * Ny] - d_v[i - Nx * Ny]) * 0.5 * rdz  \
 + (d_w[i + Nx] - d_w[i - Nx]) * 0.5 * rdy);

            S_bar = sqrt(2. * (S11 * S11 + S22 * S22 + S33 * S33 + 2. * (S12 * S12) + 2. * (S13 * S13) + 2. * (S23 * S23)));

            d_ev[i] = Cs * Cs * delta_s * delta_s * S_bar;
        }

        if (sync) {
#pragma acc wait
        }
    } // end data
}

//============================ Explicit filtering =============================
// ***************************************************************************************
/// \brief  explicitly filters variables
/// \param  out           output pointer
/// \param  in            input pointer
/// \param  sync          synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ConstSmagorinsky::ExplicitFiltering(Field *out, const Field *in, bool sync) {
}
