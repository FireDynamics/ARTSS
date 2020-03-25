/// \file 		ExplicitAdvect.cpp
/// \brief 		Explicit (BD) solver for advection equation
/// \date 		Aug 22, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "ExplicitAdvect.h"
#include "../utility/Parameters.h"
#include "../boundary/BoundaryController.h"
#include "../Domain.h"

ExplicitAdvect::ExplicitAdvect() {

    auto params = Parameters::getInstance();

    m_dt = params->getReal("physical_parameters/dt");
}

// ==================================== Advect ====================================
// ***************************************************************************************
/// \brief  solves advection \f$ \partial_t \phi_1 = - (u \cdot \nabla) \phi_0 \f$ via
///			explicit finite differences (upwind - FD in time, BD in space) with CFL \f$c < \Delta x/ \Delta t\f$
/// \param  out		output pointer
/// \param	in		input pointer
/// \param	u_vel	x -velocity
/// \param	v_vel	y -velocity
/// \param	w_vel	z -velocity
/// \param	sync	synchronous kernel launching (true, default: false)
// ***************************************************************************************
void ExplicitAdvect::advect(Field *out, Field *in, const Field *u_vel, const Field *v_vel, const Field *w_vel, bool sync) {

    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    size_t bsize = domain->GetSize(out->GetLevel());
    FieldType type = out->GetType();

    auto d_out = out->data;
    auto d_in = in->data;
    auto d_u_vel = u_vel->data;
    auto d_v_vel = v_vel->data;
    auto d_w_vel = w_vel->data;

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

#pragma acc data present(d_out[:bsize], d_in[:bsize], d_u_vel[:bsize], d_v_vel[:bsize], d_w_vel[:bsize])
    {
        const size_t Nx = domain->GetNx(out->GetLevel());    //due to unnecessary parameter passing of *this
        const size_t Ny = domain->GetNy(out->GetLevel());

        const real dx = domain->Getdx(out->GetLevel());    //due to unnecessary parameter passing of *this
        const real dy = domain->Getdy(out->GetLevel());
        const real dz = domain->Getdz(out->GetLevel());

        const real rdx = 1. / dx; //due to unnecessary parameter passing of *this
        const real rdy = 1. / dy;
        const real rdz = 1. / dz;

        const real dtx = m_dt * rdx;
        const real dty = m_dt * rdy;
        const real dtz = m_dt * rdz;

#pragma acc parallel loop independent present(d_out[:bsize], d_in[:bsize], d_u_vel[:bsize], d_v_vel[:bsize], d_w_vel[:bsize], d_iList[:bsize_i]) async
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            d_out[i] = d_in[i] - dtx * d_u_vel[i] * (d_in[i] - d_in[i - 1])\
                               - dty * d_v_vel[i] * (d_in[i] - d_in[i - Nx])\
                               - dtz * d_w_vel[i] * (d_in[i] - d_in[i - Nx * Ny]);
        }

        boundary->applyBoundary(d_out, type, sync);

        if (sync) {
#pragma acc wait
        }

    }//end data region
}
