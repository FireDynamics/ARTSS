/// \file 		SLAdvect.cpp
/// \brief 		Solves advection equation via unconditionally stable Semi-Langrangian approach
/// \date 		Aug 23, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

//==================================== Semi Lagrangian Advection ======================================
// ***************************************************************************************
/// \brief  solves advection \f$ \partial_t \phi_1 = - (u \cdot \nabla) \phi_0 \f$ via unconditionally
///			stable semi-Langrangian approach (backtrace and linear interpolation)
// ***************************************************************************************

#ifdef _OPENACC
#include <accelmath.h>
#else
#include <cmath>
#endif
#include "SLAdvect.h"
#include "../utility/Parameters.h"
#include "../boundary/BoundaryController.h"
#include "../Domain.h"

// ==================================== Constructor ====================================
// ***************************************************************************************
SLAdvect::SLAdvect() {

    auto params = Parameters::getInstance();
    m_dt = params->getReal("physical_parameters/dt");
}

// ***************************************************************************************
/// \brief  solves advection \f$ \partial_t \phi_1 = - (u \cdot \nabla) \phi_0 \f$ via unconditionally
///			stable semi-Langrangian approach (backtrace and linear interpolation)
/// \param  out		output pointer
/// \param	in		input pointer
/// \param	u_vel	x -velocity
/// \param	v_vel	y -velocity
/// \param	w_vel	z -velocity
/// \param  sync	synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void SLAdvect::advect(Field *out, Field *in, const Field *u_vel, const Field *v_vel, const Field *w_vel, bool sync) {

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

    auto bsize_i = boundary->getSize_innerList();
    size_t *d_iList = boundary->get_innerList_level_joined();

#pragma acc data present(d_out[:bsize], d_in[:bsize], d_u_vel[:bsize], d_v_vel[:bsize], d_w_vel[:bsize])
    {
        const size_t Nx = domain->GetNx(out->GetLevel());    //due to unnecessary parameter passing of *this
        const size_t Ny = domain->GetNy(out->GetLevel());

        const real dx = domain->Getdx(out->GetLevel());    //due to unnecessary parameter passing of *this
        const real dy = domain->Getdy(out->GetLevel());
        const real dz = domain->Getdz(out->GetLevel());

        const real rdx = 1. / dx;        //due to unnecessary parameter passing of *this
        const real rdy = 1. / dy;
        const real rdz = 1. / dz;

        const real dt = m_dt;

        const real dtx = dt * rdx;
        const real dty = dt * rdy;
        const real dtz = dt * rdz;

        //start indices for computational domain of inner cells
        long i_start = static_cast<long> (domain->GetIndexx1());
        long j_start = static_cast<long> (domain->GetIndexy1());
        long k_start = static_cast<long> (domain->GetIndexz1());
        long i_end = static_cast<long> (domain->GetIndexx2());
        long j_end = static_cast<long> (domain->GetIndexy2());
        long k_end = static_cast<long> (domain->GetIndexz2());

#pragma acc parallel loop independent present(d_out[:bsize], d_in[:bsize], d_u_vel[:bsize], d_v_vel[:bsize], d_w_vel[:bsize], d_iList[:bsize_i]) async
        for (size_t l = 0; l < bsize_i; ++l) {
            const size_t idx = d_iList[l];
            const long k = static_cast<long> (getCoordinateK(idx, Nx, Ny));
            const long j = static_cast<long> (getCoordinateJ(idx, Nx, Ny, k));
            const long i = static_cast<long> (getCoordinateI(idx, Nx, Ny, j, k));

            //TODO: backtracking may be outside the computational region, it is not a reasonable solution to cut the vector; idea: enlarge ghost cell to CFL*dx
            // Linear Trace Back
            real Ci = dtx * d_u_vel[idx];
            real Cj = dty * d_v_vel[idx];
            real Ck = dtz * d_w_vel[idx];

            // Calculation of horizontal indices and interpolation weights
            long int i0 = 0;
            long int i1 = 0;
            real r;

            if (Ci > 0) {
                i0 = i - (long int) (Ci);    // to round correctly
                if (i0 < i_start) {
                    i0 = i_start;
                }
                i1 = i0 - 1;
                if (i1 > i_end + 1) {
                    i1 = i_end + 1;
                }
                r = fabs(fmod(Ci, 1));
            } else {
                i1 = i - (long int) (Ci);
                if (i1 < i_start - 1) {
                    i1 = i_start - 1;
                }
                i0 = i1 + 1;
                if (i0 > i_end + 1) {
                    i0 = i_end + 1;
                }
                r = 1 - fabs(fmod(Ci, 1));
            }

            i0 = (size_t) i0;
            i1 = (size_t) i1;

            //Calculation of vertical indices and interpolation weights
            long int j0 = 0;
            long int j1 = 0;
            real s;

            if (Cj > 0) {
                j0 = j - (long int) (Cj);    // to round correctly
                if (j0 < j_start) {
                    j0 = j_start;
                }
                j1 = j0 - 1;
                if (j1 > j_end + 1) {
                    j1 = j_end + 1;
                }
                s = fabs(fmod(Cj, 1));
            } else {
                j1 = j - (long int) (Cj);
                if (j1 < j_start - 1) {
                    j1 = j_start - 1;
                }
                j0 = j1 + 1;
                if (j0 > j_end + 1) {
                    j0 = j_end + 1;
                }
                s = 1 - fabs(fmod(Cj, 1));
            }

            j0 = (size_t) j0;
            j1 = (size_t) j1;

            //Calculation of depth indices and interpolation weights
            long int k0 = 0;
            long int k1 = 0;
            real t;

            if (Ck > 0) {
                k0 = k - (long int) (Ck);    // to round correctly
                if (k0 < k_start) {
                    k0 = k_start;
                }
                k1 = k0 - 1;
                if (k1 > k_end + 1) {
                    k1 = k_end + 1;
                }
                t = fabs(fmod(Ck, 1));
            } else {
                k1 = k - (long int) (Ck);
                if (k1 < k_start - 1) {
                    k1 = k_start - 1;
                }
                k0 = k1 + 1;
                if (k0 > k_end + 1) {
                    k0 = k_end + 1;
                }
                t = 1 - fabs(fmod(Ck, 1));
            }

            k0 = (size_t) k0;
            k1 = (size_t) k1;

            // Trilinear Interpolation
            d_out[idx] = (1. - t) * ((1. - s) * ((1. - r) * d_in[IX(i0, j0, k0, Nx, Ny)] + r * d_in[IX(i1, j0, k0, Nx, Ny)])
                                     + s * ((1. - r) * d_in[IX(i0, j1, k0, Nx, Ny)] + r * d_in[IX(i1, j1, k0, Nx, Ny)]))
                         + t * ((1. - s) * ((1. - r) * d_in[IX(i0, j0, k1, Nx, Ny)] + r * d_in[IX(i1, j0, k1, Nx, Ny)])
                                + s * ((1. - r) * d_in[IX(i0, j1, k1, Nx, Ny)] + r * d_in[IX(i1, j1, k1, Nx, Ny)])); //row-major
        }

        boundary->applyBoundary(d_out, type, sync);

        if (sync) {
#pragma acc wait
        }

    }//end data region
}
