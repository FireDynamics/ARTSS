/// \file       SLAdvect.cpp
/// \brief      Solves advection equation via unconditionally stable Semi-Lagrangian approach
/// \date       Aug 23, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

//==================================== Semi Lagrangian Advection ======================================
// ***************************************************************************************
/// \brief  solves advection \f$ \partial_t \phi_1 = - (u \cdot \nabla) \phi_0 \f$ via unconditionally
///         stable semi-Lagrangian approach (backtrace and linear interpolation)
// ***************************************************************************************

#include <cmath>

#ifdef _OPENACC
#include <accelmath.h>
#endif

#include <algorithm>
#include "SLAdvect.h"
#include "../utility/Parameters.h"
#include "../boundary/BoundaryController.h"
#include "../Domain.h"

// ==================================== Constructor ====================================
// ***************************************************************************************
SLAdvect::SLAdvect() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    auto params = Parameters::getInstance();
    m_dt = params->get_real("physical_parameters/dt");
}

// ***************************************************************************************
/// \brief  solves advection \f$ \partial_t \phi_1 = - (u \cdot \nabla) \phi_0 \f$ via unconditionally
///     stable semi-Lagrangian approach (backtrace and linear interpolation)
/// \param  out   output pointer
/// \param  in    input pointer
/// \param  u_vel x -velocity
/// \param  v_vel y -velocity
/// \param  w_vel z -velocity
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void SLAdvect::advect(Field *out, Field *in, const Field *u_vel, const Field *v_vel, const Field *w_vel, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    size_t bsize = domain->get_size(out->get_level());
    FieldType type = out->get_type();

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
        const size_t Nx = domain->get_Nx(out->get_level());    // due to unnecessary parameter passing of *this
        const size_t Ny = domain->get_Ny(out->get_level());

        const real dx = domain->get_dx(out->get_level());    // due to unnecessary parameter passing of *this
        const real dy = domain->get_dy(out->get_level());
        const real dz = domain->get_dz(out->get_level());

        const real dt = m_dt;

        const real dtx = dt / dx;
        const real dty = dt / dy;
        const real dtz = dt / dz;

        const real epsilon_x = 1e-6; //dx / dt;
        const real epsilon_y = 1e-6; //dy / dt;
        const real epsilon_z = 1e-6; //dz / dt;

        // start indices for computational domain of inner cells
        auto i_start = static_cast<long int> (domain->get_index_x1());
        auto j_start = static_cast<long int> (domain->get_index_y1());
        auto k_start = static_cast<long int> (domain->get_index_z1());
        auto i_end = static_cast<long int> (domain->get_index_x2());
        auto j_end = static_cast<long int> (domain->get_index_y2());
        auto k_end = static_cast<long int> (domain->get_index_z2());

#pragma acc parallel loop independent present(d_out[:bsize], d_in[:bsize], d_u_vel[:bsize], d_v_vel[:bsize], d_w_vel[:bsize], d_iList[:bsize_i]) async
        for (size_t l = 0; l < bsize_i; ++l) {
            const size_t idx = d_iList[l];
            auto k = static_cast<long int> (getCoordinateK(idx, Nx, Ny));
            auto j = static_cast<long int> (getCoordinateJ(idx, Nx, Ny, k));
            auto i = static_cast<long int> (getCoordinateI(idx, Nx, Ny, j, k));

            // TODO: backtracking may be outside the computational region, it is not a reasonable solution to cut the vector; idea: enlarge ghost cell to CFL*dx
            // Linear Trace Back
            real Ci = dtx * d_u_vel[idx];
            real Cj = dty * d_v_vel[idx];
            real Ck = dtz * d_w_vel[idx];

            // Calculation of horizontal indices and interpolation weights
            long int i0 = i;
            long int i1 = i;
            real r = 1;

            if (Ci > epsilon_x) {
                i0 = std::max(i_start, (i - static_cast<long int>(Ci)));
                i1 = i0 - 1;
                r = fmod(Ci, 1); // closer to i1 if r is closer to 1
            } else if (Ci < -epsilon_x) {
                i0 = std::min(i_end, i - static_cast<long int>(Ci));
                i1 = i0 + 1;
                r = fmod(-Ci, 1); // closer to i1 if r is closer to 1
            }

            // Calculation of vertical indices and interpolation weights
            long int j0 = j;
            long int j1 = j;
            real s = 1;

            if (Cj > epsilon_y) {
                j0 = std::max(j_start, j - static_cast<long int>(Cj));
                j1 = j0 - 1;
                s = fmod(Cj, 1);
            } else if (Cj < -epsilon_y){
                j0 = std::min(j_end, j - static_cast<long int>(Cj));
                j1 = j0 + 1;
                s = fmod(-Cj, 1);
            }

            // Calculation of depth indices and interpolation weights
            long int k0 = k;
            long int k1 = k;
            real t = 1;

            if (Ck > epsilon_z) {
                k0 = std::max(k_start, k - static_cast<long int>(Ck));
                k1 = k0 - 1;
                t = fmod(Ck, 1);
            } else if (Ck < -epsilon_z) {
                k0 = std::min(k_end, k - static_cast<long int>(Ck));
                k1 = k0 + 1;
                t = fmod(-Ck, 1);
            }

            // Trilinear Interpolation
            size_t idx_000 = IX(i0, j0, k0, Nx, Ny);
            auto d_000 = d_in[idx_000];

            size_t idx_100 = IX(i1, j0, k0, Nx, Ny);
            auto d_100 = d_in[idx_100];

            size_t idx_010 = IX(i0, j1, k0, Nx, Ny);
            auto d_010 = d_in[idx_010];

            size_t idx_110 = IX(i1, j1, k0, Nx, Ny);
            auto d_110 = d_in[idx_110];

            size_t idx_001 = IX(i0, j0, k1, Nx, Ny);
            auto d_001 = d_in[idx_001];

            size_t idx_101 = IX(i1, j0, k1, Nx, Ny);
            auto d_101 = d_in[idx_101];

            size_t idx_011 = IX(i0, j1, k1, Nx, Ny);
            auto d_011 = d_in[idx_011];

            size_t idx_111 = IX(i1, j1, k1, Nx, Ny);
            auto d_111 = d_in[idx_111];

            auto r100 = d_000 + r * (d_100 - d_000);
            auto r110 = d_010 + r * (d_110 - d_010);
            auto r101 = d_001 + r * (d_101 - d_001);
            auto r111 = d_011 + r * (d_111 - d_011);

            auto s110 = r100 + s * (r110 - r100);
            auto s111 = r101 + s * (r111 - r101);

            auto tmp = s110 + t * (s111 - s110); // row-major
            std::vector<real> values = {d_000, d_001, d_010, d_011, d_100, d_101, d_110, d_111};
            std::vector<long int> coord_i = {i0, i0, i0, i0, i1, i1, i1, i1};
            std::vector<long int> coord_j = {j0, j0, j1, j1, j0, j0, j1, j1};
            std::vector<long int> coord_k = {k0, k1, k0, k1, k0, k1, k0, k1};
            for (long unsigned int v = 0; v < values.size(); v++){
                if (values[v] == 11111.){
                    m_logger->warn("in sladvect, cell ({}|{}|{}) was used with original cell ({}|{}|{})", coord_i[v], coord_j[v], coord_k[v], i, j, k);
                }
            }
            //if (d_000 == d_001 && d_000 == d_010 && d_000 == d_011 && d_000 == d_100 && d_000 == d_101 && d_000 == d_110 && d_000 == d_111) {
            //    tmp = d_000;
            //}
            d_out[idx] = tmp;
        }
        Utility::log_minimum(out, "in temperature advect before boundary", "advect");
        boundary->applyBoundary(d_out, type, sync);
        Utility::log_minimum(out, "in temperature advect after boundary", "advect");

        if (sync) {
#pragma acc wait
        }

    }// end data region
}
