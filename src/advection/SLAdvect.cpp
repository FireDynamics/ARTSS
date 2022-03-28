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

#include "SLAdvect.h"

#include <cmath>

#ifdef _OPENACC
#include <accelmath.h>
#endif

#include "../domain/DomainController.h"
#include "../domain/DomainData.h"

// ***************************************************************************************
/// \brief  solves advection \f$ \partial_t \phi_1 = - (u \cdot \nabla) \phi_0 \f$ via
///         unconditionally stable semi-Lagrangian approach (backtrace and linear interpolation)
/// \param  out   output pointer
/// \param  in    input pointer
/// \param  u_vel x -velocity
/// \param  v_vel y -velocity
/// \param  w_vel z -velocity
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void SLAdvect::advect(Field &out, const Field &in,
                      const Field &u_vel, const Field &v_vel, const Field &w_vel,
                      bool sync) {
    auto domain_data = DomainData::getInstance();
    auto domain_controller = DomainController::getInstance();

    auto bsize_i = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t *d_inner_list = domain_controller->get_domain_inner_list_level_joined();

#pragma acc data present(out, in, u_vel, v_vel, w_vel, d_inner_list[:bsize_i])
    {
        const size_t Nx = domain_data->get_Nx();
        const size_t Ny = domain_data->get_Ny();

        const real dx = domain_data->get_dx();
        const real dy = domain_data->get_dy();
        const real dz = domain_data->get_dz();

        const real dt = domain_data->get_physical_parameters().dt;

        const real dtx = dt / dx;
        const real dty = dt / dy;
        const real dtz = dt / dz;

        const real epsilon_x = 1e-6; //dx / dt;
        const real epsilon_y = 1e-6; //dy / dt;
        const real epsilon_z = 1e-6; //dz / dt;

#pragma acc loop independent
        for (size_t l = 0; l < bsize_i; ++l) {
            const size_t idx = d_inner_list[l];
            auto k = static_cast<long int> (getCoordinateK(idx, Nx, Ny));
            auto j = static_cast<long int> (getCoordinateJ(idx, Nx, Ny, k));
            auto i = static_cast<long int> (getCoordinateI(idx, Nx, Ny, j, k));

            // TODO: backtracking may be outside the computational region, it is not a reasonable solution to cut the vector; idea: enlarge ghost cell to CFL*dx
            // Linear Trace Back
            real Ci = dtx * u_vel[idx];
            real Cj = dty * v_vel[idx];
            real Ck = dtz * w_vel[idx];

            // Calculation of horizontal indices and interpolation weights
            auto [i0, i1, r] = calculate_backward_index(CoordinateAxis::X, i, epsilon_x, Ci);
            auto [j0, j1, s] = calculate_backward_index(CoordinateAxis::Y, j, epsilon_y, Cj);
            auto [k0, k1, t] = calculate_backward_index(CoordinateAxis::Z, k, epsilon_z, Ck);

            // Trilinear Interpolation
            size_t idx_000 = IX(i0, j0, k0, Nx, Ny);
            auto d_000 = in[idx_000];

            size_t idx_100 = IX(i1, j0, k0, Nx, Ny);
            auto d_100 = in[idx_100];

            size_t idx_010 = IX(i0, j1, k0, Nx, Ny);
            auto d_010 = in[idx_010];

            size_t idx_110 = IX(i1, j1, k0, Nx, Ny);
            auto d_110 = in[idx_110];

            size_t idx_001 = IX(i0, j0, k1, Nx, Ny);
            auto d_001 = in[idx_001];

            size_t idx_101 = IX(i1, j0, k1, Nx, Ny);
            auto d_101 = in[idx_101];

            size_t idx_011 = IX(i0, j1, k1, Nx, Ny);
            auto d_011 = in[idx_011];

            size_t idx_111 = IX(i1, j1, k1, Nx, Ny);
            auto d_111 = in[idx_111];

            auto r100 = d_000 + r * (d_100 - d_000);
            auto r110 = d_010 + r * (d_110 - d_010);
            auto r101 = d_001 + r * (d_101 - d_001);
            auto r111 = d_011 + r * (d_111 - d_011);

            auto s110 = r100 + s * (r110 - r100);
            auto s111 = r101 + s * (r111 - r101);

            auto tmp = s110 + t * (s111 - s110);  // row-major
            out[idx] = tmp;
        }
        domain_controller->apply_boundary(out, sync);

        if (sync) {
#pragma acc wait
        }
    }
}

return_backtracking_parameters SLAdvect::calculate_backward_index(CoordinateAxis axis, size_t coordinate, real epsilon, real trace_back) {
    size_t i0 = coordinate;
    size_t i1 = coordinate;
    real r = 1;

    auto domain_data = DomainData::getInstance();
    auto start = static_cast<long int>(domain_data->get_start_index_CD(axis));
    auto end = static_cast<long int>(domain_data->get_end_index_CD(axis));

    auto coord = static_cast<long int>(coordinate);

    if (trace_back > epsilon) {
        i0 = std::max(start, (coord - static_cast<long int>(trace_back)));
        i1 = i0 - 1;
        r = fmod(trace_back, 1);  // closer to i1 if r is closer to 1
    } else if (trace_back < -epsilon) {
        i0 = std::min(end, coord - static_cast<long int>(trace_back));
        i1 = i0 + 1;
        r = fmod(-trace_back, 1);  // closer to i1 if r is closer to 1
    }
    return {i0, i1, r};
}

