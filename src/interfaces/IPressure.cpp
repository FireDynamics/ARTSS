/// \file       IPressure.cpp
/// \brief      Interface for pressure method
/// \date       Sep 14, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "IPressure.h"
#include "../utility/Parameters.h"
#include "../DomainData.h"
#include "../boundary/BoundaryController.h"

//======================================== Divergence ====================================
// ***************************************************************************************
/// \brief  calculates divergence \f$ \nabla \cdot u\f$  via central finite differences
/// \param  out   output pointer (\f$ \nabla \cdot u\f$)
/// \param  in_x  input pointer (x -velocity)
/// \param  in_y  input pointer (y -velocity)
/// \param  in_z  input pointer (z -velocity)
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void IPressure::divergence(
        Field &out,
        const Field &in_x, const Field &in_y, const Field &in_z, bool sync) {
    auto domain = DomainData::getInstance();

    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();
    auto dx = domain->get_dx();
    auto dy = domain->get_dy();
    auto dz = domain->get_dz();
    auto reciprocal_dx = 1. / dx;
    auto reciprocal_dy = 1. / dy;
    auto reciprocal_dz = 1. / dz;

    auto boundary = BoundaryController::getInstance();

    size_t *d_inner_list = boundary->get_inner_list_level_joined();
    size_t *d_boundary_list = boundary->get_boundary_list_level_joined();

    auto bsize_i = boundary->get_size_inner_list();
    auto bsize_b = boundary->get_size_boundary_list();

    size_t neighbour_i = 1;
    size_t neighbour_j = Nx;
    size_t neighbour_k = Nx * Ny;
#pragma acc data present(d_inner_list[:bsize_i], d_boundary_list[:bsize_b])
#pragma acc data present(out, in_x, in_y, in_z)
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_inner_list[j];
            real value_in_x = 0.5 * reciprocal_dx * (in_x[i + neighbour_i] - in_x[i - neighbour_i]) ;
            real value_in_y = 0.5 * reciprocal_dy * (in_y[i + neighbour_j] - in_y[i - neighbour_j]) ;
            real value_in_z = 0.5 * reciprocal_dz * (in_z[i + neighbour_k] - in_z[i - neighbour_k]) ;
            out[i] = value_in_x + value_in_y + value_in_z;
        }

// boundaries
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t i = d_boundary_list[j];
            out[i] = 0.;
        }

        if (sync) {
#pragma acc wait
        }
    }
}

//======================================== Projection ====================================
// ***************************************************************************************
/// \brief  subtracts pressure gradient (\f$ \phi_4 = \phi_3 - \nabla p\f$) via finite differences
///     to make \a \f$ \phi_4 \f$ divergence-free
/// \param  out_u output pointer (x -velocity)
/// \param  out_v output pointer (y -velocity)
/// \param  out_w output pointer (z -velocity)
/// \param  in_u  input pointer (x -velocity)
/// \param  in_v  input pointer (y -velocity)
/// \param  in_w  input pointer (z -velocity)
/// \param  in_p  input pointer (pressure)
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void IPressure::projection(
        Field &out_u, Field &out_v, Field &out_w,
        const Field &in_u, const Field &in_v, const Field &in_w,
        const Field &in_p, bool sync) {
    auto domain = DomainData::getInstance();
    // local variables and parameters for GPU
    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();

    auto dx = domain->get_dx();
    auto dy = domain->get_dy();
    auto dz = domain->get_dz();

    auto reciprocal_dx = 1. / dx;
    auto reciprocal_dy = 1. / dy;
    auto reciprocal_dz = 1. / dz;

    auto boundary = BoundaryController::getInstance();

    size_t *d_inner_list = boundary->get_inner_list_level_joined();

    auto bsize_i = boundary->get_size_inner_list();

    size_t neighbour_i = 1;
    size_t neighbour_j = Nx;
    size_t neighbour_k = Nx * Ny;
#pragma acc data present(d_inner_list[:bsize_i])
#pragma acc data present(out_u, out_v, out_w, in_u, in_v, in_w, in_p)
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_inner_list[j];
            out_u[i] = in_u[i] - 0.5 * reciprocal_dx * (in_p[i + neighbour_i] - in_p[i - neighbour_i]);
            out_v[i] = in_v[i] - 0.5 * reciprocal_dy * (in_p[i + neighbour_j] - in_p[i - neighbour_j]);
            out_w[i] = in_w[i] - 0.5 * reciprocal_dz * (in_p[i + neighbour_k] - in_p[i - neighbour_k]);
        }

        // boundaries
        boundary->apply_boundary(out_u, false);
        boundary->apply_boundary(out_v, false);
        boundary->apply_boundary(out_w, false);

        if (sync) {
#pragma acc wait
        }
    }
}
