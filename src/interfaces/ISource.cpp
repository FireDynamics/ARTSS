/// \file       ISource.cpp
/// \brief      Interface for adding sources
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "ISource.h"

#include <cmath>
#ifdef _OPENACC
#include <accelmath.h>
#endif

#include "../DomainData.h"
#include "../boundary/BoundaryController.h"
#include "../utility/Parameters.h"

//======================================== Sources ====================================
//======================================== Force ======================================
// ***************************************************************************************
/// \brief  Buoyancy Force in momentum equation
/// \param  out buoyancy force
/// \param  in  temperature
/// \param  in_temperature_ambient ambient temperature
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ISource::buoyancy_force(
        Field &out,
        const Field &in, const Field &in_a,
        bool sync) {
    auto params = Parameters::getInstance();
    real beta = params->get_real("physical_parameters/beta");
    real g = params->get_real("physical_parameters/g");

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_inner_list_level_joined();
    size_t *d_bList = boundary->get_boundary_list_level_joined();

    auto bsize_i = boundary->get_size_inner_list();
    auto bsize_b = boundary->get_size_boundary_list();

#pragma acc data present(d_iList[:bsize_i], d_bList[:bsize_b], out, in, in_a)
    {
        // inner cells
#pragma acc kernels async
#pragma acc loop independent
        for (size_t i = 0; i < bsize_i; ++i) {
            const size_t idx = d_iList[i];
            out[idx] = -beta * (in[idx] - in_a[idx]) * g;
        }

        // boundary cells
#pragma acc kernels async
#pragma acc loop independent
        for (size_t i = 0; i < bsize_b; ++i) {
            const size_t idx = d_bList[i];
            out[idx] = -beta * (in[idx] - in_a[idx]) * g;
        }

        if (sync) {
#pragma acc wait
        }
    }
}

//======================================== Dissipation ====================================
// ***************************************************************************************
/// \brief  Adding dissipation (e.g. via friction) into the equation
/// \param  out   output pointer (\a T)
/// \param  in_u  input pointer (\a x -velocity)
/// \param  in_v  input pointer (\a y -velocity)
/// \param  in_w  input pointer (\a z -velocity)
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ISource::dissipate(
        Field &out,
        const Field &in_u, const Field &in_v, const Field &in_w,
        bool sync) {
    auto domain = DomainData::getInstance();
    size_t Nx = domain->get_Nx();
    size_t Ny = domain->get_Ny();

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();
    real reciprocal_dx = 1. / dx;
    real reciprocal_dy = 1. / dy;
    real reciprocal_dz = 1. / dz;

    auto params = Parameters::getInstance();

    real dt = params->get_real("physical_parameters/dt");
    real nu = params->get_real("physical_parameters/nu");

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_inner_list_level_joined();
    auto bsize_i = boundary->get_size_inner_list();

    size_t neighbour_i = 1;
    size_t neighbour_j = Nx;
    size_t neighbour_k = Nx * Ny;
#pragma acc data present(out, in_u, in_v, in_w, d_iList[:bsize_i])
    {
        // inner
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            real out_h = nu * (2 * (0.5 * reciprocal_dx * (in_u[i + neighbour_i] - in_u[i - neighbour_i]))
                                 * (0.5 * reciprocal_dx * (in_u[i + neighbour_i] - in_u[i - neighbour_i]))
                             + 2 * (0.5 * reciprocal_dy * (in_v[i + neighbour_j] - in_v[i - neighbour_j]))
                                 * (0.5 * reciprocal_dy * (in_v[i + neighbour_j] - in_v[i - neighbour_j]))
                             + 2 * (0.5 * reciprocal_dz * (in_w[i + neighbour_k] - in_w[i - neighbour_k]))
                                 * (0.5 * reciprocal_dz * (in_w[i + neighbour_k] - in_w[i - neighbour_k]))
                                + ((0.5 * reciprocal_dx * (in_v[i + neighbour_i] - in_v[i - neighbour_i]))
                                +  (0.5 * reciprocal_dy * (in_u[i + neighbour_j] - in_u[i - neighbour_j])))
                                * ((0.5 * reciprocal_dx * (in_v[i + neighbour_i] - in_v[i - neighbour_i]))
                                +  (0.5 * reciprocal_dy * (in_u[i + neighbour_j] - in_u[i - neighbour_j])))
                                + ((0.5 * reciprocal_dy * (in_w[i + neighbour_j] - in_w[i - neighbour_j]))
                                +  (0.5 * reciprocal_dz * (in_v[i + neighbour_k] - in_v[i - neighbour_k])))
                                * ((0.5 * reciprocal_dy * (in_w[i + neighbour_j] - in_w[i - neighbour_j]))
                                +  (0.5 * reciprocal_dz * (in_v[i + neighbour_k] - in_v[i - neighbour_k])))
                                + ((0.5 * reciprocal_dz * (in_u[i + neighbour_k] - in_u[i - neighbour_k]))
                                +  (0.5 * reciprocal_dx * (in_w[i + neighbour_i] - in_w[i - neighbour_i])))
                                * ((0.5 * reciprocal_dz * (in_u[i + neighbour_k] - in_u[i - neighbour_k]))
                                +  (0.5 * reciprocal_dx * (in_w[i + neighbour_i] - in_w[i - neighbour_i]))));
            out[i] += dt * out_h;
        }

        // boundaries
        boundary->apply_boundary(out, sync);

        if (sync) {
#pragma acc wait
        }
    }
}

