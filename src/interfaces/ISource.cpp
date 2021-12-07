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
#include "../boundary/DomainController.h"

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
        Settings::Settings const &settings,
        Field &out,
        const Field &in, const Field &in_a,
        bool sync) {
    real beta = settings.get_real("physical_parameters/beta");
    real g = settings.get_real("physical_parameters/g");

    auto boundary = DomainController::getInstance();
    size_t size_domain_list = boundary->get_slice_size_domain_list_level_joined(0);
    size_t *domain_list = boundary->get_domain_list_level_joined();

#pragma acc data present(domain_list[:size_domain_list], out, in, in_a)
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t i = 0; i < size_domain_list; ++i) {
            const size_t idx = domain_list[i];
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
        Settings::Settings const &settings,
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

    real dt = settings.get_real("physical_parameters/dt");
    real nu = settings.get_real("physical_parameters/nu");

    auto boundary = DomainController::getInstance();
    size_t *d_iList = boundary->get_domain_inner_list_level_joined();
    auto bsize_i = boundary->get_size_domain_inner_list_level_joined(0);

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

