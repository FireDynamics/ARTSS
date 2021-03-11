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

#include "../Domain.h"
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
void ISource::buoyancy_force(Field &out, Field const &in,
        Field const &ina, bool sync) {
    auto params = Parameters::getInstance();
    real beta = params->get_real("physical_parameters/beta");
    real g = params->get_real("physical_parameters/g");

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t *d_bList = boundary->get_boundaryList_level_joined();

    auto bsize_i = boundary->getSize_innerList();
    auto bsize_b = boundary->getSize_boundaryList();

#pragma acc data present(d_iList[:bsize_i], d_bList[:bsize_b], out, in, ina)
    {
        // inner cells
#pragma acc kernels async
#pragma acc loop independent
        for (size_t i = 0; i < bsize_i; ++i) {
            const size_t idx = d_iList[i];
            out[idx] = -beta * (in[idx] - ina[idx]) * g;
        }

        // boundary cells
#pragma acc kernels async
#pragma acc loop independent
        for (size_t i = 0; i < bsize_b; ++i) {
            const size_t idx = d_bList[i];
            out[idx] = -beta * (in[idx] - ina[idx]) * g;
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
void ISource::dissipate(Field &out,
        Field const &in_u, Field const &in_v, Field const &in_w, bool sync) {
    auto domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(out.get_level());
    size_t Ny = domain->get_Ny(out.get_level());

    real dx = domain->get_dx(out.get_level());
    real dy = domain->get_dy(out.get_level());
    real dz = domain->get_dz(out.get_level());
    auto rdx = 1. / dx;
    auto rdy = 1. / dy;
    auto rdz = 1. / dz;

    auto params = Parameters::getInstance();

    real dt = params->get_real("physical_parameters/dt");
    real nu = params->get_real("physical_parameters/nu");

    auto type = out.get_type();

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

#pragma acc data present(out, in_u, in_v, in_w, d_iList[:bsize_i])
    {
        // inner
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            real out_h = nu * (2 * (0.5 * rdx * (in_u[i + 1      ] - in_u[i - 1      ])) * (0.5 * rdx * (in_u[i + 1      ] - in_u[i - 1 ]))
                             + 2 * (0.5 * rdy * (in_v[i + Nx     ] - in_v[i - Nx     ])) * (0.5 * rdy * (in_v[i + Nx     ] - in_v[i - Nx]))
                             + 2 * (0.5 * rdz * (in_w[i + Nx * Ny] - in_w[i - Nx * Ny])) * (0.5 * rdz * (in_w[i + Nx * Ny] - in_w[i - Nx * Ny]))
                                + ((0.5 * rdx * (in_v[i + 1      ] - in_v[i - 1      ])) + (0.5 * rdy * (in_u[i + Nx     ] - in_u[i - Nx])))
                                * ((0.5 * rdx * (in_v[i + 1      ] - in_v[i - 1      ])) + (0.5 * rdy * (in_u[i + Nx     ] - in_u[i - Nx])))
                                + ((0.5 * rdy * (in_w[i + Nx     ] - in_w[i - Nx     ])) + (0.5 * rdz * (in_v[i + Nx * Ny] - in_v[i - Nx * Ny])))
                                * ((0.5 * rdy * (in_w[i + Nx     ] - in_w[i - Nx     ])) + (0.5 * rdz * (in_v[i + Nx * Ny] - in_v[i - Nx * Ny])))
                                + ((0.5 * rdz * (in_u[i + Nx * Ny] - in_u[i - Nx * Ny])) + (0.5 * rdx * (in_w[i + 1      ] - in_w[i - 1])))
                                * ((0.5 * rdz * (in_u[i + Nx * Ny] - in_u[i - Nx * Ny])) + (0.5 * rdx * (in_w[i + 1      ] - in_w[i - 1]))));
            d_out[i] += dt * out_h;
        }

        // boundaries
        boundary->applyBoundary(out.data, type, sync);

        if (sync) {
#pragma acc wait
        }
    }
}

