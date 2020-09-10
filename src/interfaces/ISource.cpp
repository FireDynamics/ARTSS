/// \file       ISource.cpp
/// \brief      Interface for adding sources
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>

#ifdef _OPENACC
#include <accelmath.h>
#endif

#include "ISource.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"

//======================================== Sources ====================================
//======================================== Force ======================================
// ***************************************************************************************
/// \brief  Buoyancy Force in momentum equation
/// \param  out buoyancy force
/// \param  in  temperature
/// \param  in_temperature_ambient ambient temperature
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ISource::buoyancy_force(Field *out, const Field *in, const Field *in_temperature_ambient, bool sync) {

    // local variables and parameters for GPU
    auto d_out = out->data;
    auto d_in = in->data;
    auto d_ina = in_temperature_ambient->data;

    auto bsize = Domain::getInstance()->get_size(out->get_level());

    auto params = Parameters::getInstance();

    real beta = params->get_real("physical_parameters/beta");

    real g = params->get_real("physical_parameters/g");

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t *d_bList = boundary->get_boundaryList_level_joined();

    auto bsize_i = boundary->getSize_innerList();
    auto bsize_b = boundary->getSize_boundaryList();

#pragma acc data present(d_iList[:bsize_i], d_bList[:bsize_b], d_out[:bsize], d_in[:bsize], d_ina[:bsize])
    {
        // inner cells
#pragma acc kernels async
#pragma acc loop independent
        for (size_t i = 0; i < bsize_i; ++i) {
            const size_t idx = d_iList[i];
            d_out[idx] = -beta * (d_in[idx] - d_ina[idx]) * g;
        }

        // boundary cells
#pragma acc kernels async
#pragma acc loop independent
        for (size_t i = 0; i < bsize_b; ++i) {
            const size_t idx = d_bList[i];
            d_out[idx] = -beta * (d_in[idx] - d_ina[idx]) * g;
        }

        if (sync) {
#pragma acc wait
        }

    }// end data region
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
void ISource::dissipate(Field *out, const Field *in_u, const Field *in_v, const Field *in_w, bool sync) {

    // local variables and parameters for GPU
    auto d_out = out->data;
    auto d_inu = in_u->data;
    auto d_inv = in_v->data;
    auto d_inw = in_w->data;

    auto domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(out->get_level());
    size_t Ny = domain->get_Ny(out->get_level());

    real dx = domain->get_dx(out->get_level());
    real dy = domain->get_dy(out->get_level());
    real dz = domain->get_dz(out->get_level());
    auto rdx = 1. / dx;
    auto rdy = 1. / dy;
    auto rdz = 1. / dz;

    auto params = Parameters::getInstance();

    real dt = params->get_real("physical_parameters/dt");
    real nu = params->get_real("physical_parameters/nu");

    auto size = Domain::getInstance()->get_size(out->get_level());
    auto type = out->get_type();

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

#pragma acc data present(d_out[:size], d_inu[:size], d_inv[:size], d_inw[:size], d_iList[:bsize_i])
    {
        // inner
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            real out_h = nu * (2 * (0.5 * rdx * (d_inu[i + 1] - d_inu[i - 1])) * (0.5 * rdx * (d_inu[i + 1] - d_inu[i - 1]))
                               + 2 * (0.5 * rdy * (d_inv[i + Nx] - d_inv[i - Nx])) * (0.5 * rdy * (d_inv[i + Nx] - d_inv[i - Nx]))
                               + 2 * (0.5 * rdz * (d_inw[i + Nx * Ny] - d_inw[i - Nx * Ny])) * (0.5 * rdz * (d_inw[i + Nx * Ny] - d_inw[i - Nx * Ny]))
                               + ((0.5 * rdx * (d_inv[i + 1] - d_inv[i - 1])) + (0.5 * rdy * (d_inu[i + Nx] - d_inu[i - Nx])))
                                 * ((0.5 * rdx * (d_inv[i + 1] - d_inv[i - 1])) + (0.5 * rdy * (d_inu[i + Nx] - d_inu[i - Nx])))
                               + ((0.5 * rdy * (d_inw[i + Nx] - d_inw[i - Nx])) + (0.5 * rdz * (d_inv[i + Nx * Ny] - d_inv[i - Nx * Ny])))
                                 * ((0.5 * rdy * (d_inw[i + Nx] - d_inw[i - Nx])) + (0.5 * rdz * (d_inv[i + Nx * Ny] - d_inv[i - Nx * Ny])))
                               + ((0.5 * rdz * (d_inu[i + Nx * Ny] - d_inu[i - Nx * Ny])) + (0.5 * rdx * (d_inw[i + 1] - d_inw[i - 1])))
                                 * ((0.5 * rdz * (d_inu[i + Nx * Ny] - d_inu[i - Nx * Ny])) + (0.5 * rdx * (d_inw[i + 1] - d_inw[i - 1]))));
            d_out[i] += dt * out_h;
        }

        // boundaries
        boundary->applyBoundary(d_out, type, sync);

        if (sync) {
#pragma acc wait
        }
    }//end data region
}
