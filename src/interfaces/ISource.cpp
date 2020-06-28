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

    auto bsize = Domain::getInstance()->get_size(out->GetLevel());

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

//===================================== Energy Source ====================================
// ***************************************************************************************
/// \brief  Manufactured (MMS) energy source in energy equation
/// \param  out   energy source
/// \param  t   time
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ISource::buoyancy_ST_MMS(Field *out, real t, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    auto d_out = out->data;
    auto bsize = domain->get_size(out->GetLevel());

    size_t Nx = domain->get_Nx(out->GetLevel());
    size_t Ny = domain->get_Ny(out->GetLevel());

    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();

    real dx = domain->get_dx(out->GetLevel());
    real dy = domain->get_dy(out->GetLevel());

    auto params = Parameters::getInstance();

    real nu = params->get_real("physical_parameters/nu");
    real beta = params->get_real("physical_parameters/beta");
    real kappa = params->get_real("physical_parameters/kappa");
    real g = params->get_real("physical_parameters/g");
    real rhoa = params->get_real("initial_conditions/rhoa");
    real rbeta = 1. / beta;
    real rg = 1. / g;
    real c_nu = 2 * nu * M_PI * M_PI - 1;
    real c_kappa = 2 * kappa * M_PI * M_PI - 1;

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t *d_bList = boundary->get_boundaryList_level_joined();

    auto bsize_i = boundary->getSize_innerList();
    auto bsize_b = boundary->getSize_boundaryList();

    size_t i, j, k;

    // inner cells
#pragma acc data present(d_iList[:bsize_i], d_bList[:bsize_b], d_out[:bsize])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = 0; l < bsize_i; ++l) {
            const size_t idx = d_iList[l];
            k = idx / (Nx * Ny);
            j = (idx - k * Nx * Ny) / Nx;
            i = idx - k * Nx * Ny - j * Nx;
            d_out[idx] = rhoa * rbeta * rg * 2 * c_nu * c_kappa * exp(-t) * std::sin(M_PI * (xi(i, X1, dx) + yj(j, Y1, dy)));
        }

        // boundary cells
#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = 0; l < bsize_b; ++l) {
            const size_t idx = d_bList[l];
            k = idx / (Nx * Ny);
            j = (idx - k * Nx * Ny) / Nx;
            i = idx - k * Nx * Ny - j * Nx;
            d_out[idx] = rhoa * rbeta * rg * 2 * c_nu * c_kappa * exp(-t) * std::sin(M_PI * (xi(i, X1, dx) + yj(j, Y1, dy)));
        }
        if (sync) {
#pragma acc wait
        }
    }// end data region
}

// ***************************************************************************************
/// \brief  Volumetric Gaussian temperature source in energy equation
/// \param  out   energy source
/// \param  HRR   total heat release rate
/// \param  cp    heat capacity
/// \param  x0    center of Gaussian (x-direction)
/// \param  y0    center of Gaussian (y-direction)
/// \param  z0    center of Gaussian (z-direction)
/// \param  sigma Radius of Gaussian
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ISource::gauss(Field *out, real HRR, real cp, real x0, real y0, real z0, real sigma_x, real sigma_y, real sigma_z, bool sync) {

    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    auto d_out = out->data;
    auto bsize = domain->get_size(out->GetLevel());

    size_t Nx = domain->get_Nx(out->GetLevel());
    size_t Ny = domain->get_Ny(out->GetLevel());

    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    real dx = domain->get_dx(out->GetLevel());
    real dy = domain->get_dy(out->GetLevel());
    real dz = domain->get_dz(out->GetLevel());

    //get parameters for Gaussian
    real sigma_x_2 = 2 * sigma_x * sigma_x;
    real r_sigma_x_2 = 1. / sigma_x_2;
    real sigma_y_2 = 2 * sigma_y * sigma_y;
    real r_sigma_y_2 = 1. / sigma_y_2;
    real sigma_z_2 = 2 * sigma_z * sigma_z;
    real r_sigma_z_2 = 1. / sigma_z_2;

    //set Gaussian to cells
    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();

    auto bsize_i = boundary->getSize_innerList();

    size_t i, j, k;
    real V = 0.;
    real HRRrV;

    // inner cells
#pragma acc data present(d_iList[:bsize_i], d_out[:bsize])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = 0; l < bsize_i; ++l) {
            const size_t idx = d_iList[l];
            k = idx / (Nx * Ny);
            j = (idx - k * Nx * Ny) / Nx;
            i = idx - k * Nx * Ny - j * Nx;

            real expr = std::exp(-(r_sigma_x_2 * ((xi(i, X1, dx) - x0) * (xi(i, X1, dx) - x0))
                                    + r_sigma_y_2 * ((yj(j, Y1, dy) - y0) * (yj(j, Y1, dy) - y0))
                                    + r_sigma_z_2 * ((zk(k, Z1, dz) - z0) * (zk(k, Z1, dz) - z0))));
            V += expr * dx * dy * dz;
        }

#pragma acc wait

        HRRrV = HRR / V;        //in case of concentration Ys*HRR
        real rcp = 1. / cp;    // to get [K/s] for energy equation (d_t T), rho:=1, otherwise *1/rho; in case of concentration 1/Hc to get kg/m^3s

#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = 0; l < bsize_i; ++l) {
            const size_t idx = d_iList[l];
            k = idx / (Nx * Ny);
            j = (idx - k * Nx * Ny) / Nx;
            i = idx - k * Nx * Ny - j * Nx;

            real expr = std::exp(-(r_sigma_x_2 * ((xi(i, X1, dx) - x0) * (xi(i, X1, dx) - x0))
                                    + r_sigma_y_2 * ((yj(j, Y1, dy) - y0) * (yj(j, Y1, dy) - y0))
                                    + r_sigma_z_2 * ((zk(k, Z1, dz) - z0) * (zk(k, Z1, dz) - z0))));
            d_out[idx] = HRRrV * rcp * expr;
        }

        if (sync) {
#pragma acc wait
        }
    } //end data region
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
    size_t Nx = domain->get_Nx(out->GetLevel());
    size_t Ny = domain->get_Ny(out->GetLevel());

    real dx = domain->get_dx(out->GetLevel());
    real dy = domain->get_dy(out->GetLevel());
    real dz = domain->get_dz(out->GetLevel());
    auto rdx = 1. / dx;
    auto rdy = 1. / dy;
    auto rdz = 1. / dz;

    auto params = Parameters::getInstance();

    real dt = params->get_real("physical_parameters/dt");
    real nu = params->get_real("physical_parameters/nu");

    auto size = Domain::getInstance()->get_size(out->GetLevel());
    auto type = out->GetType();

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
            real out_h = nu * (2 * (0.5 * rdx * (d_inu[i + 1      ] - d_inu[i - 1      ])) * (0.5 * rdx * (d_inu[i + 1      ] - d_inu[i - 1 ]))
                             + 2 * (0.5 * rdy * (d_inv[i + Nx     ] - d_inv[i - Nx     ])) * (0.5 * rdy * (d_inv[i + Nx     ] - d_inv[i - Nx]))
                             + 2 * (0.5 * rdz * (d_inw[i + Nx * Ny] - d_inw[i - Nx * Ny])) * (0.5 * rdz * (d_inw[i + Nx * Ny] - d_inw[i - Nx * Ny]))
                                + ((0.5 * rdx * (d_inv[i + 1      ] - d_inv[i - 1      ])) + (0.5 * rdy * (d_inu[i + Nx     ] - d_inu[i - Nx])))
                                * ((0.5 * rdx * (d_inv[i + 1      ] - d_inv[i - 1      ])) + (0.5 * rdy * (d_inu[i + Nx     ] - d_inu[i - Nx])))
                                + ((0.5 * rdy * (d_inw[i + Nx     ] - d_inw[i - Nx     ])) + (0.5 * rdz * (d_inv[i + Nx * Ny] - d_inv[i - Nx * Ny])))
                                * ((0.5 * rdy * (d_inw[i + Nx     ] - d_inw[i - Nx     ])) + (0.5 * rdz * (d_inv[i + Nx * Ny] - d_inv[i - Nx * Ny])))
                                + ((0.5 * rdz * (d_inu[i + Nx * Ny] - d_inu[i - Nx * Ny])) + (0.5 * rdx * (d_inw[i + 1      ] - d_inw[i - 1])))
                                * ((0.5 * rdz * (d_inu[i + Nx * Ny] - d_inu[i - Nx * Ny])) + (0.5 * rdx * (d_inw[i + 1      ] - d_inw[i - 1]))));
            d_out[i] += dt * out_h;
        }

        // boundaries
        boundary->applyBoundary(d_out, type, sync);

        if (sync) {
#pragma acc wait
        }
    }//end data region
}
