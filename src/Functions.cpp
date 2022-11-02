/// \file       Functions.cpp
/// \brief      Functions for Initialization
/// \date       Jun 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Functions.h"

#include <cmath>

#ifdef _OPENACC
#include <accelmath.h>
#endif
#include <ctime>
#include <random>

#include "domain/DomainData.h"
#include "domain/DomainController.h"
#include "utility/Utility.h"
#include "interfaces/IRandomField.h"
#include "randomField/UniformRandom.h"

const std::string FunctionNames::beltrami = "Beltrami";
const std::string FunctionNames::buoyancy_mms = "BuoyancyMMS";
const std::string FunctionNames::drift = "Drift";
const std::string FunctionNames::exp_sinus_prod = "ExpSinusProd";
const std::string FunctionNames::exp_sinus_sum = "ExpSinusSum";
const std::string FunctionNames::gauss_bubble = "GaussBubble";
const std::string FunctionNames::hat = "Hat";
const std::string FunctionNames::jet = "Jet";
const std::string FunctionNames::layers = "LayersT";
const std::string FunctionNames::mcdermott = "McDermott";
const std::string FunctionNames::sin_sin_sin = "SinSinSin";
const std::string FunctionNames::uniform = "Uniform";
const std::string FunctionNames::vortex = "Vortex";
const std::string FunctionNames::vortex_y = "VortexY";
const std::string FunctionNames::zero = "Zero";

const std::string class_name = "Functions";
namespace Functions {

// ================================ NS Test - Beltrami IC =================================
// ***************************************************************************************
/// \brief  Initial set up for NS Test - Beltrami
/// \param  out_x  x-velocity
/// \param  out_y  y-velocity
/// \param  out_z  z-velocity
/// \param  out_p  pressure
/// \param  t time
// ***************************************************************************************
void beltrami(Field &out_x, Field &out_y, Field &out_z, Field &out_p, real t, const Settings::initial_conditions::beltrami &beltrami) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;

    const real a = beltrami.a;
    const real d = beltrami.d;
    const real nu = domain_data->get_physical_parameters().nu.value();
    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x, out_y, out_z, out_p) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
        out_x[idx] = -a * (exp(a * xi(coords_i, X1, dx)) * sin(a * yj(coords_j, Y1, dy) + dz) +
                            exp(a * zk(coords_k, Z1, dz)) * cos(a * xi(coords_i, X1, dx) + dy)) * exp(-nu * d * d * t);
        out_y[idx] = -a * (exp(a * yj(coords_j, Y1, dy)) * sin(a * zk(coords_k, Z1, dz) + dx) +
                            exp(a * xi(coords_i, X1, dx)) * cos(a * yj(coords_j, Y1, dy) + dz)) * exp(-nu * d * d * t);
        out_z[idx] = -a * (exp(a * zk(coords_k, Z1, dz)) * sin(a * xi(coords_i, X1, dx) + dy) +
                            exp(a * yj(coords_j, Y1, dy)) * cos(a * zk(coords_k, Z1, dz) + dx)) * exp(-nu * d * d * t);
        out_p[idx] =
                - 0.5 * a * a * (exp(2 * a * xi(coords_i, X1, dx)) + exp(2 * a * yj(coords_j, Y1, dy)) + exp(2 * a * zk(coords_k, Z1, dz)) \
                + 2 * sin(a * xi(coords_i, X1, dx) + dy) * cos(a * zk(coords_k, Z1, dz) + dx) * exp(a * (yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz))) \
                + 2 * sin(a * yj(coords_j, Y1, dy) + dz) * cos(a * xi(coords_i, X1, dx) + dy) * exp(a * (zk(coords_k, Z1, dz) + xi(coords_i, X1, dx))) \
                + 2 * sin(a * zk(coords_k, Z1, dz) + dx) * cos(a * yj(coords_j, Y1, dy) + dz) * exp(a * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy))));
    }
}

// ================================ NS Test - Beltrami IC for p ==========================
// ***************************************************************************************
/// \brief  Initial pressure set up for NS Test - Beltrami
/// \param  out_x  pressure
// ***************************************************************************************
void beltrami_bc_p(Field &out_x, real a) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
        out_x[idx] =
                -0.5 * a * a * (exp(2 * a * xi(coords_i, X1, dx)) + exp(2 * a * yj(coords_j, Y1, dy)) + exp(2 * a * zk(coords_k, Z1, dz)) \
                + 2 * sin(a * xi(coords_i, X1, dx) + dy) * cos(a * zk(coords_k, Z1, dz) + dx) * exp(a * (yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz))) \
                + 2 * sin(a * yj(coords_j, Y1, dy) + dz) * cos(a * xi(coords_i, X1, dx) + dy) * exp(a * (zk(coords_k, Z1, dz) + xi(coords_i, X1, dx))) \
                + 2 * sin(a * zk(coords_k, Z1, dz) + dx) * cos(a * yj(coords_j, Y1, dy) + dz) * exp(a * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy))));
    }
}

// ================================ NS Test - Beltrami IC for u ==========================
// ***************************************************************************************
/// \brief  Initial x-velocity set up for NS Test - Beltrami
/// \param  out_x  x-velocity
/// \param  t time
// ***************************************************************************************
void beltrami_bc_u(Field &out_x, real t, real a, real d, real nu) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
        out_x[idx] = -a * (exp(a * xi(coords_i, X1, dx)) * sin(a * yj(coords_j, Y1, dy) + dz) +
                            exp(a * zk(coords_k, Z1, dz)) * cos(a * xi(coords_i, X1, dx) + dy)) * exp(-nu * d * d * t);
    }
}

// ================================ NS Test - Beltrami IC for v ==========================
// ***************************************************************************************
/// \brief  Initial y-velocity set up for NS Test - Beltrami
/// \param  out_y  y-velocity
/// \param  t time
// ***************************************************************************************
void beltrami_bc_v(Field &out_x, real t, real a, real d, real nu) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
        out_x[idx] = -a * (exp(a * yj(coords_j, Y1, dy)) * sin(a * zk(coords_k, Z1, dz) + dx) +
                            exp(a * xi(coords_i, X1, dx)) * cos(a * yj(coords_j, Y1, dy) + dz)) * exp(-nu * d * d * t);
    }
}

// ================================ NS Test - Beltrami IC for w ==========================
// ***************************************************************************************
/// \brief  Initial z-velocity set up for NS Test - Beltrami
/// \param  out_z  z-velocity
/// \param  t time
// ***************************************************************************************
void beltrami_bc_w(Field &out_x, real t, real a, real d, real nu) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
        out_x[idx] = -a * (exp(a * zk(coords_k, Z1, dz)) * sin(a * xi(coords_i, X1, dx) + dy) +
                            exp(a * yj(coords_j, Y1, dy)) * cos(a * zk(coords_k, Z1, dz) + dx)) * exp(-nu * d * d * t);
    }
}

// ================== NSTemp Test - MMS IC for u,v,w,p,T with buoyancy ===================
// ***************************************************************************************
/// \brief  Initial set up for NSTemp Test - MMS with buoyant force
/// \param  out_x  x-velocity
/// \param  out_y  y-velocity
/// \param  out_z  z-velocity
/// \param  out_p  pressure
/// \param  out_T  temperature
/// \param  t   time
// ***************************************************************************************
void buoyancy_mms(Field &out_x, Field &out_y, Field &out_z, Field &out_p, Field &out_T, real t) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();

    const real nu = domain_data->get_physical_parameters().nu.value();
    const real beta = domain_data->get_physical_parameters().beta;
    const real g = domain_data->get_physical_parameters().g;
    const real rhoa = domain_data->get_physical_parameters().rhoa.value();
    real rbeta = 1. / beta;
    real rg = 1. / g;
    real c = 2 * nu * M_PI * M_PI - 1;
    real rpi = 1. / M_PI;

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_k, coords_i, coords_j;

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x, out_y, out_z, out_p, out_T) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
        out_x[idx] = exp(-t) * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
        out_y[idx] = -exp(-t) * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
        out_z[idx] = 0.;
        out_p[idx] = rhoa * rpi * c * exp(-t) * cos(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
        out_T[idx] = rhoa * rbeta * rg * 2 * c * exp(-t) * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
    }
}

// ===================================== NS Test - IC for u,v,w,p ========================
// ***************************************************************************************
/// \brief  Initial set up for NS Test - Flow around cube or Channel flow with Drift
/// \param  out_x  x-velocity
/// \param  out_y  y-velocity
/// \param  out_z  z-velocity
/// \param  out_p  pressure
// ***************************************************************************************
void drift(Field &out_x, Field &out_y, Field &out_z, Field &out_p, const Settings::initial_conditions::drift &drift) {
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    const real u_lin = drift.velocity_lin[CoordinateAxis::X];
    const real v_lin = drift.velocity_lin[CoordinateAxis::Y];
    const real w_lin = drift.velocity_lin[CoordinateAxis::Z];
    const real pa = drift.pa;
    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x, out_y, out_z, out_p) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        out_x[idx] = u_lin;
        out_y[idx] = v_lin;
        out_z[idx] = w_lin;
        out_p[idx] = pa;
    }
}


// ================================ Diffusion Test - IC for u,v,w ========================
// ***************************************************************************************
/// \brief  Initial set up for Diffusion Test (c*exp*sin*sin*sin)
/// \param  out velocity
/// \param  t   time
// ***************************************************************************************
void exp_sinus_prod(Field &out, real t, const Settings::initial_conditions::exp_sinus_prod &esp) {
    auto domain_data = DomainData::getInstance();
    const real nu = domain_data->get_physical_parameters().nu.value();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    real A = 1.0;
    const real l = esp.l;
    real kpinu = 3 * l * l * M_PI * M_PI * nu;

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;

    //inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
        out[idx] = A * exp(-kpinu * t)
                * sin(l * M_PI * xi(coords_i, X1, dx))
                * sin(l * M_PI * yj(coords_j, Y1, dy))
                * sin(l * M_PI * zk(coords_k, Z1, dz));
    }
}

// ============================ Burgers Test - IC for u,v,w ==============================
// ***************************************************************************************
/// \brief  Initial set up for Burgers Test (c*exp*sin(x+y+z))
/// \param  out_x  x-velocity
/// \param  out_y  y-velocity
/// \param  out_z  z-velocity
/// \param  t   time
// ***************************************************************************************
void exp_sinus_sum(Field &out_x, Field &out_y, Field &out_z, real t) {
    auto domain_data = DomainData::getInstance();
    const real nu = domain_data->get_physical_parameters().nu.value();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();
    size_t Nz = domain_data->get_Nz();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;

    if (Nz > 3) {
        real d = 3.;                // 3D
        // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x, out_y, out_z) async
        for (size_t i = 0; i < size_domain_inner_list; i++) {
            size_t idx = domain_inner_list[i];
            coords_k = getCoordinateK(idx, Nx, Ny);
            coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
            coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);

            out_x[idx] = exp(-d * nu * t)
                    * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz));
            out_y[idx] = -0.5 * exp(-d * nu * t)
                    * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz));
            out_z[idx] = -0.5 * exp(-d * nu * t)
                    * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz));
        }
    } else {
        real d = 2.;                // 2D

        // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x, out_y, out_z) async
        for (size_t i = 0; i < size_domain_inner_list; i++) {
            size_t idx = domain_inner_list[i];
            coords_k = getCoordinateK(idx, Nx, Ny);
            coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
            coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);

            out_x[idx] = exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy));
            out_y[idx] = -exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy));
            out_z[idx] = 0.;
        }
    }
}

// ============================= Diffusion Test - IC for u,v,w ===========================
// ***************************************************************************************
/// \brief  Initial set up for Diffusion Test (c*sin*sin*sin)
/// \param  out velocity
// ***************************************************************************************
void fac_sin_sin_sin(Field &out, const Settings::initial_conditions::sin_sin_sin &sin_sin_sin) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    const real l = sin_sin_sin.l;
    real dkpi = 3 * l * l * M_PI * M_PI;
    real rdkpi = 1. / dkpi;

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
        out[idx] = -rdkpi
                * sin(l * M_PI * xi(coords_i, X1, dx))
                * sin(l * M_PI * yj(coords_j, Y1, dy))
                * sin(l * M_PI * zk(coords_k, Z1, dz));
    }
}

// ============================= Advection Test - IC for u,v,w ===========================
// ***************************************************************************************
/// \brief  Initial set up for Advection Test
/// \param  out velocity
/// \param  t time
// ***************************************************************************************
void gauss_bubble(Field &out, real t, const Settings::initial_conditions::gauss_bubble &gauss_bubble) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;

    const real u_lin = gauss_bubble.velocity_lin[CoordinateAxis::X];
    const real v_lin = gauss_bubble.velocity_lin[CoordinateAxis::Y];
    const real w_lin = gauss_bubble.velocity_lin[CoordinateAxis::Z];
    const real x_shift = gauss_bubble.shift[CoordinateAxis::X];
    const real y_shift = gauss_bubble.shift[CoordinateAxis::Y];
    const real z_shift = gauss_bubble.shift[CoordinateAxis::Z];
    const real l = gauss_bubble.l;
    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);

        real x_shift2 = ((xi(coords_i, X1, dx) - x_shift) / u_lin - t) * ((xi(coords_i, X1, dx) - x_shift) / u_lin - t);
        real y_shift2 = ((yj(coords_j, Y1, dy) - y_shift) / v_lin - t) * ((yj(coords_j, Y1, dy) - y_shift) / v_lin - t);
        real z_shift2 = 0;
        if (w_lin != 0) {
            z_shift2 = ((zk(coords_k, Z1, dz) - z_shift) / w_lin - t) * ((zk(coords_k, Z1, dz) - z_shift) / w_lin - t);
        }
        real quot = 1. / (2. * l * l);

        out[idx] = exp(-(x_shift2 + y_shift2 + z_shift2) * quot);
    }
}

// ======================== Layers (e.g. for temperature in PIV experiments) =============
// ***************************************************************************************
/// \brief  Initial set up as layers throughout the domain
/// \param  out temperature
// ***************************************************************************************
void layers(Field &out, const Settings::initial_conditions::layers_temperature &layers) {
    auto domain_data = DomainData::getInstance();

    std::vector<real> borders;
    borders.resize(layers.number_of_layers + 1);
    std::copy(layers.borders.begin(), layers.borders.end(), borders.begin() + 1);
    // layer border
    borders[0] = domain_data->get_start_coord_CD(layers.dir);
    borders[layers.number_of_layers] = domain_data->get_end_coord_CD(layers.dir);

    // set values into layers
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    auto domain_controller = DomainController::getInstance();
    size_t *domain_data_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;
    real x, y, z;

    //TODO highly inefficient
    if (layers.dir == X) {
        for (int l = 0; l < layers.number_of_layers; ++l) {
            // inner cells
#pragma acc parallel loop independent present(domain_data_list[:size_domain_inner_list], out) copyin(borders[:layers.number_of_layers+1], layers.values[:layers.number_of_layers]) async
            for (size_t i = 0; i < size_domain_inner_list; i++) {
                size_t idx = domain_data_list[i];
                coords_k = getCoordinateK(idx, Nx, Ny);
                coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
                coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
                x = xi(coords_i, X1, dx) - 0.5 * dx;
                if (borders[l] <= x && x <= borders[l + 1]) {
                    out[idx] = layers.values[l];
                }
            }
        }
    } else if (layers.dir == Y) {
        for (int l = 0; l < layers.number_of_layers; ++l) {
#pragma acc parallel loop independent present(domain_data_list[:size_domain_inner_list], out) copyin(borders[:layers.number_of_layers+1], layers.values[:layers.number_of_layers]) async
            for (size_t i = 0; i < size_domain_inner_list; i++) {
                size_t idx = domain_data_list[i];
                coords_k = getCoordinateK(idx, Nx, Ny);
                coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
                coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
                y = yj(coords_j, Y1, dy) - 0.5 * dy;
                if (borders[l] <= y && y <= borders[l + 1]) {
                    out[idx] = layers.values[l];
                }
            }
        }
    } else if (layers.dir == Z) {
        for (int l = 0; l < layers.number_of_layers; ++l) {
#pragma acc parallel loop independent present(domain_data_list[:size_domain_inner_list], out) copyin(borders[:layers.number_of_layers+1], layers.values[:layers.number_of_layers]) async
            for (size_t i = 0; i < size_domain_inner_list; i++) {
                size_t idx = domain_data_list[i];
                coords_k = getCoordinateK(idx, Nx, Ny);
                coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
                coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
                z = zk(coords_k, Z1, dz) - 0.5 * dz;
                if (borders[l] <= z && z <= borders[l + 1]) {
                    out[idx] = layers.values[l];
                }
            }
        }
    }
}


// ============================= Diffusion Test - IC for u,v,w ===========================
// ***************************************************************************************
/// \brief  Initial set up for Diffusion Test
/// \param  out velocity
// ***************************************************************************************
void hat(Field &out, const Settings::initial_conditions::hat &hat) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;

    const real start_x = hat.start_coords[CoordinateAxis::X];
    const real start_y = hat.start_coords[CoordinateAxis::Y];
    const real start_z = hat.start_coords[CoordinateAxis::Z];
    const real end_x = hat.end_coords[CoordinateAxis::X];
    const real end_y = hat.end_coords[CoordinateAxis::Y];
    const real end_z = hat.end_coords[CoordinateAxis::Z];
    const real val_in = hat.val_in;
    const real val_out = hat.val_out;
    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);

        if ((start_x <= xi(coords_i, X1, dx) && xi(coords_i, X1, dx) <= end_x) &&
            (start_y <= yj(coords_j, Y1, dy) && yj(coords_j, Y1, dy) <= end_y) &&
            (start_z <= zk(coords_k, Z1, dz) && zk(coords_k, Z1, dz) <= end_z)) {
            out[idx] = val_in;
        } else {
            out[idx] = val_out;
        }
    }
}

// ============================================== Jet ==============================================
// *************************************************************************************************
/// \brief Initial set up of a jet stream in a defined pipe
/// \param out velocity field
/// \param index_x1 starting index in x-direction
/// \param index_x2 ending index in x-direction
/// \param index_y1 starting index in y-direction
/// \param index_y2 ending index in y-direction
/// \param index_z1 starting index in z-direction
/// \param index_z2 ending index in z-direction
/// \param value velocity value to be set
// *************************************************************************************************
void jet(Field &out, const Settings::initial_conditions::jet &jet) {
    auto domain_data = DomainData::getInstance();

    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();
    Coordinate<size_t> start;
    for (size_t axis = 0; axis <= number_of_axes; axis++) {
        auto c_axis = CoordinateAxis(axis);
        start[axis]  = Utility::get_index(jet.start_coords[axis], domain_data->get_spacing(c_axis), domain_data->get_start_coord_PD(c_axis));
    }
    start[jet.dir] = domain_data->get_start_index_CD(jet.dir);
    Coordinate<size_t> end;
    for (size_t axis = 0; axis <= number_of_axes; axis++) {
        auto c_axis = CoordinateAxis(axis);
        end[axis]  = Utility::get_index(jet.end_coords[axis], domain_data->get_spacing(c_axis), domain_data->get_end_coord_PD(c_axis));
    }
    end[jet.dir] = domain_data->get_end_index_CD(jet.dir);

#pragma acc parallel loop independent present(out) async
    for (size_t i = start[CoordinateAxis::X]; i <= end[CoordinateAxis::X]; i++) {
        for (size_t j = start[CoordinateAxis::Y]; j <= end[CoordinateAxis::Y]; j++) {
            for (size_t k = start[CoordinateAxis::Z]; k <= end[CoordinateAxis::Z]; k++) {
                size_t index = IX(i, j, k, Nx, Ny);
                out[index] = jet.value;
            }
        }
    }
}

// ========================== NS Test - McDermott IC for u,v,w,p =========================
// ***************************************************************************************
/// \brief  Initial set up for NS Test - McDermott
/// \param  out_x  x-velocity
/// \param  out_y  y-velocity
/// \param  out_z  z-velocity
/// \param  out_p  pressure
/// \param  t   time
// ***************************************************************************************
void mcdermott(Field &out_x, Field &out_y, Field &out_z, Field &out_p, real t, const Settings::initial_conditions::mc_dermott &mc_dermott) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_k, coords_i, coords_j;
    const real A = mc_dermott.A;
    const real nu = domain_data->get_physical_parameters().nu.value();

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x, out_y, out_z, out_p) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);

        out_x[idx] = 1. - A * cos(xi(coords_i, X1, dx) - t) * sin(yj(coords_j, Y1, dy) - t) * exp(-2 * nu * t);
        out_y[idx] = 1. + A * sin(xi(coords_i, X1, dx) - t) * cos(yj(coords_j, Y1, dy) - t) * exp(-2 * nu * t);
        out_z[idx] = 0.;
        out_p[idx] = -0.25 * A * A * (cos(2 * (xi(coords_i, X1, dx) - t)) + cos(2 * (yj(coords_j, Y1, dy) - t))) * exp(-4 * nu * t);
    }
}

// ======== Random Function - Superposition of field values with random values ===========
// ***************************************************************************************
/// \brief  Creates random absolute/relative noise on given field
/// \param  out          field
/// \param  range        range of random numbers
/// \param  is_absolute  check if random number is relative (multiply) or absolute (additive)
/// \param  seed         custom seed if given, else seed <= 0
/// \param  step_size    interval steps of random numbers
// ***************************************************************************************
void random(Field &out, const Settings::random_parameters &random_params) {
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    IRandomField *noise_maker;
    if (random_params.custom_seed) {
        noise_maker = new UniformRandom(random_params.range, random_params.step_size, random_params.seed);
    } else {
        noise_maker = new UniformRandom(random_params.range, random_params.step_size);
    }

    out.update_host();
    auto noise = noise_maker->random_field(out.get_size());

    // inner cells
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        if (random_params.absolute) {
            out[idx] += noise[idx];
        } else {
            out[idx] *= (1 + noise[idx]);
        }
    }
    out.update_dev();
}

// ================================= Pressure Test - IC for p ============================
// ***************************************************************************************
/// \brief  Initial set up for Pressure Test (sin*sin*sin)
/// \param  out   pressure
// ***************************************************************************************
void sin_sin_sin(Field &out, const Settings::initial_conditions::sin_sin_sin &sin_sin_sin) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_i, coords_j, coords_k;
    const real l = sin_sin_sin.l;

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);

        out[idx] = sin(l * M_PI * xi(coords_i, X1, dx))
                    * sin(l * M_PI * yj(coords_j, Y1, dy))
                    * sin(l * M_PI * zk(coords_k, Z1, dz));
    }
}

// ======================= uniform distribution (eg. for force) ==========================
// ***************************************************************************************
/// \brief  Initial uniform set up
/// \param  out   force
/// \param  val   value of uniform distribution
// ***************************************************************************************
void uniform(Field &out, const Settings::initial_conditions::uniform &uniform) {
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    // inner cells
    real val = uniform.value;
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        out[idx] = val;
    }
}

// ============================= NS Test - Vortex IC for u,v,w,p =========================
// ***************************************************************************************
/// \brief  Initial set up for NS Test - Vertex
/// \param  out_x    x-velocity
/// \param  out_y    y-velocity
/// \param  out_z    z-velocity
/// \param  out_p    pressure
// ***************************************************************************************
void vortex(Field &out_x, Field &out_y, Field &out_z, Field &out_p, const Settings::initial_conditions::vortex &vortex) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();

    const real u_lin = vortex.velocity_lin[CoordinateAxis::X];
    const real w_lin = vortex.velocity_lin[CoordinateAxis::Y];
    const real v_lin = vortex.velocity_lin[CoordinateAxis::Z];
    const real rhoa = vortex.rhoa;
    const real pa = vortex.pa;

    real L = domain_data->get_lx();
    real R_c = L / 20.;
    real G = 0.04 * u_lin * R_c * sqrt(exp(1));

    real GrR_c = G / (R_c * R_c);
    real rR_c = 1. / (2. * R_c * R_c);
    real rhoGrR_c = rhoa * G * G * rR_c;

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_k, coords_i, coords_j;

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x, out_y, out_z, out_p) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);

        out_x[idx] = u_lin - GrR_c * yj(coords_j, Y1, dy) * exp(-rR_c
                * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx)
                +  yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
        out_y[idx] = v_lin + GrR_c * xi(coords_i, X1, dx) * exp(-rR_c
                * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx)
                +  yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
        out_z[idx] = 0.;
        out_p[idx] = pa - rhoGrR_c * exp(-rR_c
                * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx)
                +  yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
    }
}

void vortex_y(Field &out_x, Field &out_y, Field &out_z, Field &out_p, const Settings::initial_conditions::vortex &vortex) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();

    const real u_lin = vortex.velocity_lin[CoordinateAxis::X];
    const real w_lin = vortex.velocity_lin[CoordinateAxis::Y];
    const real v_lin = vortex.velocity_lin[CoordinateAxis::Z];
    const real rhoa = vortex.rhoa;
    const real pa = vortex.pa;

    real L = domain_data->get_ly();
    real R_c = L / 20.;
    real G = 0.04 * u_lin * R_c * sqrt(exp(1));

    real GrR_c = G / (R_c * R_c);
    real rR_c = 1. / (2. * R_c * R_c);
    real rhoGrR_c = rhoa * G * G * rR_c;

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);
    size_t coords_k, coords_i, coords_j;

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out_x, out_y, out_z, out_p) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);

        out_x[idx] = u_lin - GrR_c * yj(coords_j, Y1, dy) * exp(-rR_c
                * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx)
                +  yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
        out_y[idx] = v_lin + GrR_c * xi(coords_i, X1, dx) * exp(-rR_c
                * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx)
                +  yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
        out_z[idx] = 0.;
        out_p[idx] = pa - rhoGrR_c * exp(-rR_c
                * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx)
                +  yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
    }
}
}  // namespace Functions
