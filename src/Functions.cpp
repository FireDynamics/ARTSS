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

const std::string FunctionNames::beltrami = "Beltrami";
const std::string FunctionNames::buoyancy_mms = "BuoyancyMMS";
const std::string FunctionNames::buoyancy_st_mms = "BuoyancyST_MMS";
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
void beltrami(Field &out_x, Field &out_y, Field &out_z, Field &out_p, real t, real a, real d, real nu) {
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

// ===================================== Buoyancy Force ==================================
// ***************************************************************************************
/// \brief  Buoyancy Force
/// \param  out   force
/// \param  T   Temperature
/// \param  T_ambient    Ambient temperature
// ***************************************************************************************
void buoyancy_force(Field &out, Field &T, Field &T_ambient, real beta, real g) {
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out, T, T_ambient) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        out[idx] = -beta * (T[idx] - T_ambient[idx]) * g;
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
void buoyancy_mms(Field &out_x, Field &out_y, Field &out_z, Field &out_p, Field &out_T, real t, real nu, real beta, real g, real rhoa) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();

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

// ========== NSTemp Test - MMS source term for temperature equation with buoyancy ========
// ***************************************************************************************
/// \brief  Source term for NSTemp Test - MMS with buoyant force
/// \param  out force
/// \param  t time
// ***************************************************************************************
void buoyancy_st_mms(Field &out, real t, real nu, real beta, real kappa, real g, real rhoa) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();

    real rbeta = 1. / beta;
    real rg = 1. / g;
    real c_nu = 2 * nu * M_PI * M_PI - 1;
    real c_kappa = 2 * kappa * M_PI * M_PI - 1;

    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    size_t coords_k, coords_i, coords_j;

    // inner cells
#pragma acc parallel loop independent present(domain_inner_list[:size_domain_inner_list], out) async
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        coords_k = getCoordinateK(idx, Nx, Ny);
        coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
        coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
        out[idx] = rhoa * rbeta * rg * 2 * c_nu * c_kappa * exp(-t)
                * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
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
void drift(Field &out_x, Field &out_y, Field &out_z, Field &out_p, real u_lin, real v_lin, real w_lin, real pa) {
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

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
void exp_sinus_prod(Field &out, real t, real nu, real l) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    real A = 1.0;

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
void exp_sinus_sum(Field &out_x, Field &out_y, Field &out_z, real t, real nu) {
    auto domain_data = DomainData::getInstance();
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
void fac_sin_sin_sin(Field &out, real l) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

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
void gauss_bubble(Field &out, real t,
        real u_lin, real v_lin, real w_lin,
        real x_shift, real y_shift, real z_shift,
        real l) {
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
void layers(Field &out, int n_layers, const CoordinateAxis axis, real *borders, const real *values) {
    auto domain_data = DomainData::getInstance();

    // layer border
    borders[0] = domain_data->get_start_coord_CD(axis);
    borders[n_layers] = domain_data->get_end_coord_CD(axis);

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
    if (axis == X) {
        for (int l = 0; l < n_layers; ++l) {
            // inner cells
#pragma acc parallel loop independent present(domain_data_list[:size_domain_inner_list], out) copyin(borders[:n_layers+1], values[:n_layers]) async
            for (size_t i = 0; i < size_domain_inner_list; i++) {
                size_t idx = domain_data_list[i];
                coords_k = getCoordinateK(idx, Nx, Ny);
                coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
                coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
                x = xi(coords_i, X1, dx) - 0.5 * dx;
                if (borders[l] <= x && x <= borders[l + 1]) {
                    out[idx] = values[l];
                }
            }
        }
    } else if (axis == Y) {
        for (int l = 0; l < n_layers; ++l) {
#pragma acc parallel loop independent present(domain_data_list[:size_domain_inner_list], out) copyin(borders[:n_layers+1], values[:n_layers]) async
            for (size_t i = 0; i < size_domain_inner_list; i++) {
                size_t idx = domain_data_list[i];
                coords_k = getCoordinateK(idx, Nx, Ny);
                coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
                coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
                y = yj(coords_j, Y1, dy) - 0.5 * dy;
                if (borders[l] <= y && y <= borders[l + 1]) {
                    out[idx] = values[l];
                }
            }
        }
    } else if (axis == Z) {
        for (int l = 0; l < n_layers; ++l) {
#pragma acc parallel loop independent present(domain_data_list[:size_domain_inner_list], out) copyin(borders[:n_layers+1], values[:n_layers]) async
            for (size_t i = 0; i < size_domain_inner_list; i++) {
                size_t idx = domain_data_list[i];
                coords_k = getCoordinateK(idx, Nx, Ny);
                coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
                coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
                z = zk(coords_k, Z1, dz) - 0.5 * dz;
                if (borders[l] <= z && z <= borders[l + 1]) {
                    out[idx] = values[l];
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
void hat(Field &out,
        real start_x, real end_x,
        real start_y, real end_y,
        real start_z, real end_z,
        real val_in, real val_out) {
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
void jet(
        Field &out,
        const size_t index_x1, const size_t index_x2,
        const size_t index_y1, const size_t index_y2,
        const size_t index_z1, const size_t index_z2,
        real value) {
    auto domain_data = DomainData::getInstance();

    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

#pragma acc parallel loop independent present(out) async
    for (size_t i = index_x1; i <= index_x2; i++) {
        for (size_t j = index_y1; j <= index_y2; j++) {
            for (size_t k = index_z1; k <= index_z2; k++) {
                size_t index = IX(i, j, k, Nx, Ny);
                out[index] = value;
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
void mcdermott(Field &out_x, Field &out_y, Field &out_z, Field &out_p, real t, real nu, real A) {
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
void random(Field &out, real range, bool is_absolute, int seed, real step_size) {
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    std::mt19937 mt;
    int steps = static_cast<int>(range / step_size);
    if (seed > 0) {
        mt = std::mt19937(seed);
    } else {
        std::random_device rd;
        mt = std::mt19937(rd());
    }
    std::uniform_int_distribution<int> dist(-steps, steps);

    out.update_host();
    // inner cells
    for (size_t i = 0; i < size_domain_inner_list; i++) {
        size_t idx = domain_inner_list[i];
        // generate secret number between -range and range:
        double no = dist(mt) * step_size;
        if (is_absolute) {
            out[idx] += (no);
        } else {
            out[idx] *= (1 + no);
        }
    }
    out.update_dev();
}

// ================================= Pressure Test - IC for p ============================
// ***************************************************************************************
/// \brief  Initial set up for Pressure Test (sin*sin*sin)
/// \param  out   pressure
// ***************************************************************************************
void sin_sin_sin(Field &out, real l) {
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
void uniform(Field &out, real val) {
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    // inner cells
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
void vortex(Field &out_x, Field &out_y, Field &out_z, Field &out_p, real u_lin, real v_lin, real pa, real rhoa) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();


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

void vortex_y(Field &out_x, Field &out_y, Field &out_z, Field &out_p, real u_lin, real v_lin, real pa, real rhoa) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();

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
