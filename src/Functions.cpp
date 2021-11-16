/// \file       Functions.cpp
/// \brief      Functions for Initialization
/// \date       Jun 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>

#ifdef _OPENACC
#include <accelmath.h>
#endif
#include <ctime>
#include <random>

#include "Functions.h"
#include "utility/Parameters.h"
#include "DomainData.h"
#include "utility/Utility.h"
#include "boundary/BoundaryController.h"

const std::string FunctionNames::Beltrami = "Beltrami";
const std::string FunctionNames::BuoyancyMMS = "BuoyancyMMS";
const std::string FunctionNames::BuoyancyST_MMS = "BuoyancyST_MMS";
const std::string FunctionNames::Drift = "Drift";
const std::string FunctionNames::ExpSinusProd = "ExpSinusProd";
const std::string FunctionNames::ExpSinusSum = "ExpSinusSum";
const std::string FunctionNames::GaussBubble = "GaussBubble";
const std::string FunctionNames::Hat = "Hat";
const std::string FunctionNames::Jet = "Jet";
const std::string FunctionNames::McDermott = "McDermott";
const std::string FunctionNames::RandomC = "RandomC";
const std::string FunctionNames::SinSinSin = "SinSinSin";
const std::string FunctionNames::Uniform = "Uniform";
const std::string FunctionNames::Vortex = "Vortex";
const std::string FunctionNames::VortexY = "VortexY";
const std::string FunctionNames::Zero = "Zero";

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
    void Beltrami(Field &out_x, Field &out_y, Field &out_z, Field &out_p, real t) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto params = Parameters::getInstance();

        real a = params->get_real("initial_conditions/a");  // 0.25 * M_PI;
        real d = params->get_real("initial_conditions/d");  // 0.5 * M_PI;
        real nu = params->get_real("physical_parameters/nu");  // 1;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x, out_y, out_z, out_p) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void BeltramiBC_p(Field &out_x) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto params = Parameters::getInstance();

        real a = params->get_real("initial_conditions/a");  // 0.25 * M_PI;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void BeltramiBC_u(Field &out_x, real t) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto params = Parameters::getInstance();

        real a = params->get_real("initial_conditions/a");  // 0.25 * M_PI;
        real d = params->get_real("initial_conditions/d");  // 0.5 * M_PI;
        real nu = params->get_real("physical_parameters/nu");  // 1.;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_list_level_joined();
        size_t size_domain_list = boundary->get_slice_size_domain_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void BeltramiBC_v(Field &out_x, real t) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto params = Parameters::getInstance();

        real a = params->get_real("initial_conditions/a");  // 0.25 * M_PI;
        real d = params->get_real("initial_conditions/d");  // 0.5 * M_PI;
        real nu = params->get_real("physical_parameters/nu");  // 1.;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_list_level_joined();
        size_t size_domain_list = boundary->get_slice_size_domain_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void BeltramiBC_w(Field &out_x, real t) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto params = Parameters::getInstance();

        real a = params->get_real("initial_conditions/a");  // 0.25 * M_PI;
        real d = params->get_real("initial_conditions/d");  // 0.25 * M_PI;
        real nu = params->get_real("physical_parameters/nu");

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_list_level_joined();
        size_t size_domain_list = boundary->get_slice_size_domain_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void BuoyancyForce(Field &out, Field &T, Field &T_ambient) {
        auto params = Parameters::getInstance();
        real beta = params->get_real("physical_parameters/beta");
        real g = params->get_real("physical_parameters/g");  // -9.81;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out, T, T_ambient) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void BuoyancyMMS(Field &out_x, Field &out_y, Field &out_z, Field &out_p, Field &out_T, real t) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();

        auto params = Parameters::getInstance();
        real nu = params->get_real("physical_parameters/nu");
        real beta = params->get_real("physical_parameters/beta");
        real g = params->get_real("physical_parameters/g");
        real rhoa = params->get_real("initial_conditions/rhoa");
        real rbeta = 1. / beta;
        real rg = 1. / g;
        real c = 2 * nu * M_PI * M_PI - 1;
        real rpi = 1. / M_PI;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_k, coords_i, coords_j;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x, out_y, out_z, out_p, out_T) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void BuoyancyST_MMS(Field &out, real t) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();

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
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);

        size_t coords_k, coords_i, coords_j;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void Drift(Field &out_x, Field &out_y, Field &out_z, Field &out_p) {
        auto params = Parameters::getInstance();

        real u_lin = params->get_real("initial_conditions/u_lin");
        real v_lin = params->get_real("initial_conditions/v_lin");
        real w_lin = params->get_real("initial_conditions/w_lin");
        real pa = params->get_real("initial_conditions/pa");

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x, out_y, out_z, out_p) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void ExpSinusProd(Field &out, real t) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto params = Parameters::getInstance();

        real nu = params->get_real("physical_parameters/nu");
        real l = params->get_real("initial_conditions/l");
        real A = 1.0;

        real kpinu = 3 * l * l * M_PI * M_PI * nu;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;

        //inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void ExpSinusSum(Field &out_x, Field &out_y, Field &out_z, real t) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();
        size_t Nz = domain->get_Nz();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto params = Parameters::getInstance();

        real nu = params->get_real("physical_parameters/nu");

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;

        if (Nz > 3) {
            real d = 3.;                // 3D
            // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x, out_y, out_z) async
            for (size_t i = 0; i < size_domain_list; i++) {
                size_t idx = domain_list[i];
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
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x, out_y, out_z) async
            for (size_t i = 0; i < size_domain_list; i++) {
                size_t idx = domain_list[i];
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
    void FacSinSinSin(Field &out) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto params = Parameters::getInstance();
        real l = params->get_real("initial_conditions/l"); //2;

        real dkpi = 3 * l * l * M_PI * M_PI;
        real rdkpi = 1. / dkpi;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void GaussBubble(Field &out, real t) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto params = Parameters::getInstance();
        real u_lin = params->get_real("initial_conditions/u_lin");
        real v_lin = params->get_real("initial_conditions/v_lin");
        real w_lin = params->get_real("initial_conditions/w_lin");
        real x_shift = params->get_real("initial_conditions/x_shift");
        real y_shift = params->get_real("initial_conditions/y_shift");
        real z_shift = params->get_real("initial_conditions/z_shift");
        real l = params->get_real("initial_conditions/l");

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
            coords_k = getCoordinateK(idx, Nx, Ny);
            coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
            coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);

            real x_shift2 = ((xi(coords_i, X1, dx) - x_shift) / u_lin - t) * ((xi(coords_i, X1, dx) - x_shift) / u_lin - t);
            real y_shift2 = ((yj(coords_j, Y1, dy) - y_shift) / v_lin - t) * ((yj(coords_j, Y1, dy) - y_shift) / v_lin - t);
            real z_shift2 = ((zk(coords_k, Z1, dz) - z_shift) / w_lin - t) * ((zk(coords_k, Z1, dz) - z_shift) / w_lin - t);
            real quot = 1. / (2. * l * l);

            out[idx] = exp(-(x_shift2 + y_shift2 + z_shift2) * quot);
        }
    }

// ======================== Layers (e.g. for temperature in PIV experiments) =============
// ***************************************************************************************
/// \brief  Initial set up as layers throughout the domain
/// \param  out temperature
// ***************************************************************************************
    void Layers(Field &out) {
        auto domain = DomainData::getInstance();
        auto params = Parameters::getInstance();
        int n_layers = params->get_int("initial_conditions/n_layers");

        // layer border
        real *bord = new real[n_layers + 1];
        real val_bord;

        for (int l = 1; l < n_layers; ++l) {
            std::string val_bord_l = "initial_conditions/border_";
            val_bord_l += std::to_string(l);
            val_bord = params->get_real(val_bord_l);
            bord[l] = val_bord;
        }

        std::string dir = params->get("initial_conditions/dir"); //x,y,z

        if (dir == "x") {
            real x1 = domain->get_x1();
            real x2 = domain->get_x2();
            bord[0] = x1;
            bord[n_layers] = x2;
        } else if (dir == "y") {
            real y1 = domain->get_y1();
            real y2 = domain->get_y2();
            bord[0] = y1;
            bord[n_layers] = y2;
        } else if (dir == "z") {
            real z1 = domain->get_z1();
            real z2 = domain->get_z2();
            bord[0] = z1;
            bord[n_layers] = z2;
        } else {
#ifndef BENCHMARKING
            auto m_logger = Utility::create_logger("Functions");
            m_logger->error("No distance for layers specified!");
#endif
            //TODO(issue 6) Error handling
        }

        // get values in layers
        // layer values
        real *val = new real[n_layers];
        real val_out;

        for (int l = 0; l < n_layers; ++l) {
            std::string val_out_l = "initial_conditions/value_";
            val_out_l += std::to_string(l + 1);
            val_out = params->get_real(val_out_l);
            val[l] = val_out;
        }

        // set values into layers
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;
        real x, y, z;

        if (dir == "x") {
            for (int l = 0; l < n_layers; ++l) {
                // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out) copyin(bord[:n_layers+1], val[:n_layers]) async
                for (size_t i = 0; i < size_domain_list; i++) {
                    size_t idx = domain_list[i];
                    coords_k = getCoordinateK(idx, Nx, Ny);
                    coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
                    coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
                    x = xi(coords_i, X1, dx) - 0.5 * dx;
                    if (bord[l] <= x && x <= bord[l + 1]) {
                        out[idx] = val[l];
                    }
                }
            }

        } else if (dir == "y") {
            for (int l = 0; l < n_layers; ++l) {
                // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out) copyin(bord[:n_layers+1], val[:n_layers]) async
                for (size_t i = 0; i < size_domain_list; i++) {
                    size_t idx = domain_list[i];
                    coords_k = getCoordinateK(idx, Nx, Ny);
                    coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
                    coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
                    y = yj(coords_j, Y1, dy) - 0.5 * dy;
                    if (bord[l] <= y && y <= bord[l + 1]) {
                        out[idx] = val[l];
                    }
                }
            }

        } else if (dir == "z") {
            for (int l = 0; l < n_layers; ++l) {
                // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out) copyin(bord[:n_layers+1], val[:n_layers]) async
                for (size_t i = 0; i < size_domain_list; i++) {
                    size_t idx = domain_list[i];
                    coords_k = getCoordinateK(idx, Nx, Ny);
                    coords_j = getCoordinateJ(idx, Nx, Ny, coords_k);
                    coords_i = getCoordinateI(idx, Nx, Ny, coords_j, coords_k);
                    z = zk(coords_k, Z1, dz) - 0.5 * dz;
                    if (bord[l] <= z && z <= bord[l + 1]) {
                        out[idx] = val[l];
                    }
                }
            }

        } else {
#ifndef BENCHMARKING
            auto m_logger = Utility::create_logger("Functions");
            m_logger->error("No distance for layers specified!");
#endif
            //TODO(issue 6) Error handling
        }
    }


// ============================= Diffusion Test - IC for u,v,w ===========================
// ***************************************************************************************
/// \brief  Initial set up for Diffusion Test
/// \param  out velocity
// ***************************************************************************************
    void Hat(Field &out) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto params = Parameters::getInstance();
        real start_x = params->get_real("initial_conditions/x1");
        real end_x = params->get_real("initial_conditions/x2");
        real start_y = params->get_real("initial_conditions/y1");
        real end_y = params->get_real("initial_conditions/y2");
        real start_z = params->get_real("initial_conditions/z1");
        real end_z = params->get_real("initial_conditions/z2");
        real val_in = params->get_real("initial_conditions/val_in");
        real val_out = params->get_real("initial_conditions/val_out");

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void Jet(
            Field &out,
            const size_t index_x1, const size_t index_x2,
            const size_t index_y1, const size_t index_y2,
            const size_t index_z1, const size_t index_z2,
            real value) {
        auto domain = DomainData::getInstance();

        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

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
    void McDermott(Field &out_x, Field &out_y, Field &out_z, Field &out_p, real t) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();

        auto params = Parameters::getInstance();
        real nu = params->get_real("physical_parameters/nu");

        real A = params->get_real("initial_conditions/A"); //2;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_k, coords_i, coords_j;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x, out_y, out_z, out_p) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void Random(Field &out, real range, bool is_absolute, int seed, real step_size) {
        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);

        std::mt19937 mt;
        int steps = static_cast<int>(range / step_size);
        if (seed > 0) {
            mt = std::mt19937(seed);
        } else {
            std::random_device rd;
            mt = std::mt19937(rd());
        }
        std::uniform_int_distribution<int> dist(-steps, steps);

        // inner cells
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void SinSinSin(Field &out) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();
        real dz = domain->get_dz();

        auto params = Parameters::getInstance();
        real l = params->get_real("initial_conditions/l"); //2;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_i, coords_j, coords_k;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void Uniform(Field &out, real val) {
        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_list_level_joined();
        size_t size_domain_list = boundary->get_slice_size_domain_list_level_joined(0);

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
    void Vortex(Field &out_x, Field &out_y, Field &out_z, Field &out_p) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();

        auto params = Parameters::getInstance();

        real u_lin = params->get_real("initial_conditions/u_lin");
        real v_lin = params->get_real("initial_conditions/v_lin");

        real L = domain->get_lx();
        real R_c = L / 20.;
        real G = 0.04 * u_lin * R_c * sqrt(exp(1));
        real pa = params->get_real("initial_conditions/pa");
        real rhoa = params->get_real("initial_conditions/rhoa");

        real GrR_c = G / (R_c * R_c);
        real rR_c = 1. / (2. * R_c * R_c);
        real rhoGrR_c = rhoa * G * G * rR_c;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_k, coords_i, coords_j;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x, out_y, out_z, out_p) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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

    void VortexY(Field &out_x, Field &out_y, Field &out_z, Field &out_p) {
        auto domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();

        real dx = domain->get_dx();
        real dy = domain->get_dy();

        auto params = Parameters::getInstance();

        real u_lin = params->get_real("initial_conditions/u_lin");
        real v_lin = params->get_real("initial_conditions/v_lin");

        real L = domain->get_ly();
        real R_c = L / 20.;
        real G = 0.04 * u_lin * R_c * sqrt(exp(1));
        real pa = params->get_real("initial_conditions/pa");
        real rhoa = params->get_real("initial_conditions/rhoa");

        real GrR_c = G / (R_c * R_c);
        real rR_c = 1. / (2. * R_c * R_c);
        real rhoGrR_c = rhoa * G * G * rR_c;

        auto boundary = BoundaryController::getInstance();
        size_t *domain_list = boundary->get_domain_inner_list_level_joined();
        size_t size_domain_list = boundary->get_size_domain_inner_list_level_joined(0);
        size_t coords_k, coords_i, coords_j;

        // inner cells
#pragma acc parallel loop independent present(domain_list[:size_domain_list], out_x, out_y, out_z, out_p) async
        for (size_t i = 0; i < size_domain_list; i++) {
            size_t idx = domain_list[i];
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
