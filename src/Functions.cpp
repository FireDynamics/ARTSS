/// \file       Functions.cpp
/// \brief      Functions for Initialization
/// \date       Jun 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>
#include <ctime>
#include <random>

#include "Functions.h"
#include "utility/Parameters.h"
#include "Domain.h"
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
const std::string FunctionNames::McDermott = "McDermott";
const std::string FunctionNames::RandomC = "RandomC";
const std::string FunctionNames::RampTanh = "RampTanh";
const std::string FunctionNames::SinSinSin = "SinSinSin";
const std::string FunctionNames::Uniform = "Uniform";
const std::string FunctionNames::Vortex = "Vortex";
const std::string FunctionNames::VortexY = "VortexY";
const std::string FunctionNames::Zero = "Zero";

namespace Functions {

// ================================ NS Test - Beltrami IC =================================
// ***************************************************************************************
/// \brief  Initial set up for NS Test - Beltrami
/// \param  outx  x-velocity
/// \param  outy  y-velocity
/// \param  outz  z-velocity
/// \param  outp  pressure
/// \param  t time
// ***************************************************************************************
    void Beltrami(Field *outx, Field *outy, Field *outz, Field *outp, real t) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(outx->GetLevel());
        size_t Ny = domain->get_Ny(outx->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(outx->GetLevel());
        real dy = domain->get_dy(outx->GetLevel());
        real dz = domain->get_dz(outx->GetLevel());

        auto params = Parameters::getInstance();

        real a = params->get_real("initial_conditions/a"); //0.25 * M_PI;
        real d = params->get_real("initial_conditions/d"); //0.5 * M_PI;
        real nu = params->get_real("physical_parameters/nu"); //1;

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            outx->data[idx] = -a * (exp(a * xi(coords_i, X1, dx)) * sin(a * yj(coords_j, Y1, dy) + dz) +
                                    exp(a * zk(coords_k, Z1, dz)) * cos(a * xi(coords_i, X1, dx) + dy)) * exp(-nu * d * d * t);
            outy->data[idx] = -a * (exp(a * yj(coords_j, Y1, dy)) * sin(a * zk(coords_k, Z1, dz) + dx) +
                                    exp(a * xi(coords_i, X1, dx)) * cos(a * yj(coords_j, Y1, dy) + dz)) * exp(-nu * d * d * t);
            outz->data[idx] = -a * (exp(a * zk(coords_k, Z1, dz)) * sin(a * xi(coords_i, X1, dx) + dy) +
                                    exp(a * yj(coords_j, Y1, dy)) * cos(a * zk(coords_k, Z1, dz) + dx)) * exp(-nu * d * d * t);
            outp->data[idx] =
                    -0.5 * a * a * (exp(2 * a * xi(coords_i, X1, dx)) + exp(2 * a * yj(coords_j, Y1, dy)) + exp(2 * a * zk(coords_k, Z1, dz)) \
 + 2 * sin(a * xi(coords_i, X1, dx) + dy) * cos(a * zk(coords_k, Z1, dz) + dx) * exp(a * (yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz))) \
 + 2 * sin(a * yj(coords_j, Y1, dy) + dz) * cos(a * xi(coords_i, X1, dx) + dy) * exp(a * (zk(coords_k, Z1, dz) + xi(coords_i, X1, dx))) \
 + 2 * sin(a * zk(coords_k, Z1, dz) + dx) * cos(a * yj(coords_j, Y1, dy) + dz) * exp(a * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy))));
        }

        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            outx->data[idx] = -a * (exp(a * xi(coords_i, X1, dx)) * sin(a * yj(coords_j, Y1, dy) + dz) +
                                    exp(a * zk(coords_k, Z1, dz)) * cos(a * xi(coords_i, X1, dx) + dy)) * exp(-nu * d * d * t);
            outy->data[idx] = -a * (exp(a * yj(coords_j, Y1, dy)) * sin(a * zk(coords_k, Z1, dz) + dx) +
                                    exp(a * xi(coords_i, X1, dx)) * cos(a * yj(coords_j, Y1, dy) + dz)) * exp(-nu * d * d * t);
            outz->data[idx] = -a * (exp(a * zk(coords_k, Z1, dz)) * sin(a * xi(coords_i, X1, dx) + dy) +
                                    exp(a * yj(coords_j, Y1, dy)) * cos(a * zk(coords_k, Z1, dz) + dx)) * exp(-nu * d * d * t);
            outp->data[idx] =
                    -0.5 * a * a * (exp(2 * a * xi(coords_i, X1, dx)) + exp(2 * a * yj(coords_j, Y1, dy)) + exp(2 * a * zk(coords_k, Z1, dz)) \
 + 2 * sin(a * xi(coords_i, X1, dx) + dy) * cos(a * zk(coords_k, Z1, dz) + dx) * exp(a * (yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz))) \
 + 2 * sin(a * yj(coords_j, Y1, dy) + dz) * cos(a * xi(coords_i, X1, dx) + dy) * exp(a * (zk(coords_k, Z1, dz) + xi(coords_i, X1, dx))) \
 + 2 * sin(a * zk(coords_k, Z1, dz) + dx) * cos(a * yj(coords_j, Y1, dy) + dz) * exp(a * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy))));
        }

    }

// ================================ NS Test - Beltrami IC for p ==========================
// ***************************************************************************************
/// \brief  Initial pressure set up for NS Test - Beltrami
/// \param  outx  pressure
// ***************************************************************************************
    void BeltramiBC_p(Field *outx) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(outx->GetLevel());
        size_t Ny = domain->get_Ny(outx->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(outx->GetLevel());
        real dy = domain->get_dy(outx->GetLevel());
        real dz = domain->get_dz(outx->GetLevel());

        auto params = Parameters::getInstance();

        real a = params->get_real("initial_conditions/a"); //0.25 * M_PI;

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            outx->data[idx] =
                    -0.5 * a * a * (exp(2 * a * xi(coords_i, X1, dx)) + exp(2 * a * yj(coords_j, Y1, dy)) + exp(2 * a * zk(coords_k, Z1, dz)) \
 + 2 * sin(a * xi(coords_i, X1, dx) + dy) * cos(a * zk(coords_k, Z1, dz) + dx) * exp(a * (yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz))) \
 + 2 * sin(a * yj(coords_j, Y1, dy) + dz) * cos(a * xi(coords_i, X1, dx) + dy) * exp(a * (zk(coords_k, Z1, dz) + xi(coords_i, X1, dx))) \
 + 2 * sin(a * zk(coords_k, Z1, dz) + dx) * cos(a * yj(coords_j, Y1, dy) + dz) * exp(a * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy))));
        }
        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            outx->data[idx] =
                    -0.5 * a * a * (exp(2 * a * xi(coords_i, X1, dx)) + exp(2 * a * yj(coords_j, Y1, dy)) + exp(2 * a * zk(coords_k, Z1, dz)) \
 + 2 * sin(a * xi(coords_i, X1, dx) + dy) * cos(a * zk(coords_k, Z1, dz) + dx) * exp(a * (yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz))) \
 + 2 * sin(a * yj(coords_j, Y1, dy) + dz) * cos(a * xi(coords_i, X1, dx) + dy) * exp(a * (zk(coords_k, Z1, dz) + xi(coords_i, X1, dx))) \
 + 2 * sin(a * zk(coords_k, Z1, dz) + dx) * cos(a * yj(coords_j, Y1, dy) + dz) * exp(a * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy))));
        }
    }

// ================================ NS Test - Beltrami IC for u ==========================
// ***************************************************************************************
/// \brief  Initial x-velocity set up for NS Test - Beltrami
/// \param  outx  x-velocity
/// \param  t time
// ***************************************************************************************
    void BeltramiBC_u(Field *outx, real t) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(outx->GetLevel());
        size_t Ny = domain->get_Ny(outx->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(outx->GetLevel());
        real dy = domain->get_dy(outx->GetLevel());
        real dz = domain->get_dz(outx->GetLevel());

        auto params = Parameters::getInstance();

        real a = params->get_real("initial_conditions/a"); //0.25 * M_PI;
        real d = params->get_real("initial_conditions/d"); //0.5 * M_PI;
        real nu = params->get_real("physical_parameters/nu"); //1.;

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            outx->data[idx] = -a * (exp(a * xi(coords_i, X1, dx)) * sin(a * yj(coords_j, Y1, dy) + dz) +
                                    exp(a * zk(coords_k, Z1, dz)) * cos(a * xi(coords_i, X1, dx) + dy)) * exp(-nu * d * d * t);
        }
        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            outx->data[idx] = -a * (exp(a * xi(coords_i, X1, dx)) * sin(a * yj(coords_j, Y1, dy) + dz) +
                                    exp(a * zk(coords_k, Z1, dz)) * cos(a * xi(coords_i, X1, dx) + dy)) * exp(-nu * d * d * t);
        }
    }

// ================================ NS Test - Beltrami IC for v ==========================
// ***************************************************************************************
/// \brief  Initial y-velocity set up for NS Test - Beltrami
/// \param  outy  y-velocity
/// \param  t time
// ***************************************************************************************
    void BeltramiBC_v(Field *outx, real t) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(outx->GetLevel());
        size_t Ny = domain->get_Ny(outx->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(outx->GetLevel());
        real dy = domain->get_dy(outx->GetLevel());
        real dz = domain->get_dz(outx->GetLevel());

        auto params = Parameters::getInstance();

        real a = params->get_real("initial_conditions/a");//0.25 * M_PI;
        real d = params->get_real("initial_conditions/d");//0.5 * M_PI;
        real nu = params->get_real("physical_parameters/nu"); //1.;

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            outx->data[idx] = -a * (exp(a * yj(coords_j, Y1, dy)) * sin(a * zk(coords_k, Z1, dz) + dx) +
                                    exp(a * xi(coords_i, X1, dx)) * cos(a * yj(coords_j, Y1, dy) + dz)) * exp(-nu * d * d * t);
        }
        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            outx->data[idx] = -a * (exp(a * yj(coords_j, Y1, dy)) * sin(a * zk(coords_k, Z1, dz) + dx) +
                                    exp(a * xi(coords_i, X1, dx)) * cos(a * yj(coords_j, Y1, dy) + dz)) * exp(-nu * d * d * t);
        }
    }

// ================================ NS Test - Beltrami IC for w ==========================
// ***************************************************************************************
/// \brief  Initial z-velocity set up for NS Test - Beltrami
/// \param  outz  z-velocity
/// \param  t time
// ***************************************************************************************
    void BeltramiBC_w(Field *outx, real t) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(outx->GetLevel());
        size_t Ny = domain->get_Ny(outx->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(outx->GetLevel());
        real dy = domain->get_dy(outx->GetLevel());
        real dz = domain->get_dz(outx->GetLevel());

        auto params = Parameters::getInstance();

        real a = params->get_real("initial_conditions/a");//0.25 * M_PI;
        real d = params->get_real("initial_conditions/d");//0.25 * M_PI;
        real nu = params->get_real("physical_parameters/nu");

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            outx->data[idx] = -a * (exp(a * zk(coords_k, Z1, dz)) * sin(a * xi(coords_i, X1, dx) + dy) +
                                    exp(a * yj(coords_j, Y1, dy)) * cos(a * zk(coords_k, Z1, dz) + dx)) * exp(-nu * d * d * t);
        }
        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            outx->data[idx] = -a * (exp(a * zk(coords_k, Z1, dz)) * sin(a * xi(coords_i, X1, dx) + dy) +
                                    exp(a * yj(coords_j, Y1, dy)) * cos(a * zk(coords_k, Z1, dz) + dx)) * exp(-nu * d * d * t);
        }
    }

// ===================================== Buoyancy Force ==================================
// ***************************************************************************************
/// \brief  Buoyancy Force
/// \param  out   force
/// \param  T   Temperature
/// \param  Ta    Ambient temperature
// ***************************************************************************************
    void BuoyancyForce(Field *out, Field *T, Field *Ta) {

        auto d_out = out->data;
        auto d_T = T->data;
        auto d_Ta = Ta->data;

        auto params = Parameters::getInstance();

        real beta = params->get_real("physical_parameters/beta");

        real g = params->get_real("physical_parameters/g"); //-9.81;

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            d_out[idx] = -beta * (d_T[idx] - d_Ta[idx]) * g;
        }

        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            d_out[idx] = -beta * (d_T[idx] - d_Ta[idx]) * g;
        }
    }

// ================== NSTemp Test - MMS IC for u,v,w,p,T with buoyancy ===================
// ***************************************************************************************
/// \brief  Initial set up for NSTemp Test - MMS with buoyant force
/// \param  outx  x-velocity
/// \param  outy  y-velocity
/// \param  outz  z-velocity
/// \param  outp  pressure
/// \param  outT  temperature
/// \param  t   time
// ***************************************************************************************
    void BuoyancyMMS(Field *outx, Field *outy, Field *outz, Field *outp, Field *outT, real t) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(outx->GetLevel());
        size_t Ny = domain->get_Ny(outx->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();

        real dx = domain->get_dx(outx->GetLevel());
        real dy = domain->get_dy(outx->GetLevel());

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
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            outx->data[idx] = exp(-t) * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
            outy->data[idx] = -exp(-t) * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
            outz->data[idx] = 0.;
            outp->data[idx] = rhoa * rpi * c * exp(-t) * cos(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
            outT->data[idx] = rhoa * rbeta * rg * 2 * c * exp(-t) * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
        }

        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            outx->data[idx] = exp(-t) * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
            outy->data[idx] = -exp(-t) * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
            outz->data[idx] = 0.;
            outp->data[idx] = rhoa * rpi * c * exp(-t) * cos(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
            outT->data[idx] = rhoa * rbeta * rg * 2 * c * exp(-t) * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
        }
    }

// ========== NSTemp Test - MMS source term for temperature equation with buoyancy ========
// ***************************************************************************************
/// \brief  Source term for NSTemp Test - MMS with buoyant force
/// \param  out force
/// \param  t time
// ***************************************************************************************
    void BuoyancyST_MMS(Field *out, real t) {

        auto domain = Domain::getInstance();
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
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();

        std::vector<size_t> coords;
        size_t coords_i, coords_j;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            out->data[idx] =
                    rhoa * rbeta * rg * 2 * c_nu * c_kappa * exp(-t) * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
        }

        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            out->data[idx] =
                    rhoa * rbeta * rg * 2 * c_nu * c_kappa * exp(-t) * sin(M_PI * (xi(coords_i, X1, dx) + yj(coords_j, Y1, dy)));
        }
    }

// ===================================== NS Test - IC for u,v,w,p ========================
// ***************************************************************************************
/// \brief  Initial set up for NS Test - Flow around cube or Channel flow with Drift
/// \param  outx  x-velocity
/// \param  outy  y-velocity
/// \param  outz  z-velocity
/// \param  outp  pressure
// ***************************************************************************************
    void Drift(Field *outx, Field *outy, Field *outz, Field *outp) {

        auto params = Parameters::getInstance();

        real u_lin = params->get_real("initial_conditions/u_lin");
        real v_lin = params->get_real("initial_conditions/v_lin");
        real w_lin = params->get_real("initial_conditions/w_lin");
        real pa = params->get_real("initial_conditions/pa");


        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            outx->data[idx] = u_lin;
            outy->data[idx] = v_lin;
            outz->data[idx] = w_lin;
            outp->data[idx] = pa;
        }

        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            outx->data[idx] = u_lin;
            outy->data[idx] = v_lin;
            outz->data[idx] = w_lin;
            outp->data[idx] = pa;
        }
    }


// ================================ Diffusion Test - IC for u,v,w ========================
// ***************************************************************************************
/// \brief  Initial set up for Diffusion Test (c*exp*sin*sin*sin)
/// \param  out velocity
/// \param  t   time
// ***************************************************************************************
    void ExpSinusProd(Field *out, real t) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(out->GetLevel());
        size_t Ny = domain->get_Ny(out->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(out->GetLevel());
        real dy = domain->get_dy(out->GetLevel());
        real dz = domain->get_dz(out->GetLevel());

        auto params = Parameters::getInstance();

        real nu = params->get_real("physical_parameters/nu");
        real l = params->get_real("initial_conditions/l");
        real A = 1.0;

        real kpinu = 3 * l * l * M_PI * M_PI * nu;

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;

        //inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            out->data[idx] = A * exp(-kpinu * t) * sin(l * M_PI * xi(coords_i, X1, dx)) * sin(l * M_PI * yj(coords_j, Y1, dy)) *
                             sin(l * M_PI * zk(coords_k, Z1, dz));
        }
        //boundary
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            out->data[idx] = A * exp(-kpinu * t) * sin(l * M_PI * xi(coords_i, X1, dx)) * sin(l * M_PI * yj(coords_j, Y1, dy)) *
                             sin(l * M_PI * zk(coords_k, Z1, dz));
        }
    }

// ============================ Burgers Test - IC for u,v,w ==============================
// ***************************************************************************************
/// \brief  Initial set up for Burgers Test (c*exp*sin(x+y+z))
/// \param  outx  x-velocity
/// \param  outy  y-velocity
/// \param  outz  z-velocity
/// \param  t   time
// ***************************************************************************************
    void ExpSinusSum(Field *outx, Field *outy, Field *outz, real t) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(outx->GetLevel());
        size_t Ny = domain->get_Ny(outx->GetLevel());
        size_t Nz = domain->get_Nz(outx->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(outx->GetLevel());
        real dy = domain->get_dy(outx->GetLevel());
        real dz = domain->get_dz(outx->GetLevel());

        auto params = Parameters::getInstance();

        real nu = params->get_real("physical_parameters/nu");

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;

        if (Nz != 3) {
            real d = 3.;                // 3D

            //inner cells
            for (size_t i = 0; i < size_iList; i++) {
                size_t idx = iList[i];
                coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
                coords_i = coords[0];
                coords_j = coords[1];
                coords_k = coords[2];

                outx->data[idx] = exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz));
                outy->data[idx] = -0.5 * exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz));
                outz->data[idx] = -0.5 * exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz));
            }
            //boundary
            for (size_t i = 0; i < size_bList; i++) {
                size_t idx = bList[i];
                coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
                coords_i = coords[0];
                coords_j = coords[1];
                coords_k = coords[2];

                outx->data[idx] = exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz));
                outy->data[idx] = -0.5 * exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz));
                outz->data[idx] = -0.5 * exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) + zk(coords_k, Z1, dz));
            }

        } else {
            real d = 2.;                // 2D

            //inner cells
            for (size_t i = 0; i < size_iList; i++) {
                size_t idx = iList[i];
                coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
                coords_i = coords[0];
                coords_j = coords[1];

                outx->data[idx] = exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy));
                outy->data[idx] = -exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy));
                outz->data[idx] = 0.;
            }
            //boundary
            for (size_t i = 0; i < size_bList; i++) {
                size_t idx = bList[i];
                coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
                coords_i = coords[0];
                coords_j = coords[1];

                outx->data[idx] = exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy));
                outy->data[idx] = -exp(-d * nu * t) * sin(xi(coords_i, X1, dx) + yj(coords_j, Y1, dy));
                outz->data[idx] = 0.;
            }
        }
    }

// ============================= Diffusion Test - IC for u,v,w ===========================
// ***************************************************************************************
/// \brief  Initial set up for Diffusion Test (c*sin*sin*sin)
/// \param  out velocity
// ***************************************************************************************
    void FacSinSinSin(Field *out) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(out->GetLevel());
        size_t Ny = domain->get_Ny(out->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(out->GetLevel());
        real dy = domain->get_dy(out->GetLevel());
        real dz = domain->get_dz(out->GetLevel());

        auto params = Parameters::getInstance();

        real l = params->get_real("initial_conditions/l"); //2;
        real dkpi = 3 * l * l * M_PI * M_PI;
        real rdkpi = 1. / dkpi;

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            out->data[idx] = -rdkpi * sin(l * M_PI * xi(coords_i, X1, dx)) * sin(l * M_PI * yj(coords_j, Y1, dy)) *
                             sin(l * M_PI * zk(coords_k, Z1, dz));
        }
        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            out->data[idx] = -rdkpi * sin(l * M_PI * xi(coords_i, X1, dx)) * sin(l * M_PI * yj(coords_j, Y1, dy)) *
                             sin(l * M_PI * zk(coords_k, Z1, dz));
        }
    }

// ============================= Advection Test - IC for u,v,w ===========================
// ***************************************************************************************
/// \brief  Initial set up for Advection Test
/// \param  out velocity
/// \param  t time
// ***************************************************************************************
    void GaussBubble(Field *out, real t) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(out->GetLevel());
        size_t Ny = domain->get_Ny(out->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(out->GetLevel());
        real dy = domain->get_dy(out->GetLevel());
        real dz = domain->get_dz(out->GetLevel());

        auto params = Parameters::getInstance();

        real u_lin = params->get_real("initial_conditions/u_lin");
        real v_lin = params->get_real("initial_conditions/v_lin");
        real w_lin = params->get_real("initial_conditions/w_lin");
        real xshift = params->get_real("initial_conditions/xshift");
        real yshift = params->get_real("initial_conditions/yshift");
        real zshift = params->get_real("initial_conditions/zshift");
        real l = params->get_real("initial_conditions/l");

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];

            real xshift2 = ((xi(coords_i, X1, dx) - xshift) / u_lin - t) * ((xi(coords_i, X1, dx) - xshift) / u_lin - t);
            real yshift2 = ((yj(coords_j, Y1, dy) - yshift) / v_lin - t) * ((yj(coords_j, Y1, dy) - yshift) / v_lin - t);
            real zshift2 = ((zk(coords_k, Z1, dz) - zshift) / w_lin - t) * ((zk(coords_k, Z1, dz) - zshift) / w_lin - t);
            real quot = 1. / (2. * l * l);

            out->data[idx] = exp(-(xshift2 + yshift2 + zshift2) * quot);
        }
        // boundary
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];

            real xshift2 = ((xi(coords_i, X1, dx) - xshift) / u_lin - t) * ((xi(coords_i, X1, dx) - xshift) / u_lin - t);
            real yshift2 = ((yj(coords_j, Y1, dy) - yshift) / v_lin - t) * ((yj(coords_j, Y1, dy) - yshift) / v_lin - t);
            real zshift2 = ((zk(coords_k, Z1, dz) - zshift) / w_lin - t) * ((zk(coords_k, Z1, dz) - zshift) / w_lin - t);
            real quot = 1. / (2. * l * l);

            out->data[idx] = exp(-(xshift2 + yshift2 + zshift2) * quot);
        }
    }

// ======================== Layers (e.g. for temperature in PIV experiments) =============
// ***************************************************************************************
/// \brief  Initial set up as layers throughout the domain
/// \param  out temperature
// ***************************************************************************************
    void Layers(Field *out) {

        auto domain = Domain::getInstance();
        auto params = Parameters::getInstance();
        size_t n_layers = static_cast<size_t> (params->get_int("initial_conditions/n_layers"));

        // layer border
        real *bord = new real[n_layers + 1];
        real val_bord;

        for (size_t l = 1; l < n_layers; ++l) {
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
        }
        //TODO Error handling

        // get values in layers
        // layer values
        real *val = new real[n_layers];
        real val_out;

        for (size_t l = 0; l < n_layers; ++l) {
            std::string val_out_l = "initial_conditions/value_";
            val_out_l += std::to_string(l + 1);
            val_out = params->get_real(val_out_l);
            val[l] = val_out;
        }

        //set values into layers
        size_t Nx = domain->get_Nx(out->GetLevel());
        size_t Ny = domain->get_Ny(out->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(out->GetLevel());
        real dy = domain->get_dy(out->GetLevel());
        real dz = domain->get_dz(out->GetLevel());

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        size_t *oList = boundary->get_obstacleList();
        size_t size_oList = boundary->getSize_obstacleList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;
        real x, y, z;

        // set values into layers
        for (size_t l = 0; l < n_layers; ++l) {
            //inner cells
            for (size_t i = 0; i < size_iList; i++) {
                size_t idx = iList[i];
                coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
                coords_i = coords[0];
                coords_j = coords[1];
                coords_k = coords[2];

                if (dir == "x") {
                    x = xi(coords_i, X1, dx) - 0.5 * dx;
                    if (bord[l] <= x && x <= bord[l + 1]) out->data[idx] = val[l];

                } else if (dir == "y") {
                    y = yj(coords_j, Y1, dy) - 0.5 * dy;
                    if (bord[l] <= y && y <= bord[l + 1]) out->data[idx] = val[l];

                } else if (dir == "z") {
                    z = zk(coords_k, Z1, dz) - 0.5 * dz;
                    if (bord[l] <= z && z <= bord[l + 1]) out->data[idx] = val[l];
                } else {
#ifndef BENCHMARKING
                    auto m_logger = Utility::create_logger("Functions");
                    m_logger->error("No distance for layers specified!");
#endif
                }
                //TODO Error handling
            }

            //boundary
            for (size_t i = 0; i < size_bList; i++) {
                size_t idx = bList[i];
                coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
                coords_i = coords[0];
                coords_j = coords[1];
                coords_k = coords[2];

                if (dir == "x") {
                    x = xi(coords_i, X1, dx) - 0.5 * dx;
                    if (bord[l] <= x && x <= bord[l + 1]) out->data[idx] = val[l];
                    if (x < bord[0]) out->data[idx] = val[0];
                } else if (dir == "y") {
                    y = yj(coords_j, Y1, dy) - 0.5 * dy;
                    if (bord[l] <= y && y <= bord[l + 1]) out->data[idx] = val[l];
                    if (y < bord[0]) out->data[idx] = val[0];
                } else if (dir == "z") {
                    z = zk(coords_k, Z1, dz) - 0.5 * dz;
                    if (bord[l] <= z && z <= bord[l + 1]) out->data[idx] = val[l];
                    if (z < bord[0]) out->data[idx] = val[0];
                } else {
#ifndef BENCHMARKING
                    auto m_logger = Utility::create_logger("Functions");
                    m_logger->error("No distance for layers specified!");
#endif
                }
                //TODO Error handling
            }
            //obstacles
            for (size_t i = 0; i < size_oList; i++) {
                size_t idx = oList[i];
                coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
                coords_i = coords[0];
                coords_j = coords[1];
                coords_k = coords[2];

                if (dir == "x") {
                    x = xi(coords_i, X1, dx) - 0.5 * dx;
                    if (bord[l] <= x && x <= bord[l + 1]) out->data[idx] = val[l];
                    if (x < bord[0]) out->data[idx] = val[0];
                } else if (dir == "y") {
                    y = yj(coords_j, Y1, dy) - 0.5 * dy;
                    if (bord[l] <= y && y <= bord[l + 1]) out->data[idx] = val[l];
                    if (y < bord[0]) out->data[idx] = val[0];
                } else if (dir == "z") {
                    z = zk(coords_k, Z1, dz) - 0.5 * dz;
                    if (bord[l] <= z && z <= bord[l + 1]) out->data[idx] = val[l];
                    if (z < bord[0]) out->data[idx] = val[0];
                } else {
#ifndef BENCHMARKING
                    auto m_logger = Utility::create_logger("Functions");
                    m_logger->error("No distance for layers specified!");
#endif
                }
                //TODO Error handling
            }

        } //end layer loop
    }


// ============================= Diffusion Test - IC for u,v,w ===========================
// ***************************************************************************************
/// \brief  Initial set up for Diffusion Test
/// \param  out velocity
// ***************************************************************************************
    void Hat(Field *out) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(out->GetLevel());
        size_t Ny = domain->get_Ny(out->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(out->GetLevel());
        real dy = domain->get_dy(out->GetLevel());
        real dz = domain->get_dz(out->GetLevel());

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
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;

//inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];

            if ((start_x <= xi(coords_i, X1, dx) && xi(coords_i, X1, dx) <= end_x) &&
                (start_y <= yj(coords_j, Y1, dy) && yj(coords_j, Y1, dy) <= end_y) &&
                (start_z <= zk(coords_k, Z1, dz) && zk(coords_k, Z1, dz) <= end_z)) {
                out->data[idx] = val_in;
            } else {
                out->data[idx] = val_out;
            }
        }

//boundary
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];

            if ((start_x <= xi(coords_i, X1, dx) && xi(coords_i, X1, dx) <= end_x) &&
                (start_y <= yj(coords_j, Y1, dy) && yj(coords_j, Y1, dy) <= end_y) &&
                (start_z <= zk(coords_k, Z1, dz) && zk(coords_k, Z1, dz) <= end_z)) {
                out->data[idx] = val_in;
            } else out->data[idx] = val_out;
        }
    }

// ========================== NS Test - McDermott IC for u,v,w,p =========================
// ***************************************************************************************
/// \brief  Initial set up for NS Test - McDermott
/// \param  outx  x-velocity
/// \param  outy  y-velocity
/// \param  outz  z-velocity
/// \param  outp  pressure
/// \param  t   time
// ***************************************************************************************
    void McDermott(Field *outx, Field *outy, Field *outz, Field *outp, real t) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(outx->GetLevel());
        size_t Ny = domain->get_Ny(outx->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();

        real dx = domain->get_dx(outx->GetLevel());
        real dy = domain->get_dy(outx->GetLevel());

        auto params = Parameters::getInstance();
        real nu = params->get_real("physical_parameters/nu");

        real A = params->get_real("initial_conditions/A"); //2;

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            outx->data[idx] = 1. - A * cos(xi(coords_i, X1, dx) - t) * sin(yj(coords_j, Y1, dy) - t) * exp(-2 * nu * t);
            outy->data[idx] = 1. + A * sin(xi(coords_i, X1, dx) - t) * cos(yj(coords_j, Y1, dy) - t) * exp(-2 * nu * t);
            outz->data[idx] = 0.;
            outp->data[idx] =
                    -0.25 * A * A * (cos(2 * (xi(coords_i, X1, dx) - t)) + cos(2 * (yj(coords_j, Y1, dy) - t))) * exp(-4 * nu * t);
        }

        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            outx->data[idx] = 1. - A * cos(xi(coords_i, X1, dx) - t) * sin(yj(coords_j, Y1, dy) - t) * exp(-2 * nu * t);
            outy->data[idx] = 1. + A * sin(xi(coords_i, X1, dx) - t) * cos(yj(coords_j, Y1, dy) - t) * exp(-2 * nu * t);
            outz->data[idx] = 0.;
            outp->data[idx] =
                    -0.25 * A * A * (cos(2 * (xi(coords_i, X1, dx) - t)) + cos(2 * (yj(coords_j, Y1, dy) - t))) * exp(-4 * nu * t);
        }
    }

// ============================= Ramp up function for HRR source =========================
// ***************************************************************************************
/// \brief  Ramp up function (in time) for Gaussian temperature source in energy equation
/// \param  t time
// ***************************************************************************************
    real RampTanh(real t) {
        auto params = Parameters::getInstance();
        real tau = params->get_real("solver/temperature/source/tau");

        real rt = t / tau;
        real result = tanh(rt);

        return result;
    }

// === Random Function - Superposition of field values with random values ===
// ***************************************************************************************
/// \brief  Creates random absolute/relative noise on given field
/// \param  out     field
/// \param  range   range of random numbers
/// \param  is_absolute    Check if random number is relative (multiply) or absolute (additive)
/// \param  seed    custom seed if given, else seed <= 0
/// \param  step_size interval steps of random numbers
// ***************************************************************************************
    void Random(Field *out, real range, bool is_absolute, int seed, real step_size) {
        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        size_t *oList = boundary->get_obstacleList();
        size_t size_oList = boundary->getSize_obstacleList();

        std::mt19937 mt;
        double steps = range/step_size;
        if (seed > 0) {
          mt = std::mt19937(seed);
        } else {
          std::random_device rd;
          mt = std::mt19937(rd());
        }
        std::uniform_int_distribution<int> dist(-steps, steps);

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            // generate secret number between -range and range:
            if (is_absolute) {
                out->data[idx] += (dist(mt)*step_size);
            } else {
                out->data[idx] *= (1 + dist(mt)*step_size);
            }
        }
        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            // generate secret number between -range and range:
            if (is_absolute) {
                out->data[idx] += (dist(mt)*step_size);
            } else {
                out->data[idx] *= (1 + dist(mt)*step_size);
            }
        }
        // obstacles
        for (size_t i = 0; i < size_oList; i++) {
            size_t idx = oList[i];
            // generate secret number between -range and range:
            if (is_absolute) {
                out->data[idx] += (dist(mt)*step_size);
            } else {
                out->data[idx] *= (1 + dist(mt)*step_size);
            }
        }
    }

// ================================= Pressure Test - IC for p ============================
// ***************************************************************************************
/// \brief  Initial set up for Pressure Test (sin*sin*sin)
/// \param  out   pressure
// ***************************************************************************************
    void SinSinSin(Field *out) {
        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(out->GetLevel());
        size_t Ny = domain->get_Ny(out->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();
        real Z1 = domain->get_Z1();

        real dx = domain->get_dx(out->GetLevel());
        real dy = domain->get_dy(out->GetLevel());
        real dz = domain->get_dz(out->GetLevel());

        auto params = Parameters::getInstance();
        real l = params->get_real("initial_conditions/l"); //2;

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j, coords_k;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            out->data[idx] =
                    sin(l * M_PI * xi(coords_i, X1, dx)) * sin(l * M_PI * yj(coords_j, Y1, dy)) * sin(l * M_PI * zk(coords_k, Z1, dz));
        }
        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            coords_k = coords[2];
            out->data[idx] =
                    sin(l * M_PI * xi(coords_i, X1, dx)) * sin(l * M_PI * yj(coords_j, Y1, dy)) * sin(l * M_PI * zk(coords_k, Z1, dz));
        }
    }

// ======================= uniform distribution (eg. for force) ==========================
// ***************************************************************************************
/// \brief  Initial uniform set up
/// \param  out   force
/// \param  val   value of uniform distribution
// ***************************************************************************************
    void Uniform(Field *out, real val) {

        auto boundary = BoundaryController::getInstance();
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            out->data[idx] = val;
        }
        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            out->data[idx] = val;
        }
    }

// ============================= NS Test - Vortex IC for u,v,w,p =========================
// ***************************************************************************************
/// \brief  Initial set up for NS Test - Vertex
/// \param  outx    x-velocity
/// \param  outy    y-velocity
/// \param  outz    z-velocity
/// \param  outp    pressure
// ***************************************************************************************
    void Vortex(Field *outx, Field *outy, Field *outz, Field *outp) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(outx->GetLevel());
        size_t Ny = domain->get_Ny(outx->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();

        real dx = domain->get_dx(outx->GetLevel());
        real dy = domain->get_dy(outx->GetLevel());

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
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            outx->data[idx] = u_lin - GrR_c * yj(coords_j, Y1, dy) *
                                      exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
            outy->data[idx] = v_lin + GrR_c * xi(coords_i, X1, dx) *
                                      exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
            outz->data[idx] = 0.;
            outp->data[idx] =
                    pa - rhoGrR_c * exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
        }

        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            outx->data[idx] = u_lin - GrR_c * yj(coords_j, Y1, dy) *
                                      exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
            outy->data[idx] = v_lin + GrR_c * xi(coords_i, X1, dx) *
                                      exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
            outz->data[idx] = 0.;
            outp->data[idx] =
                    pa - rhoGrR_c * exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
        }
    }

    void VortexY(Field *outx, Field *outy, Field *outz, Field *outp) {

        auto domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(outx->GetLevel());
        size_t Ny = domain->get_Ny(outx->GetLevel());

        real X1 = domain->get_X1();
        real Y1 = domain->get_Y1();

        real dx = domain->get_dx(outx->GetLevel());
        real dy = domain->get_dy(outx->GetLevel());

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
        size_t *iList = boundary->get_innerList_level_joined();
        size_t size_iList = boundary->getSize_innerList();
        size_t *bList = boundary->get_boundaryList_level_joined();
        size_t size_bList = boundary->getSize_boundaryList();
        std::vector<size_t> coords;
        size_t coords_i, coords_j;

        // inner cells
        for (size_t i = 0; i < size_iList; i++) {
            size_t idx = iList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            outx->data[idx] = u_lin - GrR_c * yj(coords_j, Y1, dy) *
                                      exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
            outy->data[idx] = v_lin + GrR_c * xi(coords_i, X1, dx) *
                                      exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
            outz->data[idx] = 0.;
            outp->data[idx] =
                    pa - rhoGrR_c * exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
        }

        // boundary cells
        for (size_t i = 0; i < size_bList; i++) {
            size_t idx = bList[i];
            coords = Utility::coordinateFromLinearIndex(idx, Nx, Ny);
            coords_i = coords[0];
            coords_j = coords[1];
            outx->data[idx] = u_lin - GrR_c * yj(coords_j, Y1, dy) *
                                      exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
            outy->data[idx] = v_lin + GrR_c * xi(coords_i, X1, dx) *
                                      exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
            outz->data[idx] = 0.;
            outp->data[idx] =
                    pa - rhoGrR_c * exp(-rR_c * (xi(coords_i, X1, dx) * xi(coords_i, X1, dx) + yj(coords_j, Y1, dy) * yj(coords_j, Y1, dy)));
        }
    }

    void Zero(Field *field, size_t *arr_idx, size_t arr_idx_size) {
        auto data = field->data;
        for (size_t idx = 0; idx < arr_idx_size; idx++) {
            *(data + arr_idx[idx]) = 0;
        }
    }
}//end of namespace
