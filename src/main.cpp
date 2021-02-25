/// \file       main.cpp
/// \brief      Hybrid (CPU and GPU) Version for solving Navier-Stokes equations in 3D
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include "TimeIntegration.h"
#include "utility/tinyxml2.h"
#include "utility/Parameters.h"
#include "solver/SolverController.h"

#ifdef _OPENACC
    #include <openacc.h>
#endif

#include <iostream>
#include "source/GaussFunction.h"
#include "boundary/BoundaryController.h"

void test() {
    auto params = Parameters::getInstance();
    auto domain = Domain::getInstance();
    auto out = new Field(FieldType::UNKNOWN_FIELD, 0.0);

    real HRR = params->get_real("solver/temperature/source/HRR");    // heat release rate in [kW]
    real cp = params->get_real("solver/temperature/source/cp");        // specific heat capacity in [kJ/ kg K]
    real x0 = params->get_real("solver/temperature/source/x0");
    real y0 = params->get_real("solver/temperature/source/y0");
    real z0 = params->get_real("solver/temperature/source/z0");
    real sigma_x = params->get_real("solver/temperature/source/sigma_x");
    real sigma_y = params->get_real("solver/temperature/source/sigma_y");
    real sigma_z = params->get_real("solver/temperature/source/sigma_z");
    real tau = params->get_real("solver/temperature/source/tau");
    auto m_source_function_temperature = new GaussFunction(HRR, cp, x0, y0, z0, sigma_x, sigma_y, sigma_z, tau);

    m_source_function_temperature->update_source(out, 1.0);
    
    std::cout << domain->get_size() << std::endl;
    std::ofstream file("test.csv");

    auto boundary = BoundaryController::getInstance();

    const auto multigrid = boundary->getMultigrid();
    const auto obst_size = multigrid->getSize_NumberOfObstacles();
    const auto obst_list = multigrid->getObstacles();

    const auto Nx = domain->get_Nx();
    const auto Ny = domain->get_Ny();
    file << "index" << "," << "i" << ","
        << "j" << ","
        << "k" << ","
        << "data" << ","
        << "blocked" << std::endl;
    std::cout << obst_size << std::endl;
    for (int i=0; i < domain->get_size(); ++i) {
        auto coords_k = getCoordinateK(i, Nx, Ny);
        auto coords_j = getCoordinateJ(i, Nx, Ny, coords_k);
        auto coords_i = getCoordinateI(i, Nx, Ny, coords_j, coords_k);
        auto blocked = false;
        for (size_t obst_id=0; obst_id < obst_size; ++obst_id) {
            blocked = blocked || obst_list[0][obst_id]->isObstacleCell(coords_i, coords_j, coords_k);
        }
        file << i << "," << coords_i << ","
            << coords_j << ","
            << coords_k << ","
            << out->data[i] << ","
            << static_cast<int>(blocked) << std::endl;
    }
}

int main(int argc, char **argv) {
    // Initialisation
    // Parameters
    std::string XML_filename;
    auto params = Parameters::getInstance();
    if (argc > 1) {
        XML_filename.assign(argv[1]);
        params->parse(XML_filename);
    } else {
        std::cerr << "XML file missing" << std::endl;
        std::exit(1);
    }

    // SolverController *sc = new SolverController();

#ifdef _OPENACC
    // Initialise GPU
    acc_device_t dev_type = acc_get_device_type();
    acc_init( dev_type );
#endif

    // Integrate over time and solve numerically
    // Time integration
    // TimeIntegration ti(sc);
    // ti.run();
    test();

    // Clean up
    // delete sc;
    return 0;
}
