/// \file       main.cpp
/// \brief      Hybrid (CPU and GPU) Version for solving Navier-Stokes equations in 3D
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include "TimeIntegration.h"
#include "Domain.h"
#include "boundary/BoundaryController.h"
#include "utility/tinyxml2.h"
#include "utility/Settings.h"
#include "solver/SolverController.h"

#ifdef _OPENACC
    #include <openacc.h>
#endif

int main(int argc, char **argv) {
    // Initialisation
    // Parameters
    std::string XML_filename;
    if (argc <= 1) {
        std::cerr << "XML file missing" << std::endl;
        std::exit(1);
    }

    Settings settings(argv[1]);
    // Domain::getInstance(settings);
    // BoundaryController::getInstance(settings);

    settings.print_config();
    return 0;

    SolverController *sc = new SolverController(settings);

#ifdef _OPENACC
    // Initialise GPU
    acc_device_t dev_type = acc_get_device_type();
    acc_init(dev_type);
#endif

    // Integrate over time and solve numerically
    // Time integration
    TimeIntegration ti(settings, sc);
    ti.run();

    // Clean up
    delete sc;
    return 0;
}
