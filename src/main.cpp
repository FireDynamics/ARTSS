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

int Field::counter = 0;
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

#ifdef _OPENACC
    // Initialise GPU
    acc_device_t dev_type = acc_get_device_type();
    acc_init( dev_type );
#endif

    SolverController *sc = new SolverController();

    // Integrate over time and solve numerically
    // Time integration
    TimeIntegration ti(sc);
    ti.run();

    // Clean up
    delete sc;
    return 0;
}
