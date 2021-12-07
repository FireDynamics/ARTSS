/// \file       main.cpp
/// \brief      Hybrid (CPU and GPU) Version for solving Navier-Stokes equations in 3D
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include "TimeIntegration.h"
#include "DomainData.h"
#include "boundary/BoundaryController.h"
#include "solver/SolverController.h"
#include "utility/settings/Settings.h"
#include "visualisation/VTKWriter.h"

#ifdef _OPENACC
    #include <openacc.h>
#endif

int VTKWriter::vtk_counter = 0;
int Field::counter = 0;
int main(int argc, char **argv) {
    // Initialisation
    // Parameters
    std::string XML_filename;
    if (argc <= 1) {
        std::cerr << "XML file missing" << std::endl;
        std::exit(1);
    }

#ifdef _OPENACC
    // Initialise GPU
    acc_device_t dev_type = acc_get_device_type();
    acc_init(dev_type);
#endif
    Settings::Settings settings(argv[1]);
    DomainData::getInstance(settings);
    BoundaryController::getInstance(settings);

    SolverController *sc = new SolverController(settings);

    // Integrate over time and solve numerically
    // Time integration
    TimeIntegration ti(settings, sc);
    ti.run();

    // Clean up
    delete sc;
    return 0;
}
