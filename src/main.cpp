/// \file       main.cpp
/// \brief      Hybrid (CPU and GPU) Version for solving Navier-Stokes equations in 3D
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include "TimeIntegration.h"
#include "domain/DomainData.h"
#include "domain/DomainController.h"
#include "solver/SolverController.h"
#include "utility/settings/Settings.h"

#ifdef _OPENACC
    #include <openacc.h>
#endif

int main(int argc, char **argv) {
    try {
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

        Settings::Settings_new settings_new = Settings::parse_settings_from_file(argv[1]);
        Settings::Settings settings(argv[1]);
        DomainData::init(settings_new);
        DomainController::init(settings);

        SolverController *sc = new SolverController(settings_new);

        // Integrate over time and solve numerically
        // Time integration
        TimeIntegration ti(settings_new, sc);
        ti.run();

        // Clean up
        delete sc;
    } catch (const std::exception &ex) {
        std::cerr << ex.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception encountered" << std::endl;
    }
    return 0;
}
