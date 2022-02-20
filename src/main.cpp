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

#include <iostream>
#include "source/GaussFunction.h"
#include "boundary/BoundaryController.h"


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

        Settings::Settings settings = Settings::parse_settings_from_file(argv[1]);
        size_t multigrid_level = 0;

        const auto &solver = settings.solver_parameters.description;
        if (solver.find("NS") != std::string::npos || solver == SolverTypes::PressureSolver) {
            multigrid_level = settings.solver_parameters.pressure.solver.n_level;
        }
        DomainData::init(settings.physical_parameters, settings.domain_parameters, multigrid_level);
        DomainController::init(settings);

        SolverController *sc = new SolverController(settings);

        // Integrate over time and solve numerically
        // Time integration
        TimeIntegration ti(settings, sc);
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
