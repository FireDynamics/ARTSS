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
#ifdef ASSIMILATION
    #include <mpi.h>
    #define PORT 7777

    void server();
#endif

int main(int argc, char **argv) {
    try {
#ifdef ASSIMILATION
        MPI_Init(nullptr, nullptr);
        int num_procs;
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

        int comm_size, comm_rank;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

        if (num_procs < 2 && comm_rank == 0) {
            // TODO Error Handling
            std::cout << "Data assimilation has to be run with at least 2 processes, "
                         "currently available: " << num_procs << " processes" << std::endl;
            std::exit(1);
        }

        if (comm_rank == 0) {
#endif
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
#ifdef ASSIMILATION
        MPI_Request request;
        MPI_Isend(settings.logging_parameters.level.c_str(), static_cast<int>(settings.logging_parameters.level.size()), MPI_CHAR, 1, 0, MPI_COMM_WORLD, &request);
#endif
        size_t multigrid_level = 0;

        auto solver = settings.solver_parameters.description;
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
#ifdef ASSIMILATION
        }
        if (comm_rank == 1) {
            server();
        }
    MPI_Finalize();
#endif

    } catch (const std::exception &ex) {
        std::cerr << ex.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown exception encountered" << std::endl;
    }
    return 0;
}
