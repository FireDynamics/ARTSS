/// \file       main.cpp
/// \brief      Hybrid (CPU and GPU) Version for solving Navier-Stokes equations
///             in 3D
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <unistd.h>
#include <getopt.h>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <iostream>

#include "interfaces/ISolver.h"
#include "solver/DiffusionTurbSolver.h"
#include "solver/DiffusionSolver.h"
#include "solver/AdvectionSolver.h"
#include "solver/AdvectionDiffusionSolver.h"
#include "solver/NSSolver.h"
#include "solver/NSTurbSolver.h"
#include "solver/NSTempSolver.h"
#include "solver/NSTempTurbSolver.h"
#include "solver/NSTempConSolver.h"
#include "solver/NSTempTurbConSolver.h"
#include "solver/PressureSolver.h"
#include "TimeIntegration.h"
#include "utility/tinyxml2.h"
#include "utility/Parameters.h"
#include "analysis/Analysis.h"
#include "utility/Utility.h"
#include "utility/MPI.h"

#ifdef _OPENACC
    #include <openacc.h>
#endif


int main(int argc, char **argv) {

    boost::mpi::environment MPIENV;
    boost::mpi::communicator MPIWORLD;


    // Initialization
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

    const int MPIX{params->get_int("domain_parameters/MESHX")};
    const int MPIY{params->get_int("domain_parameters/MESHY")};
    const int MPIZ{params->get_int("domain_parameters/MESHZ")};

    if (MPIWORLD.size() != MPIX*MPIY*MPIZ) {
       if (MPIWORLD.rank() == 0) std::cerr << "Number of meshes does not match the number of defined MPI processes!" << std::endl;
       MPIENV.abort(1);
    }

    boost::mpi::cartesian_dimension MPIDIM[] = {{MPIX, false}, {MPIY, false}, {MPIZ, false}};
    boost::mpi::cartesian_communicator MPICART(MPIWORLD, boost::mpi::cartesian_topology(MPIDIM), true);

    auto mpi = MPI::getInstance(MPIWORLD, MPICART);



    // Solver
    ISolver *solver;
    std::string string_solver = params->get("solver/description");
    if (string_solver == SolverTypes::DiffusionSolver) solver = new DiffusionSolver();
    else if (string_solver == SolverTypes::AdvectionSolver) solver = new AdvectionSolver();
    else if (string_solver == SolverTypes::AdvectionDiffusionSolver) solver = new AdvectionDiffusionSolver();
    else if (string_solver == SolverTypes::PressureSolver) solver = new PressureSolver();
    else if (string_solver == SolverTypes::DiffusionTurbSolver) solver = new DiffusionTurbSolver();
    else if (string_solver == SolverTypes::NSSolver) solver = new NSSolver();
    else if (string_solver == SolverTypes::NSTurbSolver) solver = new NSTurbSolver();
    else if (string_solver == SolverTypes::NSTempSolver) solver = new NSTempSolver();
    else if (string_solver == SolverTypes::NSTempTurbSolver) solver = new NSTempTurbSolver();
    else if (string_solver == SolverTypes::NSTempConSolver) solver = new NSTempConSolver();
    else if (string_solver == SolverTypes::NSTempTurbConSolver) solver = new NSTempTurbConSolver();

    else {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger("main");
        m_logger->error("Solver not yet implemented! Simulation stopped!");
#endif
        std::exit(1);
    }

#ifdef _OPENACC
    // Initialize GPU
    acc_device_t dev_type = acc_get_device_type();
    acc_init( dev_type );
#endif
// //
// //     // Integrate over time and solve numerically
// //     // Time integration
// //     TimeIntegration ti(solver);
// //     ti.run();
// //
//     // Clean up
//     delete solver;

    return 0;
}
