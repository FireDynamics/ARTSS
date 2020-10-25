/// \file       main.cpp
/// \brief      Hybrid (CPU and GPU) Version for solving Navier-Stokes equations
///             in 3D
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <unistd.h>
#include <getopt.h>

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

#ifdef USEMPI
    #include <mpi.h>
    #include "utility/MPIHandler.h"
#endif

#ifdef _OPENACC
    #include <openacc.h>
#endif


int main(int argc, char **argv) {

#ifdef USEMPI
    int MPIWORLD;
    int MPIRANK;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &MPIWORLD);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIRANK);
#endif


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

#ifdef USEMPI
    int MPIX, MPIY, MPIZ;

    char * temp_x = Utility::get_flag(argv, argv + argc, "-x");
    char * temp_y = Utility::get_flag(argv, argv + argc, "-y");
    char * temp_z = Utility::get_flag(argv, argv + argc, "-z");

    if ( !temp_x && !temp_y && !temp_z ) {
        MPIX = params->get_int("domain_parameters/MESHX");
        MPIY = params->get_int("domain_parameters/MESHY");
        MPIZ = params->get_int("domain_parameters/MESHZ");
    } else {  
        MPIX = (temp_x) ? atoi(temp_x) : 1;
        MPIY = (temp_y) ? atoi(temp_y) : 1;
        MPIZ = (temp_z) ? atoi(temp_z) : 1;
    }

    if (MPIWORLD != MPIX*MPIY*MPIZ) {
       if (MPIRANK == 0) std::cerr << "Number of meshes does not match the number of defined MPI processes!" << std::endl;
       MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Comm MPICART;

    int dimensions[3] = {MPIX, MPIY, MPIZ};
    int periodic[3];

    periodic[0] = 1;
    periodic[1] = 1;
    periodic[2] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 3, dimensions, periodic, 1, &MPICART);

    auto mpi_handler = MPIHandler::getInstance(MPICART, dimensions);
#endif


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

    // Integrate over time and solve numerically
    // Time integration
    TimeIntegration ti(solver);
    ti.run();

    // Clean up
    delete solver;

#ifdef USEMPI
    MPI_Finalize();
#endif

    return 0;
}
