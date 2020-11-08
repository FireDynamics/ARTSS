/// \file       main.cpp
/// \brief      Hybrid (CPU and GPU) Version for solving Navier-Stokes equations in 3D
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include "TimeIntegration.h"
#include "utility/tinyxml2.h"
#include "utility/Utility.h"
#include "utility/Parameters.h"
#include "solver/SolverController.h"

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

#ifdef USEMPI
    int MPIX, MPIY, MPIZ;

    char* temp_x = Utility::get_flag(argv, argv + argc, "-x");
    char* temp_y = Utility::get_flag(argv, argv + argc, "-y");
    char* temp_z = Utility::get_flag(argv, argv + argc, "-z");
    
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

    SolverController *sc = new SolverController();

#ifdef _OPENACC
    // Initialise GPU
    acc_device_t dev_type = acc_get_device_type();
    acc_init( dev_type );
#endif

    // Integrate over time and solve numerically
    // Time integration
    TimeIntegration ti(sc);
    ti.run();

    delete sc;

#ifdef USEMPI
    MPI_Finalize();
#endif
    return 0;
}
