/// \file       MPI.h
/// \brief      Custom MPI handler
/// \date       Sep 01, 2020
/// \author     Max Joseph Böhler
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "MPI.h"
#include "Parameters.h"


MPI *MPI::single = nullptr;

// Singleton
MPI *MPI::getInstance() {
    return single;
}

// Singleton
MPI *MPI::getInstance(boost::mpi::communicator& MPIWORLD,
                      boost::mpi::cartesian_communicator& MPICART) {
    if (single == nullptr) {
        single = new MPI(MPIWORLD, MPICART);
    }
    return single;
}

MPI::MPI(boost::mpi::communicator& MPIWORLD,
         boost::mpi::cartesian_communicator& MPICART)
         : MPIWORLD(MPIWORLD), MPICART(MPICART) {}
