/// \file       MPI.h
/// \brief      Custom MPI handler
/// \date       Sep 01, 2020
/// \author     Max Joseph Böhler
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "MPIHandler.h"
#include "Parameters.h"



MPIHandler *MPIHandler::single = nullptr;

// Singleton
MPIHandler *MPIHandler::getInstance() {
    return single;
}

// Singleton
MPIHandler *MPIHandler::getInstance(boost::mpi::communicator& MPIWORLD,
                      boost::mpi::cartesian_communicator& MPICART) {
    if (single == nullptr) {
        single = new MPIHandler(MPIWORLD, MPICART);
    }
    return single;
}

MPIHandler::MPIHandler(boost::mpi::communicator& MPIWORLD,
                       boost::mpi::cartesian_communicator& MPICART) :
                       m_MPIWORLD(MPIWORLD), m_MPICART(MPICART), m_dimensions(MPICART.topology().stl()) {
                         #ifndef BENCHMARKING
                             m_logger = Utility::create_logger(typeid(this).name());
                         #endif
                       }

void MPIHandler::convert_domain(real& x1, real& x2, int direction) {
    int  dimension      { m_dimensions[direction].size };
    real local_size     { (x2 - x1) /  dimension};
    int  rank_dimension { m_MPICART.coordinates(m_MPICART.rank())[direction] };

    x1 =+ rank_dimension * local_size;
    x2 =+ (rank_dimension + 1) * local_size;
}

int MPIHandler::convert_grid(std::string param, int direction) {
    auto params = Parameters::getInstance();
    int local_size = params->get_real(param) / m_dimensions[direction].size;

    return local_size;
}
