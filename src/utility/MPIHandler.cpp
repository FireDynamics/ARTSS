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

bool MPIHandler::convert_obstacle(real& x1, real& x2, int direction) {
    auto domain   = Domain::getInstance();
    real domainX1, domainX2;
    if (direction == 0) {
        domainX1 = domain->get_X1();
        domainX2 = domain->get_X2();
    } else if (direction == 1) {
        domainX1 = domain->get_Y1();
        domainX2 = domain->get_Y2();
    } else {
        domainX1 = domain->get_Z1();
        domainX2 = domain->get_Z2();
    }

    if (x1 < domainX1 && x2 > domainX2 ) {
      x1 = domainX1;
      x2 = domainX2;
      return true;
  } else if (x1 < domainX1 && x2 > domainX1 && x2 < domainX2) {
      x1 = domainX1;
      return true;
  } else if (x1 > domainX1 && x1 < domainX2 && x2 > domainX2) {
      x2 = domainX2;
      return true;
  } else if (x1 >= domainX1 && x2 <= domainX2) {
      return true;
  } else {
      return false;
  }
}

bool MPIHandler::has_obstacle(real& ox1, real& ox2, real& oy1, real& oy2, real& oz1, real& oz2) {
    bool checkX, checkY, checkZ;
    checkX = convert_obstacle(ox1, ox2, 0);
    checkY = convert_obstacle(oy1, oy2, 1);
    checkZ = convert_obstacle(oz1, oz2, 2);

    if (checkX == false || checkY == false || checkZ == false) return false;

    return true;
}

int MPIHandler::convert_grid(std::string param, int direction) {
    auto params = Parameters::getInstance();
    int local_size = params->get_real(param) / m_dimensions[direction].size;

    return local_size;
}
