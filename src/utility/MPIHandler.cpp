/// \file       MPI.h
/// \brief      Custom MPI handler
/// \date       Sep 01, 2020
/// \author     Max Joseph Böhler
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "MPIHandler.h"
#include "Parameters.h"

#include "GlobalMacrosTypes.h"



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
                       m_MPIWORLD(MPIWORLD), m_MPICART(MPICART), m_dimensions(MPICART.topology().stl()), m_mpi_neighbour(6,0)  {
                           m_Xdim = m_dimensions.at(0).size; 
                           m_Ydim = m_dimensions.at(1).size; 
                           m_Zdim = m_dimensions.at(2).size; 
                           
                           check_mpi_neighbour();
                       }

void MPIHandler::exchange_data(real *data_field, size_t direction, size_t* d_patch) {
    std::vector< size_t > idx_inner;
    if      (direction == 0) idx_inner = m_inner_boundary_x1;
    else if (direction == 1) idx_inner = m_inner_boundary_x2;
    else if (direction == 2) idx_inner = m_inner_boundary_y1;
    else if (direction == 3) idx_inner = m_inner_boundary_y2;
    else if (direction == 4) idx_inner = m_inner_boundary_z1;
    else if (direction == 5) idx_inner = m_inner_boundary_z2;

   for (size_t i = 0; i < idx_inner.size(); i++)
    {
        data_field[d_patch[i]] = data_field[idx_inner.at(i)];
    }
}

void MPIHandler::get_inner_index(){
    Domain *domain = Domain::getInstance();
    size_t i1 = domain->get_index_x1() - 1;
    size_t j1 = domain->get_index_y1() - 1;
    size_t k1 = domain->get_index_z1() - 1;

    size_t i2 = domain->get_index_x2() + 1;
    size_t j2 = domain->get_index_y2() + 1;
    size_t k2 = domain->get_index_z2() + 1;
   
    size_t Nx = domain->get_Nx();
    size_t Ny = domain->get_Ny();


    for (size_t k = k1; k <= k2; ++k) {
        for (size_t j = j1; j <= j2; ++j) {
            m_inner_boundary_x1.push_back(IX(i1 + 1, j, k, Nx, Ny));
            m_inner_boundary_x2.push_back(IX(i2 - 1, j, k, Nx, Ny));
        }
    }

    for (size_t k = k1; k <= k2; ++k) {
        for (size_t i = i1; i <= i2; ++i) {
            m_inner_boundary_y1.push_back(IX(i, j1 + 1, k, Nx, Ny));
            m_inner_boundary_y2.push_back(IX(i, j2 - 1, k, Nx, Ny));
        }
    }

    for (size_t j = j1; j <= j2; ++j) {
        for (size_t i = i1; i <= i2; ++i) {
            m_inner_boundary_z1.push_back(IX(i, j, k1 + 1, Nx, Ny));
            m_inner_boundary_z2.push_back(IX(i, j, k2 - 1, Nx, Ny));
        }
    }
}

void MPIHandler::check_mpi_neighbour() {
    std::vector<int> rank_coordinates;
    rank_coordinates = get_coords();

    m_mpi_neighbour.at(0) = (rank_coordinates.at(0) == 0) ? 0 : 1;
    m_mpi_neighbour.at(1) = (rank_coordinates.at(0) == m_Xdim - 1) ? 0 : 1;

    m_mpi_neighbour.at(2) = (rank_coordinates.at(1) == 0) ? 0 : 1;
    m_mpi_neighbour.at(3) = (rank_coordinates.at(1) == m_Ydim - 1) ? 0 : 1;

    m_mpi_neighbour.at(4) = (rank_coordinates.at(2) == 0) ? 0 : 1;
    m_mpi_neighbour.at(5) = (rank_coordinates.at(2) == m_Zdim - 1) ? 0 : 1;
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
    real domain_start, domain_end;
    if (direction == 0) {
        domain_start = domain->get_X1();
        domain_end = domain->get_X2();
    } else if (direction == 1) {
        domain_start = domain->get_Y1();
        domain_end = domain->get_Y2();
    } else {
        domain_start = domain->get_Z1();
        domain_end = domain->get_Z2();
    }

    if (x1 < domain_start && x2 > domain_end ) {
      x1 = domain_start;
      x2 = domain_end;
      return true;
  } else if (x1 < domain_start && x2 > domain_start && x2 < domain_end) {
      x1 = domain_start;
      return true;
  } else if (x1 > domain_start && x1 < domain_end && x2 > domain_end) {
      x2 = domain_end;
      return true;
  } else if (x1 >= domain_start && x2 <= domain_end) {
      return true;
  } else if (x1 == domain_start && x2 > domain_end) {
      x1 = domain_start;
      x2 = domain_end;
      return true;
  } else if (x1 < domain_start && x2 == domain_end) {
      x1 = domain_start;
      x2 = domain_end;
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
