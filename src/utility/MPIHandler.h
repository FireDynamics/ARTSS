/// \file       MPI.h
/// \brief      Custom MPI handler
/// \date       Sep 01, 2020
/// \author     Max Joseph Böhler
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_MPI_H
#define ARTSS_UTILITY_MPI_H

#include <iostream>
#include <string>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/cartesian_communicator.hpp>

#include "GlobalMacrosTypes.h"
#include "Parameters.h"
#include "Utility.h"
#include "../Domain.h"


class MPIHandler {
 public:
    static MPIHandler* getInstance();
    static MPIHandler* getInstance(boost::mpi::communicator& MPIWORLD,
                            boost::mpi::cartesian_communicator& MPICART);

    // Getter
    int get_rank() { return m_MPICART.rank(); };
    std::vector<int> get_coords() { return m_MPICART.coordinates(m_MPICART.rank()); }
    std::vector<int> get_mpi_neighbour() { return m_mpi_neighbour; }

    void convert_domain(real& x1, real& x2, int direction);
    int  convert_grid(std::string param, int direction);
    bool convert_obstacle(real& x1, real& x2, int direction);
    bool has_obstacle(real& ox1, real& ox2, real& oy1, real& oy2, real& oz1, real& oz2);
    void get_inner_index();
    void exchange_data(real *data_field, size_t direction, size_t* d_patch);

private:
    static MPIHandler* single;

    MPIHandler(boost::mpi::communicator& MPIWORLD, boost::mpi::cartesian_communicator& MPICART);
    boost::mpi::communicator m_MPIWORLD;
    boost::mpi::cartesian_communicator m_MPICART;
    std::vector< boost::mpi::cartesian_dimension > m_dimensions;
    std::vector< size_t > m_inner_boundary_x1;
    std::vector< size_t > m_inner_boundary_x2;
    std::vector< size_t > m_inner_boundary_y1;
    std::vector< size_t > m_inner_boundary_y2;
    std::vector< size_t > m_inner_boundary_z1;
    std::vector< size_t > m_inner_boundary_z2;

    int m_Xdim;
    int m_Ydim;
    int m_Zdim;

    void check_mpi_neighbour();

    std::vector<int> m_mpi_neighbour;

};

#endif /* ARTSS_UTILITY_MPI_H */
