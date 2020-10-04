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
#include <boost/serialization/vector.hpp>

#include "GlobalMacrosTypes.h"
#include "Parameters.h"
#include "Utility.h"
#include "../Domain.h"
#include "../boundary/BoundaryData.h"


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
    void calc_inner_index(size_t level);
    void exchange_data(real *data_field, Patch p, size_t* d_patch, const size_t patch_starts);

    void set_barrier() { return m_MPICART.barrier(); }

private:
    static MPIHandler* single;

    MPIHandler(boost::mpi::communicator& MPIWORLD, boost::mpi::cartesian_communicator& MPICART);
    boost::mpi::communicator m_MPIWORLD;
    boost::mpi::cartesian_communicator m_MPICART;
    std::vector< boost::mpi::cartesian_dimension > m_dimensions;
    std::vector< size_t > m_inner_left;
    std::vector< size_t > m_inner_right;
    std::vector< size_t > m_inner_bottom;
    std::vector< size_t > m_inner_top;
    std::vector< size_t > m_inner_front;
    std::vector< size_t > m_inner_back;

    int m_Xdim;
    int m_Ydim;
    int m_Zdim;

    void check_mpi_neighbour();

    std::vector<int> m_mpi_neighbour;
    std::vector< std::pair <int, int > > m_mpi_neighbour_rank_offset;

};

#endif /* ARTSS_UTILITY_MPI_H */
