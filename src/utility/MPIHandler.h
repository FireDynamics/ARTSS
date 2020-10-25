/// \file       MPIHandler.h
/// \brief      Custom MPI handler
/// \date       October, 2020
/// \author     Max Joseph Böhler
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_MPIHANDLER_H
#define ARTSS_UTILITY_MPIHANDLER_H

#include <iostream>
#include <string>

#include <mpi.h>

#include "GlobalMacrosTypes.h"
#include "Parameters.h"
#include "Utility.h"
#include "../Domain.h"
#include "../boundary/BoundaryData.h"


class MPIHandler {
 public:
    static MPIHandler* getInstance();
    static MPIHandler* getInstance(MPI_Comm MPICART, int* dimensions);

    // Getter
    int get_rank() { return m_CARTRANK; }
    std::vector<int> get_coords() { return {m_XRANK, m_YRANK, m_ZRANK}; }
    std::vector<int> get_mpi_neighbour() { return m_mpi_neighbour; }

    void convert_domain(real& x1, real& x2, int direction);
    void convert_grid(size_t& n, int direction);
    bool convert_obstacle(real& x1, real& x2, int direction);
    bool has_obstacle(real& ox1, real& ox2, real& oy1, real& oy2, real& oz1, real& oz2);
    void calc_inner_index();
    void exchange_data(real *data_field, size_t** index_fields, const size_t *patch_starts, size_t level);
    void exchange_data(real *data_field, size_t** index_fields, const size_t *patch_starts, size_t level, std::vector<bool> periodic);

    int set_barrier() { return MPI_Barrier(m_MPICART); }

    double get_max_val(double val);

private:
    static MPIHandler* single;

    MPIHandler(MPI_Comm MPICART, int* dimensions);
    MPI_Comm m_MPICART;
    MPI_Comm m_XCOM;
    MPI_Comm m_YCOM;
    MPI_Comm m_ZCOM;


    std::vector< std::vector<size_t> > m_inner_left;
    std::vector< std::vector<size_t> > m_inner_right;
    std::vector< std::vector<size_t> > m_inner_bottom;
    std::vector< std::vector<size_t> > m_inner_top;
    std::vector< std::vector<size_t> > m_inner_front;
    std::vector< std::vector<size_t> > m_inner_back;
    std::vector<int> m_mpi_neighbour;
    std::vector<int> m_mpi_neighbour_periodic;
    std::vector<int> m_dimensions;

    int m_sendrecv_ctr;
    int m_Xdim;
    int m_Ydim;
    int m_Zdim;
    int m_CARTRANK;
    int m_XRANK;
    int m_YRANK;
    int m_ZRANK;

    void check_mpi_neighbour();
};

#endif /* ARTSS_UTILITY_MPIHANDLER_H */
