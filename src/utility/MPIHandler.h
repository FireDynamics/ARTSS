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


class MPIHandler {
 public:
    static MPIHandler* getInstance();
    static MPIHandler* getInstance(boost::mpi::communicator& MPIWORLD,
                            boost::mpi::cartesian_communicator& MPICART);

    // Getter
    int get_rank() { return m_MPICART.rank(); };
    std::vector<int> get_coords() { return m_MPICART.coordinates(m_MPICART.rank()); }

    void convert_domain(real& x1, real& x2, int direction);
    int  convert_grid(std::string param, int direction);

private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    static MPIHandler* single;

    MPIHandler(boost::mpi::communicator& MPIWORLD, boost::mpi::cartesian_communicator& MPICART);
    boost::mpi::communicator m_MPIWORLD;
    boost::mpi::cartesian_communicator m_MPICART;
    std::vector< boost::mpi::cartesian_dimension > m_dimensions;
};

#endif /* ARTSS_UTILITY_MPI_H */
