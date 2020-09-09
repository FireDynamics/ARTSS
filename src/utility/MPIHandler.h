/// \file       MPI.h
/// \brief      Custom MPI handler
/// \date       Sep 01, 2020
/// \author     Max Joseph Böhler
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_MPI_H
#define ARTSS_UTILITY_MPI_H

#include <iostream>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/cartesian_communicator.hpp>


class MPIHandler {

private:
    static MPIHandler* single;

    MPIHandler(boost::mpi::communicator& MPIWORLD, boost::mpi::cartesian_communicator& MPICART);
    boost::mpi::communicator MPIWORLD;
    boost::mpi::cartesian_communicator MPICART;

public:
    static MPIHandler* getInstance();
    static MPIHandler* getInstance(boost::mpi::communicator& MPIWORLD,
                            boost::mpi::cartesian_communicator& MPICART);

    // Getter
    int get_rank() { return MPICART.rank(); };
    std::vector<int> get_coords() { return MPICART.coordinates(MPICART.rank()); }

};

#endif /* ARTSS_UTILITY_MPI_H */
