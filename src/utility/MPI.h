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


class MPI {

private:
    static MPI* single;

    MPI(boost::mpi::communicator& MPIWORLD, boost::mpi::cartesian_communicator& MPICART);
    boost::mpi::communicator MPIWORLD;
    boost::mpi::cartesian_communicator MPICART;

public:
    static MPI* getInstance();
    static MPI* getInstance(boost::mpi::communicator& MPIWORLD,
                            boost::mpi::cartesian_communicator& MPICART);

    // Getter
    int getRank() { return MPICART.rank(); };

};

#endif /* ARTSS_UTILITY_MPI_H */
