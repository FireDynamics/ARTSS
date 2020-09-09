/// \file       DiffusionTurbSolver.h
/// \brief      Defines the steps to solve the turbulent diffusion equation
/// \date       Aug 18, 2016
/// \author     Suryanarayana Maddu
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_TURBDIFFUSIONSOLVER_H_
#define ARTSS_SOLVER_TURBDIFFUSIONSOLVER_H_

#include "../Domain.h"
#include "../interfaces/ISolver.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/ITurbulence.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Parameters.h"
#include "SolverSelection.h"
#include "../utility/Utility.h"


class DiffusionTurbSolver: public ISolver {
public:
    DiffusionTurbSolver(boost::mpi::cartesian_communicator& MPICART);
    ~DiffusionTurbSolver() override;

    void do_step(real t, bool sync) override;

private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    IDiffusion*   dif;
    ITurbulence*  mu_tub;

    real m_nu;

    static void control();
};

#endif /* ARTSS_SOLVER_TURBDIFFUSIONSOLVER_H_ */
