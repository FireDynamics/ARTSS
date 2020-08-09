/// \file     AdvectionDiffusionSolver.h
/// \brief    Defines the steps to solve the advection and diffusion equation
/// \date     May 20, 2016
/// \author   Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_
#define ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../utility/GlobalMacrosTypes.h"
#ifndef PROFILING
#include <spdlog/logger.h>
#endif

class AdvectionDiffusionSolver : public ISolver {
public:
  AdvectionDiffusionSolver();
  ~AdvectionDiffusionSolver() override;

  void do_step(real t, bool sync) override;

private:
  IAdvection* adv;
  IDiffusion* dif;

  real m_nu;

    static void control();
#ifndef PROFILING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_ */
