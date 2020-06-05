/// \file 		AdvectionDiffusionSolver.h
/// \brief 		Defines the steps to solve the advection and diffusion equation
/// \date 		May 20, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_
#define ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_

#include "../interfaces/SolverI.h"
#include "../interfaces/AdvectionI.h"
#include "../interfaces/DiffusionI.h"
#include "../utility/GlobalMacrosTypes.h"
#ifndef PROFILING
#include <spdlog/logger.h>
#endif

class AdvectionDiffusionSolver : public SolverI {
public:
	AdvectionDiffusionSolver();
	~AdvectionDiffusionSolver() override;

	void DoStep(real t, bool sync) override;

private:
	AdvectionI* adv;
	DiffusionI* dif;

	real m_nu;

    static void control();
#ifndef PROFILING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_ */
