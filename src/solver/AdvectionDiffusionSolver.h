/// \file 		AdvectionDiffusionSolver.h
/// \brief 		Defines the steps to solve the advection and diffusion equation
/// \date 		May 20, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_
#define ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../utility/GlobalMacrosTypes.h"

class AdvectionDiffusionSolver : public ISolver {
public:
	AdvectionDiffusionSolver();
	~AdvectionDiffusionSolver() override;

	void DoStep(real t, bool sync) override;

private:
	IAdvection* adv;
	IDiffusion* dif;

	real m_nu;

    static void control();
};

#endif /* ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_ */
