/// \file 		DiffusionSolver.h
/// \brief 		Defines the steps to solve the diffusion equation
/// \date 		May 20, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_DIFFUSIONSOLVER_H_
#define ARTSS_SOLVER_DIFFUSIONSOLVER_H_

#include "../interfaces/SolverI.h"
#include "../interfaces/DiffusionI.h"

class DiffusionSolver: public SolverI {
public:
	DiffusionSolver();
	~DiffusionSolver() override;

	void DoStep(real t, bool sync) override;

private:
	DiffusionI* dif;
	real m_nu;

    static void control();
};

#endif /* ARTSS_SOLVER_DIFFUSIONSOLVER_H_ */
