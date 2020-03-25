/// \file 		DiffusionTurbSolver.h
/// \brief 		Defines the steps to solve the turbulent diffusion equation
/// \date 		August 18, 2016
/// \author 	Suryanarayana Maddu
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_TURBDIFFUSIONSOLVER_H_
#define ARTSS_SOLVER_TURBDIFFUSIONSOLVER_H_

#include "../interfaces/SolverI.h"
#include "../interfaces/DiffusionI.h"
#include "../interfaces/TurbulenceI.h"
#include "../utility/GlobalMacrosTypes.h"

class DiffusionTurbSolver: public SolverI {
public:
	DiffusionTurbSolver();
	~DiffusionTurbSolver() override;

	void DoStep(real t, bool sync) override;

private:
	DiffusionI*   dif;
	TurbulenceI*  mu_tub;

	real m_nu;

    static void control();
};

#endif /* ARTSS_SOLVER_TURBDIFFUSIONSOLVER_H_ */
