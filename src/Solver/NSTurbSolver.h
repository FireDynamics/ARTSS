/// \file 		NSTurbSolver.h
/// \brief 		Defines the steps to solve advection, diffusion, pressure and add sources (with LES turbulence)
/// \date 		Feb 15, 2017
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_NSTURBSOLVER_H_
#define ARTSS_SOLVER_NSTURBSOLVER_H_

#include "../Interfaces/SolverI.h"
#include "../Interfaces/AdvectionI.h"
#include "../Interfaces/DiffusionI.h"
#include "../Interfaces/PressureI.h"
#include "../Interfaces/SourceI.h"
#include "../Interfaces/TurbulenceI.h"
#include "../Utility/GlobalMacrosTypes.h"

class NSTurbSolver:public SolverI {
public:
	NSTurbSolver();
	~NSTurbSolver() override;

	void DoStep(real t, bool sync) override;

private:
	AdvectionI* adv_vel;
	DiffusionI* dif_vel;
	PressureI* pres;
	SourceI* sou_vel;
	TurbulenceI*  mu_tub;

	real m_nu;

  static void control();

  std::string m_forceFct;
};

#endif /* ARTSS_SOLVER_NSTURBSOLVER_H_ */
