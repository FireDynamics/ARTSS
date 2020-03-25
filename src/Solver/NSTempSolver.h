///	\file 		Navier-Stokes Solver with force f(T)
/// \brief 		Defines the steps to solve advection, diffusion, pressure and add sources (dependent on T), solves temperature equation
/// \date 		Dec 2, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_NSTEMPSOLVER_H_
#define ARTSS_SOLVER_NSTEMPSOLVER_H_


#include "../Interfaces/SolverI.h"
#include "../Interfaces/AdvectionI.h"
#include "../Interfaces/DiffusionI.h"
#include "../Interfaces/PressureI.h"
#include "../Interfaces/SourceI.h"
#include "../Utility/GlobalMacrosTypes.h"

class NSTempSolver:public SolverI {
public:
	NSTempSolver();
	~NSTempSolver() override;

	void DoStep(real t, bool sync) override;

private:
	AdvectionI* adv_vel;
	DiffusionI* dif_vel;
	AdvectionI* adv_temp;
	DiffusionI* dif_temp;
	PressureI* pres;
	SourceI* sou_vel;
	SourceI* sou_temp;

	real m_nu;
	real m_kappa;
	std::string m_dir_vel = "";

  static void control();

  std::string m_forceFct;
  bool m_hasDissipation;
  std::string m_tempFct;
};

#endif /* ARTSS_SOLVER_NSTEMPSOLVER_H_ */
