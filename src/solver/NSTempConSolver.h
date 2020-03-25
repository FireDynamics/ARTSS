/// \file 		Navier-Stokes Solver with force f(T)
/// \brief 		Defines the steps to solve advection, diffusion, pressure and add sources (dependent on T), solves temperature and concentration equation
/// \date 		Sep 27, 2017
/// \author 	Küsters
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_SOLVER_NSTEMPCONSOLVER_H_
#define ARTSS_SOLVER_NSTEMPCONSOLVER_H_


#include "../interfaces/SolverI.h"
#include "../interfaces/AdvectionI.h"
#include "../interfaces/DiffusionI.h"
#include "../interfaces/PressureI.h"
#include "../interfaces/SourceI.h"
#include "../utility/GlobalMacrosTypes.h"

class NSTempConSolver:public SolverI {
public:
	NSTempConSolver();
	~NSTempConSolver() override;

	void DoStep(real t, bool sync) override;

private:
	AdvectionI* adv_vel;
	DiffusionI* dif_vel;
	AdvectionI* adv_temp;
	DiffusionI* dif_temp;
	AdvectionI* adv_con;
	DiffusionI* dif_con;
	PressureI* pres;
	SourceI* sou_vel;
	SourceI* sou_temp;
	SourceI* sou_con;

	real m_nu;
	real m_kappa;
	real m_gamma;
	std::string m_dir_vel;

  static void control();

  std::string m_forceFct;
  bool m_hasDissipation;
  std::string m_tempFct;
  std::string m_conFct;
};

#endif /* ARTSS_SOLVER_NSTEMPCONSOLVER_H_ */
