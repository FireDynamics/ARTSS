/*
 * NSTempTurbSolver.h
 *
 *  Created on: Feb 15, 2017
 *      Author: Severt
 */

/// \file 		NSTempTurbSolver.h
/// \brief 		Defines the steps to solve advection, diffusion, pressure and add sources (dependent on T), solves temperature equation and turbulence
/// \date 		Feb 15, 2017
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.
///

#ifndef ARTSS_SOLVER_NSTEMPTURBSOLVER_H_
#define ARTSS_SOLVER_NSTEMPTURBSOLVER_H_

#include "../interfaces/SolverI.h"
#include "../interfaces/AdvectionI.h"
#include "../interfaces/DiffusionI.h"
#include "../interfaces/PressureI.h"
#include "../interfaces/SourceI.h"
#include "../interfaces/TurbulenceI.h"
#include "../utility/GlobalMacrosTypes.h"

class NSTempTurbSolver:public SolverI {
public:
	NSTempTurbSolver();
	~NSTempTurbSolver() override;

	void DoStep(real t, bool sync) override;

private:
	AdvectionI* adv_vel;
	DiffusionI* dif_vel;
	AdvectionI* adv_temp;
	DiffusionI* dif_temp;
	PressureI* pres;
	SourceI* sou_vel;
	SourceI* sou_temp;
	TurbulenceI*  mu_tub;

	real m_nu;
	real m_kappa;

	std::string m_dir_vel = "";

    static void control();

    bool m_hasTurbulence;
    bool m_hasDissipation;
    std::string m_forceFct;
    std::string m_tempFct;
};

#endif /* ARTSS_SOLVER_NSTEMPTURBSOLVER_H_ */
