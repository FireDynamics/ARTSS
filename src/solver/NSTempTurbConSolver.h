/*
 * NSTempTurbConSolver.h
 *
 *  Created on: Oct 02, 2017
 *      Author: Küsters
 */

/// \file 		NSTempTurbConSolver.h
/// \brief 		Defines the (fractional) steps to solve the incompressible Navier-Stokes equations with force f(T), turbulence and concentration C
/// \date 		Sep 27, 2017
/// \author 	Küsters
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef NSTEMPTURBCONSOLVER_H_
#define NSTEMPTURBCONSOLVER_H_

#include "../interfaces/SolverI.h"
#include "../interfaces/AdvectionI.h"
#include "../interfaces/DiffusionI.h"
#include "../interfaces/PressureI.h"
#include "../interfaces/SourceI.h"
#include "../interfaces/TurbulenceI.h"
#include "../utility/GlobalMacrosTypes.h"

class NSTempTurbConSolver:public SolverI {
public:
	NSTempTurbConSolver();
	~NSTempTurbConSolver() override;

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
	TurbulenceI*  mu_tub;

	real m_nu;
	real m_kappa;
	real m_gamma;
	std::string m_dir_vel = "";

    static void control();

    bool m_hasTurbulenceTemperature;
    bool m_hasTurbulenceConcentration;
    bool m_hasDissipation;
    std::string m_forceFct;
    std::string m_tempFct;
    std::string m_conFct;
};

#endif /* NSTEMPTURBCONSOLVER_H_ */
