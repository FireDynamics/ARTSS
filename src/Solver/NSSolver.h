/// \file 		NSSolver.h
/// \brief 		Defines the (fractional) steps to solve the incompressible Navier-Stokes equations
/// \date 		Dec 2, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved
#ifndef ARTSS_SOLVER_NSSOLVER_H_
#define ARTSS_SOLVER_NSSOLVER_H_

#include "../Interfaces/SolverI.h"
#include "../Interfaces/AdvectionI.h"
#include "../Interfaces/DiffusionI.h"
#include "../Interfaces/PressureI.h"
#include "../Interfaces/SourceI.h"
#include "../Utility/GlobalMacrosTypes.h"

class NSSolver : public SolverI {
public:
	NSSolver();

	~NSSolver() override;

	void DoStep(real t, bool sync) override;

private:
	AdvectionI *adv_vel;
	DiffusionI *dif_vel;
	PressureI *pres;
	SourceI *sou;

	real m_nu;

	static void control();

	std::string m_sourceFct;
};

#endif /* ARTSS_SOLVER_NSSOLVER_H_ */
