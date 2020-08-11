/// \file       PressureSolver.h
/// \brief      Defines the steps to solve the pressure Poisson equation
/// \date       Sep 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_PRESSURESOLVER_H_
#define ARTSS_SOLVER_PRESSURESOLVER_H_

#include "../pressure/VCycleMG.h"
#include "../utility/Parameters.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../Domain.h"
#include "SolverSelection.h"
#include "../boundary/BoundaryData.h"

#include "../interfaces/ISolver.h"
#include "../interfaces/IPressure.h"

class PressureSolver : public ISolver {
public:
    PressureSolver();
    ~PressureSolver() override;
    void do_step(real t, bool sync) override;

private:
    IPressure *pres;
    static void control();
};

#endif /* ARTSS_SOLVER_PRESSURESOLVER_H_ */
