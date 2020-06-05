/// \file       PressureSolver.h
/// \brief      Defines the steps to solve the pressure Poisson equation
/// \date       Sep 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_PRESSURESOLVER_H_
#define ARTSS_SOLVER_PRESSURESOLVER_H_


#include <spdlog/spdlog.h>
#include <iostream>

#include "../pressure/VCycleMG.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "SolverSelection.h"
#include "../boundary/BoundaryData.h"

#include "../interfaces/SolverI.h"
#include "../interfaces/PressureI.h"
#include "../utility/GlobalMacrosTypes.h"

class PressureSolver: public SolverI {
 public:
    PressureSolver();
    ~PressureSolver() override;

    void DoStep(real t, bool sync) override;

 private:
    PressureI* pres;

    static void control();
};

#endif /* ARTSS_SOLVER_PRESSURESOLVER_H_ */
