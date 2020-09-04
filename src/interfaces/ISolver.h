/// \file       ISolver.h
/// \brief      Interface holds solvers for solving governing equations
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACES_ISOLVER_H_
#define ARTSS_INTERFACES_ISOLVER_H_

#include <string>

#include "../utility/Utility.h"
#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"
#include "ISource.h"

struct SolverTypes {
    inline static const std::string AdvectionSolver = "AdvectionSolver";
    inline static const std::string AdvectionDiffusionSolver = "AdvectionDiffusionSolver";
    inline static const std::string DiffusionSolver = "DiffusionSolver";
    inline static const std::string DiffusionTurbSolver = "DiffusionTurbSolver";
    inline static const std::string NSSolver = "NSSolver";
    inline static const std::string NSTempSolver = "NSTempSolver";
    inline static const std::string NSTempConSolver = "NSTempConSolver";
    inline static const std::string NSTempTurbSolver = "NSTempTurbSolver";
    inline static const std::string NSTempTurbConSolver = "NSTempTurbConSolver";
    inline static const std::string NSTurbSolver = "NSTurbSolver";
    inline static const std::string PressureSolver = "PressureSolver";
};

class ISolver {
 public:
    virtual ~ISolver() = default;
    virtual void do_step(real t, bool sync) = 0;
};

#endif /* ARTSS_INTERFACES_ISOLVER_H_ */
