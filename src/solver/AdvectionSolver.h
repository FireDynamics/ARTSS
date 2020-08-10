/// \file       AdvectionSolver.h
/// \brief      Defines the steps to solve the advection equation
/// \date       August 22, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_ADVECTIONSOLVER_H_
#define ARTSS_SOLVER_ADVECTIONSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../utility/GlobalMacrosTypes.h"

class AdvectionSolver : public ISolver {
public:
    AdvectionSolver();
    ~AdvectionSolver() override;

    void do_step(real t, bool sync) override;

private:
    IAdvection* adv;
    Field* u_linm;
    Field* v_linm;
    Field* w_linm;

    static void control();
};

#endif /* ARTSS_SOLVER_ADVECTIONSOLVER_H_ */
