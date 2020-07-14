/// \file       NSTurbSolver.h
/// \brief      Defines the steps to solve advection, diffusion, pressure and add sources (with LES turbulence)
/// \date       Feb 15, 2017
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_NSTURBSOLVER_H_
#define ARTSS_SOLVER_NSTURBSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/ISource.h"
#include "../interfaces/ITurbulence.h"
#include "../utility/GlobalMacrosTypes.h"

class NSTurbSolver : public ISolver {
public:
    NSTurbSolver();
    ~NSTurbSolver() override;

    void do_step(real t, bool sync) override;

private:
    IAdvection *adv_vel;
    IDiffusion *dif_vel;
    IPressure *pres;
    ISource *sou_vel;
    ITurbulence *mu_tub;

    real m_nu;

    static void control();

    std::string m_force_function;
};

#endif /* ARTSS_SOLVER_NSTURBSOLVER_H_ */
