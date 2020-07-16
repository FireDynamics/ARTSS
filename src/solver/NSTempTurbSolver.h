/// \file       NSTempTurbSolver.h
/// \brief      Defines the steps to solve advection, diffusion, pressure and add sources (dependent on T), solves temperature equation and turbulence
/// \date       Feb 15, 2017
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_NSTEMPTURBSOLVER_H_
#define ARTSS_SOLVER_NSTEMPTURBSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/ISource.h"
#include "../interfaces/ITurbulence.h"
#include "../utility/GlobalMacrosTypes.h"

class NSTempTurbSolver : public ISolver {
public:
    NSTempTurbSolver();
    ~NSTempTurbSolver() override;

    void do_step(real t, bool sync) override;

private:
    IAdvection *adv_vel;
    IDiffusion *dif_vel;
    IAdvection *adv_temp;
    IDiffusion *dif_temp;
    IPressure *pres;
    ISource *sou_vel;
    ISource *sou_temp;
    ITurbulence *mu_tub;

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
