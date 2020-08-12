/// \file 		NSTempTurbConSolver.h
/// \brief 		Defines the (fractional) steps to solve the incompressible Navier-Stokes equations with force f(T), turbulence and concentration C
/// \date 		Oct 02, 2017
/// \author 	Küsters
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_SOLVER_NSTEMPTURBCONSOLVER_H_
#define ARTSS_SOLVER_NSTEMPTURBCONSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/ISource.h"
#include "../interfaces/ITurbulence.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../field/FieldController.h"

class NSTempTurbConSolver : public ISolver {
public:
    NSTempTurbConSolver(FieldController *field_controller);
    ~NSTempTurbConSolver();

    void do_step(real t, bool sync) override;

private:
    FieldController *m_field_controller;

    IAdvection *adv_vel;
    IDiffusion *dif_vel;
    IAdvection *adv_temp;
    IDiffusion *dif_temp;
    IAdvection *adv_con;
    IDiffusion *dif_con;
    IPressure *pres;
    ISource *sou_vel;
    ISource *sou_temp;
    ISource *sou_con;
    ITurbulence *mu_tub;

    real m_nu;
    real m_kappa;
    real m_gamma;
    std::string m_dir_vel = "";

    static void control();

    bool m_has_turbulence_temperature;
    bool m_has_turbulence_concentration;
    bool m_has_dissipation;
    std::string m_forceFct;
    std::string m_tempFct;
    std::string m_conFct;
};

#endif /* ARTSS_SOLVER_NSTEMPTURBCONSOLVER_H_ */
