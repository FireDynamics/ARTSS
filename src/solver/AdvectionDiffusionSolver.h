/// \file     AdvectionDiffusionSolver.h
/// \brief    Defines the steps to solve the advection and diffusion equation
/// \date     May 20, 2016
/// \author   Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_
#define ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_

#include "../field/FieldController.h"
#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"

class AdvectionDiffusionSolver : public ISolver {
public:
    AdvectionDiffusionSolver(Settings const &settings, FieldController *field_controller);
    ~AdvectionDiffusionSolver();

    void do_step(real t, bool sync) override;

private:
    Settings const &m_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    IAdvection *adv;
    IDiffusion *dif;

    FieldController *m_field_controller;

    real m_nu;

    void control();
};

#endif /* ARTSS_SOLVER_ADVECTIONDIFFUSIONSOLVER_H_ */
