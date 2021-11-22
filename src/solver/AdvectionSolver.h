/// \file       AdvectionSolver.h
/// \brief      Defines the steps to solve the advection equation
/// \date       August 22, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_ADVECTIONSOLVER_H_
#define ARTSS_SOLVER_ADVECTIONSOLVER_H_

#include "../field/FieldController.h"
#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"

class AdvectionSolver : public ISolver {
 public:
    AdvectionSolver(
        Settings const &settings,
        FieldController *field_controller,
        real u_lin, real v_lin, real w_lin);
    AdvectionSolver(Settings const &settings, FieldController *field_controller);
    ~AdvectionSolver();

    void do_step(real t, bool sync) override;

 private:
    Settings const &m_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    FieldController *m_field_controller;

    IAdvection *adv;
    Field m_u_lin;
    Field m_v_lin;
    Field m_w_lin;

    void control();
};

#endif /* ARTSS_SOLVER_ADVECTIONSOLVER_H_ */
