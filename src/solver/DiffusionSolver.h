/// \file       DiffusionSolver.h
/// \brief      Defines the steps to solve the diffusion equation
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOLVER_DIFFUSIONSOLVER_H_
#define ARTSS_SOLVER_DIFFUSIONSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IDiffusion.h"
#include "../field/FieldController.h"

#ifndef BENCHMARKING
#include "../utility/Utility.h"
#endif

class DiffusionSolver: public ISolver {
 public:
    DiffusionSolver(FieldController *field_controller);
    ~DiffusionSolver();

    void do_step(real t, bool sync) override;

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    FieldController *m_field_controller;
    IDiffusion *dif;
    real m_nu;

    static void control();
};

#endif /* ARTSS_SOLVER_DIFFUSIONSOLVER_H_ */
