/// \file       NSSolver.h
/// \brief      Defines the (fractional) steps to solve the incompressible Navier-Stokes equations
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved
#ifndef ARTSS_SOLVER_NSSOLVER_H_
#define ARTSS_SOLVER_NSSOLVER_H_

#include "../interfaces/ISolver.h"
#include "../interfaces/IAdvection.h"
#include "../interfaces/IDiffusion.h"
#include "../interfaces/IPressure.h"
#include "../interfaces/ISource.h"
#include "../utility/GlobalMacrosTypes.h"

class NSSolver : public ISolver {
public:
  NSSolver();

  ~NSSolver() override;

  void do_step(real t, bool sync) override;

private:
  IAdvection *adv_vel;
  IDiffusion *dif_vel;
  IPressure *pres;
  ISource *sou;

  real m_nu;

  static void control();

  std::string m_sourceFct;
};

#endif /* ARTSS_SOLVER_NSSOLVER_H_ */
