/// \file       SLAdvect.h
/// \brief      Solves advection equation via unconditionally stable Semi-Lagrangian approach
/// \date       Aug 23, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ADVECTION_SLADVECT_H_
#define ARTSS_ADVECTION_SLADVECT_H_

#include "../interfaces/IAdvection.h"
#include "../field/Field.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"

class SLAdvect : public IAdvection {
 public:
    SLAdvect() = default;
    ~SLAdvect() override = default;

    void advect(Field &out, const Field &in, const Field &u_vel, const Field &v_vel, const Field &w_vel, bool sync) override;

 private:
};

#endif /* ARTSS_ADVECTION_SLADVECT_H_ */

