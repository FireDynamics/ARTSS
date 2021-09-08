/// \file       SLAdvect.h
/// \brief      Solves advection equation via unconditionally stable Semi-Lagrangian approach
/// \date       Aug 23, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ADVECTION_SLADVECT_H_
#define ARTSS_ADVECTION_SLADVECT_H_

#include "../interfaces/IAdvection.h"
#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"

class SLAdvect : public IAdvection {
 public:
    SLAdvect();
    ~SLAdvect() override = default;

    void advect(
            Field &out, Field const &in,
            Field const &u_vel, Field const &v_vel, Field const &w_vel,
            bool sync) override;

 private:
    real m_dt;
};

#endif /* ARTSS_ADVECTION_SLADVECT_H_ */

