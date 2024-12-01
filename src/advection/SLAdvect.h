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

using return_backtracking_parameters  = std::tuple<size_t, size_t, real>;
class SLAdvect : public IAdvection {
 public:
    SLAdvect() = default;
    ~SLAdvect() override = default;

    void advect(Field &out, const Field &in, const Field &u_vel, const Field &v_vel, const Field &w_vel, bool sync) override;
    static return_backtracking_parameters calculate_backward_index(CoordinateAxis axis, size_t coordinate, real epsilon, real trace_back);
 private:
};

#endif /* ARTSS_ADVECTION_SLADVECT_H_ */

