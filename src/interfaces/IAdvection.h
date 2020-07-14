/// \file       IAdvection.h
/// \brief      Interface for advection methods
/// \date       Aug 22, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACES_IADVECTION_H_
#define ARTSS_INTERFACES_IADVECTION_H_

#include "../Field.h"

class IAdvection {
public:
    IAdvection() = default;
    virtual ~IAdvection() = default;
    virtual void advect(Field *out, Field *in, const Field *u_vel, const Field *v_vel, const Field *w_vel, bool sync) = 0;
};

#endif /* ARTSS_INTERFACES_IADVECTION_H_ */
