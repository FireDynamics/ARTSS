/// \file       ITurbulence.h
/// \brief      Interface for Turbulence method
/// \date       August 18, 2016
/// \author     Suryanarayana Maddu
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACES_ITURBULENCE_H_
#define ARTSS_INTERFACES_ITURBULENCE_H_

#include "../field/Field.h"

class ITurbulence {

public:
    ITurbulence() = default;
    virtual ~ITurbulence() = default;

    virtual void calc_turbulent_viscosity(
            Field &ev,
            const Field &in_u, const Field &in_v, const Field &in_w,
            bool sync) = 0;
    virtual void explicit_filtering(Field &out, const Field &in, bool sync) = 0;
};

#endif /* ARTSS_INTERFACES_ITURBULENCE_H_ */
