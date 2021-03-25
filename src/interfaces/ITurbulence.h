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

    virtual void CalcTurbViscosity(
            Field &ev,
            Field const &in_u, Field const &in_v, Field const &in_w, bool sync) = 0;
    virtual void ExplicitFiltering(Field &out, Field const &in, bool sync) = 0;
};

#endif /* ARTSS_INTERFACES_ITURBULENCE_H_ */
