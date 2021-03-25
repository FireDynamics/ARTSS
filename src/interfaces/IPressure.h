/// \file       IPressure.h
/// \brief      Interface for pressure method
/// \date       Sep 14, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACES_IPRESSURE_H_
#define ARTSS_INTERFACES_IPRESSURE_H_

#include "../field/Field.h"
#include "ISolver.h"

class IPressure {
 public:
    IPressure() = default;
    virtual ~IPressure() = default;
    virtual void pressure(Field &out, Field const &b, real t, bool sync) = 0;

    void divergence(
            Field &out,
            Field const &in_x, Field const &in_y, Field const &in_z, bool sync);
    void projection(
            Field &out_u, Field &out_v, Field &out_w,
            Field const &in_u, Field const &in_v, Field const &in_w,
            Field const &in_p, bool sync);
};

#endif /* ARTSS_INTERFACES_IPRESSURE_H_ */
