/// \file       IPressure.h
/// \brief      Interface for pressure method
/// \date       Sep 14, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACES_IPRESSURE_H_
#define ARTSS_INTERFACES_IPRESSURE_H_

#include "../Field.h"

class IPressure {
public:
    IPressure() = default;
    virtual ~IPressure() = default;
    virtual void pressure(Field *out, Field *b, real t, bool sync) = 0;

    void divergence(Field *out, const Field *in_x, const Field *in_y, const Field *in_z, bool sync);
    void projection(Field *out_u, Field *out_v, Field *out_w, const Field *in_u, const Field *in_v, const Field *in_w, const Field *in_p, bool sync);
};

#endif /* ARTSS_INTERFACES_IPRESSURE_H_ */
