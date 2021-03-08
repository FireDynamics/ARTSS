/// \file       IDiffusion.h
/// \brief      Interface for diffusion methods
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACES_IDIFFUSION_H_
#define ARTSS_INTERFACES_IDIFFUSION_H_

#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"

class IDiffusion {

public:
    virtual void diffuse(Field *out, const Field &in, const Field &b,
            const Field &u, const Field &v, const Field &w,
            real D, bool sync) = 0;
    virtual void diffuse(Field *out, const Field &in, const Field &b,
            const Field &u, const Field &v, const Field &w,
            real D, const Field &ev, bool sync) = 0;
};

#endif /* ARTSS_INTERFACES_IDIFFUSION_H_ */
