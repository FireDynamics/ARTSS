/// \file       ISource.h
/// \brief      Interface for adding sources
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACES_ISOURCE_H_
#define ARTSS_INTERFACES_ISOURCE_H_

#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"

class ISource {
 public:
    virtual ~ISource() = default;

    virtual void add_source(Field &out_x, Field &out_y, Field &out_z,
            Field const &S_x, Field const &S_y, Field const &S_z, bool sync) = 0;
    virtual void add_source(Field &out, Field const &S, bool sync) = 0;

    void buoyancy_force(Field &out,
            Field const &in, Field const &in_temperature_ambient, bool sync = true);
    void dissipate(Field &out,
            Field const &in_u, Field const &in_v, Field const &in_w, bool sync = true);
};

#endif /* ARTSS_INTERFACES_ISOURCE_H_ */

