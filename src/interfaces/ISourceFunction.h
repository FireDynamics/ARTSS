/// \file       ISourceFunction.h
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_INTERFACES_ISOURCEFUNCTION_H
#define ARTSS_INTERFACES_ISOURCEFUNCTION_H

#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"

class ISourceFunction {
public:
    virtual void update_source(Field *out, real t_cur) = 0;
    virtual void set_noise(real range, int seed, real step_size) { }
};


#endif /* ARTSS_INTERFACES_ISOURCEFUNCTION_H */
