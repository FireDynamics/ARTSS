/// \file       ISourceFunction.h
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_INTERFACES_ISOURCEFUNCTION_H
#define ARTSS_INTERFACES_ISOURCEFUNCTION_H

#include "IRandomField.h"
#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"
#include <string>

class ISourceFunction {
 public:
    void set_noise(IRandomField *noise_maker) {
        m_has_noise = true;
        m_noise_maker = noise_maker;
    }

    virtual void update_source(Field &out, real t_cur) = 0;

 protected:
    bool m_has_noise = false;
    IRandomField *m_noise_maker;
};


#endif /* ARTSS_INTERFACES_ISOURCEFUNCTION_H */
