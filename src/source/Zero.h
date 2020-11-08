/// \file       Zero.h
/// \brief      
/// \date       Sep 09, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_ZERO_H
#define ARTSS_ZERO_H


#include "../interfaces/ISourceFunction.h"

class Zero: public ISourceFunction {
public:
    void update_source(Field *out, real t_cur) override { }
};


#endif /* ARTSS_ZERO_H */
