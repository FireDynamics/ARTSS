/// \file       Zero.h
/// \brief      
/// \date       Sep 09, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_SOURCE_ZERO_H
#define ARTSS_SOURCE_ZERO_H


#include "../interfaces/ISourceFunction.h"

class Zero: public ISourceFunction {
public:
    void update_source(Field*, real) override { }
    void set_noise(real, int, real) override { };
};


#endif /* ARTSS_SOURCE_ZERO_H */
