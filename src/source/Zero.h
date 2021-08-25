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
    void read_header_part(std::string &header) override { }
    std::string write_header_part() override { return {}; }
};

#endif /* ARTSS_ZERO_H */
