/// \file       Cube.h
/// \brief      
/// \date       Sep 29, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_CUBE_H
#define ARTSS_CUBE_H


#include <random>

#include "../interfaces/ISourceFunction.h"

class Cube: public ISourceFunction {
 public:
    Cube(real value,
            real x_start, real y_start, real z_start,
            real x_end, real y_end, real z_end);
    ~Cube();
    void update_source(Field *out, real t_cur) override;

 private:
    Field *m_source_field;
    void set_up(real value,
            real x_start, real y_start, real z_start,
            real x_end, real y_end, real z_end);
};


#endif /* ARTSS_CUBE_H */
