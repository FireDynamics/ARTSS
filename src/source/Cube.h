/// \file       Cube.h
/// \brief      
/// \date       Sep 29, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_CUBE_H
#define ARTSS_CUBE_H


#include "../interfaces/ISourceFunction.h"
#include <random>

class Cube: public ISourceFunction {
public:
    Cube(real value, real x0, real y0, real z0, real sigma_x, real sigma_y, real sigma_z);
    ~Cube();
    void update_source(Field *out, real t_cur) override;
    void set_noise(real range, int seed, real step_size) override;
private:
    Field *m_source_field;
    int m_has_noise = false;
    std::mt19937 m_mt;
    double m_steps;
    real m_step_size;
    void set_up(real value, real x0, real y0, real z0, real sigma_x, real sigma_y, real sigma_z);
};


#endif /* ARTSS_CUBE_H */
