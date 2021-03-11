/// \file       GaussFunction.h
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_GAUSSFUNCTION_H
#define ARTSS_GAUSSFUNCTION_H


#include "../field/Field.h"
#include "../interfaces/ISourceFunction.h"

class GaussFunction: public ISourceFunction {
public:
    GaussFunction(real HRR, real cp, real x0, real y0, real z0, real sigma_x, real sigma_y, real sigma_z, real tau);
    ~GaussFunction();
    void update_source(Field &out, real t_cur) override;
private:
    void create_spatial_values(real HRR, real cp, real x0, real y0, real z0, real sigma_x, real sigma_y, real sigma_z);
    Field m_field_spatial_values;
    real m_tau;
    real get_time_value(real t_cur);
};


#endif /* ARTSS_GAUSSFUNCTION_H */
