/// \file       GaussFunction.h
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_SOURCE_GAUSSFUNCTION_H_
#define ARTSS_SOURCE_GAUSSFUNCTION_H_

#include <math.h>
#include <memory>

#include "../field/Field.h"
#include "../boundary/Obstacle.h"
#include "../interfaces/ISourceFunction.h"
#include "../utility/GlobalMacrosTypes.h"

class GaussFunction: public ISourceFunction {
 public:
    GaussFunction(
            real HRR, real cp,
            real x0, real y0, real z0,
            real sigma_x, real sigma_y, real sigma_z,
            real tau);

    GaussFunction(real HRR, real cp);
    ~GaussFunction();

    void update_source(Field &out, real t_cur) override;

 private:
    real get_time_value(real t_cur);
    void create_spatial_values();
    Field m_field_spatial_values;
    real m_tau;
    real m_HRR, m_cp;
    real m_x0, m_y0, m_z0;
    real m_sigma_x, m_sigma_y, m_sigma_z;
};

#endif /* ARTSS_SOURCE_GAUSSFUNCTION_H_ */

