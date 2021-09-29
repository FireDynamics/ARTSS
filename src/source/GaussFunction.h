/// \file       GaussFunction.h
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_SOURCE_GAUSSFUNCTION_H_
#define ARTSS_SOURCE_GAUSSFUNCTION_H_


#include <string>
#include "../field/Field.h"
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
    void set_HRR(real HRR) { m_HRR = HRR; }
    void set_x0(real x0) { m_x0 = x0; }
    void set_y0(real y0) { m_y0 = y0; }
    void set_z0(real z0) { m_z0 = z0; }
    void set_sigma_x(real sigma_x) {  m_sigma_x = sigma_x; }
    void set_sigma_y(real sigma_y) {  m_sigma_y = sigma_y; }
    void set_sigma_z(real sigma_z) {  m_sigma_z = sigma_z; }

    ~GaussFunction();

    void update_source(Field &out, real t_cur) override;

    void read_header_part(std::string &header) override;
    std::string write_header_part() override;

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

