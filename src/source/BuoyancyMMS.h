/// \file       BuoyancyMMS.h
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_SOURCE_BUOYANCYMMS_H_
#define ARTSS_SOURCE_BUOYANCYMMS_H_


#include "../interfaces/ISourceFunction.h"

class BuoyancyMMS: public ISourceFunction {
public:
    BuoyancyMMS();
    ~BuoyancyMMS();
    void update_source(Field *out, real t_cur) override;
    void set_noise(bool has_noise) { m_has_noise = has_noise; }
private:
    void set_up();
    Field *m_source_field;
    bool m_has_noise = false;
};


#endif /* ARTSS_SOURCE_BUOYANCYMMS_H_ */
