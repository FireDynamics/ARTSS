/// \file       BuoyancyMMS.h
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_SOURCE_BUOYANCYMMS_H_
#define ARTSS_SOURCE_BUOYANCYMMS_H_


#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"
#include "../interfaces/ISourceFunction.h"

class BuoyancyMMS: public ISourceFunction {
 public:
    BuoyancyMMS();
    ~BuoyancyMMS() = default;
    void update_source(Field &out, real t_cur) override;
 private:
    void set_up();
    Field m_source_field;
};


#endif /* ARTSS_SOURCE_BUOYANCYMMS_H_ */
