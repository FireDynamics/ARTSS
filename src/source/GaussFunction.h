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
#include "../utility/settings/Settings.h"
#include "../utility/Utility.h"

class GaussFunction: public ISourceFunction {
 public:
    explicit GaussFunction(const Settings::solver::sources::gauss &settings);

    void update_source(Field &out, real t_cur) override;

 private:
    real get_time_value(real t_cur);
    void create_spatial_values();
    Field m_field_spatial_values;
    const Settings::solver::sources::gauss m_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_SOURCE_GAUSSFUNCTION_H_ */

