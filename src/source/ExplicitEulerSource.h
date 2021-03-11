/// \file       ExplicitSource.h
/// \brief      Adding source via Explicit Euler
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOURCE_EXPLICITEULERSOURCE_H_
#define ARTSS_SOURCE_EXPLICITEULERSOURCE_H_

#include <memory>
#include <string>
#include "../interfaces/ISource.h"
#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"

class ExplicitEulerSource : public ISource {
 public:
    ExplicitEulerSource();

    void add_source(Field &out_x, Field &out_y, Field &out_z,
            Field const &s_x, Field const &s_y, Field const &s_z, bool sync) override;
    void add_source(Field &out, Field const &s, bool sync) override;

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    real m_dt;
    std::string m_dir_vel;
};

#endif /* ARTSS_SOURCE_EXPLICITEULERSOURCE_H_ */
