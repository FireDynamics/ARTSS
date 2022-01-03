/// \file       ExplicitSource.h
/// \brief      Adding source via Explicit Euler
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOURCE_EXPLICITEULERSOURCE_H_
#define ARTSS_SOURCE_EXPLICITEULERSOURCE_H_

#include "../interfaces/ISource.h"
#include "../field/Field.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"

class ExplicitEulerSource : public ISource {
 public:
    explicit ExplicitEulerSource(const std::vector<CoordinateAxis> &dir);

    void add_source(Field &out_x, Field &out_y, Field &out_z,
                    const Field &S_x, const Field &S_y, const Field &S_z,
                    bool sync) override;
    void add_source(Field &out, const Field &S, bool sync) override;

 private:
    const std::vector<CoordinateAxis> &m_dir;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    std::string m_dir_vel;
};

#endif /* ARTSS_SOURCE_EXPLICITEULERSOURCE_H_ */
