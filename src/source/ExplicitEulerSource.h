/// \file       ExplicitSource.h
/// \brief      Adding source via Explicit Euler
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOURCE_EXPLICITEULERSOURCE_H_
#define ARTSS_SOURCE_EXPLICITEULERSOURCE_H_

#include "../interfaces/ISource.h"
#include "../Field.h"
#include "../utility/GlobalMacrosTypes.h"

class ExplicitEulerSource: public ISource {
public:
  ExplicitEulerSource();

  void add_source(Field* out_x, Field* out_y, Field* out_z, Field* S_x, Field* S_y, Field* S_z, bool sync) override;
  void add_source(Field* out, Field* S, bool sync) override;

private:
  real m_dt;
  std::string m_dir_vel ="";
};

#endif /* ARTSS_SOURCE_EXPLICITEULERSOURCE_H_ */
