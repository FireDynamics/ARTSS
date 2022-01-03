/// \file       Functions.h
/// \brief      Functions for Initialization
/// \date       Jun 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_FUNCTIONS_H_
#define ARTSS_FUNCTIONS_H_

#include <string>

#include "domain/DomainData.h"
#include "field/Field.h"
#include "utility/settings/Settings.h"

struct FunctionNames{
    static const std::string beltrami;
    static const std::string buoyancy_mms;
    static const std::string buoyancy_st_mms;
    static const std::string drift;
    static const std::string exp_sinus_prod;
    static const std::string exp_sinus_sum;
    static const std::string gauss_bubble;
    static const std::string hat;
    static const std::string jet;
    static const std::string layers;
    static const std::string mcdermott;
    static const std::string sin_sin_sin;
    static const std::string uniform;
    static const std::string vortex;
    static const std::string vortex_y;
    static const std::string zero;
};

namespace Functions {  // alphabetically ordered

  void beltrami(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
                real t, const Settings::initial_conditions::beltrami &beltrami);
  void beltrami_bc_p(Field &out_x, real a);
  void beltrami_bc_u(Field &out_x, real t, real a, real d);
  void beltrami_bc_v(Field &out_x, real t, real a, real d);
  void beltrami_bc_w(Field &out_x, real t, real a, real d);
  void buoyancy_mms(Field &out_x, Field &out_y, Field &out_z, Field &out_p, Field &out_T, real t);
  void buoyancy_st_mms(Field &out, real t, real beta, real kappa, real g, real rhoa);

  void drift(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
             const Settings::initial_conditions::drift &drift);

  void exp_sinus_prod(Field &out, real t, const Settings::initial_conditions::exp_sinus_prod &exp_sinus_prod);
  void exp_sinus_sum(Field &out_x, Field &out_y, Field &out_z, real t);

  void fac_sin_sin_sin(Field &out, const Settings::initial_conditions::sin_sin_sin &sin_sin_sin);

  void gauss_bubble(Field &out, real t, const Settings::initial_conditions::gauss_bubble &gauss_bubble);

  void hat(Field &out, const Settings::initial_conditions::hat &hat);

  void jet(Field &out, const Settings::initial_conditions::jet &jet);

  void layers(Field &out, const Settings::initial_conditions::layers_temperature &layers);

  void mcdermott(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
          real t, const Settings::initial_conditions::mc_dermott &mc_dermott);

  void random(Field &out, const Settings::random_parameters &random_params);

  void sin_sin_sin(Field &out, const Settings::initial_conditions::sin_sin_sin &sin_sin_sin);

  void uniform(Field &out, const Settings::initial_conditions::uniform &uniform);

  void vortex(Field &out_x, Field &out_y, Field &out_z, Field &out_p, const Settings::initial_conditions::vortex &vortex);

  void vortex_y(Field &out_x, Field &out_y, Field &out_z, Field &out_p, const Settings::initial_conditions::vortex &vortex);
};

#endif /* ARTSS_FUNCTIONS_H_ */
