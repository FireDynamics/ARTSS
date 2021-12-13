/// \file       Functions.h
/// \brief      Functions for Initialization
/// \date       Jun 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_FUNCTIONS_H_
#define ARTSS_FUNCTIONS_H_

#include <string>

#include "DomainData.h"
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
    static const std::string mcdermott;
    static const std::string sin_sin_sin;
    static const std::string uniform;
    static const std::string vortex;
    static const std::string vortex_y;
    static const std::string zero;
};

namespace Functions {  // alphabetically ordered

  void beltrami(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
          real t, real a, real d, real nu);
  void beltrami_bc_p(Field &out_x, real a);
  void beltrami_bc_u(Field &out_x, real t, real a, real d, real nu);
  void beltrami_bc_v(Field &out_x, real t, real a, real d, real nu);
  void beltrami_bc_w(Field &out_x, real t, real a, real d, real nu);
  void buoyancy_force(Field &out, Field &T, Field &T_ambient, real beta, real g);
  void buoyancy_mms(Field &out_x, Field &out_y, Field &out_z, Field &out_p, Field &out_T,
          real t, real nu, real beta, real g, real rhoa);
  void buoyancy_st_mms(Field &out, real t, real nu, real beta, real kappa, real g, real rhoa);

  void drift(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
          real u_lin, real v_lin, real w_lin, real pa);

  void exp_sinus_prod(Field &out, real t, real nu, real l);
  void exp_sinus_sum(Field &out_x, Field &out_y, Field &out_z, real t, real nu);

  void fac_sin_sin_sin(Field &out, real l);

  void gauss_bubble(Field &out, real t,
          real u_lin, real v_lin, real w_lin,
          real x_shift, real y_shift, real z_shift,
          real l);

  void hat(Field &out,
          real start_x, real end_x,
          real start_y, real end_y,
          real start_z, real end_z,
          real val_in, real val_out);

  void jet(
          Field &out,
          size_t index_x1, size_t index_x2,
          size_t index_y1, size_t index_y2,
          size_t index_z1, size_t index_z2,
          real value);

  void layers(Field &out, int n_layers, CoordinateAxis axis,
              real *borders, const real *values);

  void mcdermott(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
          real t, real nu, real A);

  void random(Field &out,
          real range, bool is_absolute, int seed, real step_size);

  void sin_sin_sin(Field &out, real l);

  void uniform(Field &out, real val);

  void vortex(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
          real u_lin, real v_lin, real pa, real rhoa);
  void vortex_y(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
          real u_lin, real v_lin, real pa, real rhoa);
};

#endif /* ARTSS_FUNCTIONS_H_ */
