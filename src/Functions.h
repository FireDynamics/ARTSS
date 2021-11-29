/// \file       Functions.h
/// \brief      Functions for Initialization
/// \date       Jun 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_FUNCTIONS_H_
#define ARTSS_FUNCTIONS_H_

#include <string>

#include "Domain.h"
#include "field/Field.h"
#include "utility/settings/Settings.h"

struct FunctionNames{
    static const std::string Beltrami;
    static const std::string BuoyancyMMS;
    static const std::string BuoyancyST_MMS;
    static const std::string Drift;
    static const std::string ExpSinusProd;
    static const std::string ExpSinusSum;
    static const std::string GaussBubble;
    static const std::string Hat;
    static const std::string Jet;
    static const std::string McDermott;
    static const std::string RandomC;
    static const std::string SinSinSin;
    static const std::string Uniform;
    static const std::string Vortex;
    static const std::string VortexY;
    static const std::string Zero;
};

namespace Functions {  // alphabetically ordered

  void Beltrami(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
          real t, real a, real d, real nu);
  void BeltramiBC_p(Field &out_x, real a);
  void BeltramiBC_u(Field &out_x, real t, real a, real d, real nu);
  void BeltramiBC_v(Field &out_x, real t, real a, real d, real nu);
  void BeltramiBC_w(Field &out_x, real t, real a, real d, real nu);
  void BuoyancyForce(Field &out, Field &T, Field &T_ambient, real beta, real g);
  void BuoyancyMMS(Field &out_x, Field &out_y, Field &out_z, Field &out_p, Field &out_T,
          real t, real nu, real beta, real g, real rhoa);
  void BuoyancyST_MMS(Field &out, real t, real nu, real beta, real kappa, real g, real rhoa);

  void Drift(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
          real u_lin, real v_lin, real w_lin, real pa);

  void ExpSinusProd(Field &out, real t, real nu, real l);
  void ExpSinusSum(Field &out_x, Field &out_y, Field &out_z, real t, real nu);

  void FacSinSinSin(Field &out, real l);

  void GaussBubble(Field &out, real t,
          real u_lin, real v_lin, real w_lin,
          real x_shift, real y_shift, real z_shift,
          real l);

  void Hat(Field &out,
          real start_x, real end_x,
          real start_y, real end_y,
          real start_z, real end_z,
          real val_in, real val_out);

  void Jet(
          Field &out,
          size_t index_x1, size_t index_x2,
          size_t index_y1, size_t index_y2,
          size_t index_z1, size_t index_z2,
          real value);

  void Layers(std::string const log_level, std::string const log_file, Field &out,
              int n_layers, std::string const dir, real *borders, real *values);

  void McDermott(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
          real t, real nu, real A);

  void Random(Field &out,
          real range, bool is_absolute, int seed, real step_size);

  void SinSinSin(Field &out, real l);

  void Uniform(Field &out, real val);

  void Vortex(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
          real u_lin, real v_lin, real pa, real rhoa);
  void VortexY(Field &out_x, Field &out_y, Field &out_z, Field &out_p,
          real u_lin, real v_lin, real pa, real rhoa);
};

#endif /* ARTSS_FUNCTIONS_H_ */
