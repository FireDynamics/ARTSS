/// \file       Functions.h
/// \brief      Functions for Initialization
/// \date       Jun 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_FUNCTIONS_H_
#define ARTSS_FUNCTIONS_H_

#include "field/Field.h"
#include "Domain.h"
#include <string>

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

  void Beltrami(Field &out_x, Field &out_y, Field &out_z, Field &out_p, real t);
  void BeltramiBC_p(Field &out_x);
  void BeltramiBC_u(Field &out_x, real t);
  void BeltramiBC_v(Field &out_x, real t);
  void BeltramiBC_w(Field &out_x, real t);
  void BuoyancyForce(Field &out, Field &T, Field &T_ambient);
  void BuoyancyMMS(Field &out_x, Field &out_y, Field &out_z, Field &out_p, Field &out_T, real t);
  void BuoyancyST_MMS(Field &out, real t);

  void Drift(Field &out_x, Field &out_y, Field &out_z, Field &out_p);

  void ExpSinusProd(Field &out, real t);
  void ExpSinusSum(Field &out_x, Field &out_y, Field &out_z, real t);

  void FacSinSinSin(Field &out);

  void GaussBubble(Field &out, real t);

  void Hat(Field &out);

  void Jet(
          Field &out,
          size_t index_x1, size_t index_x2,
          size_t index_y1, size_t index_y2,
          size_t index_z1, size_t index_z2,
          real value);

  void Layers(Field &out);

  void McDermott(Field &out_x, Field &out_y, Field &out_z, Field &out_p, real t);

  void Random(Field &out, real range, bool is_absolute, int seed, real step_size);

  void SinSinSin(Field &out);

  void Uniform(Field &out, real val);

  void Vortex(Field &out_x, Field &out_y, Field &out_z, Field &out_p);
  void VortexY(Field &out_x, Field &out_y, Field &out_z, Field &out_p);

  void Zero(Field &field, size_t *arr_idx, size_t arr_idx_size);
};

#endif /* ARTSS_FUNCTIONS_H_ */
