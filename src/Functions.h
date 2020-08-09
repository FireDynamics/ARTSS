/// \file       Functions.h
/// \brief      Functions for Initialization
/// \date       Jun 13, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_FUNCTIONS_H_
#define ARTSS_FUNCTIONS_H_

#include <string>
#include "Field.h"

struct FunctionNames{
    static const std::string Beltrami;
    static const std::string BuoyancyMMS;
    static const std::string BuoyancyST_MMS;
    static const std::string Drift;
    static const std::string ExpSinusProd;
    static const std::string ExpSinusSum;
    static const std::string GaussBubble;
    static const std::string Hat;
    static const std::string McDermott;
    static const std::string RandomC;
    static const std::string RampTanh;
    static const std::string SinSinSin;
    static const std::string Uniform;
    static const std::string Vortex;
    static const std::string VortexY;
    static const std::string Zero;
};

namespace Functions {  // alphabetically ordered

  void Beltrami(Field* outx, Field* outy, Field* outz, Field* outp, real t);
  void BeltramiBC_p(Field* outx);
  void BeltramiBC_u(Field* outx, real t);
  void BeltramiBC_v(Field* outx, real t);
  void BeltramiBC_w(Field* outx, real t);
  void BuoyancyForce(Field* out, Field* T, Field* Ta);
  void BuoyancyMMS(Field* outx, Field* outy, Field* outz, Field* outp, Field* outT, real t);
  void BuoyancyST_MMS(Field* out, real t);

  void Drift(Field* outx, Field* outy, Field* outz, Field* outp);

  void ExpSinusProd(Field* out, real t);
  void ExpSinusSum(Field* outx, Field* outy, Field* outz, real t);

  void FacSinSinSin(Field* out);

  void GaussBubble(Field* out, real t);

  void Hat(Field* out);

  void Layers(Field* out);

  void McDermott(Field* outx, Field* outy, Field* outz, Field* outp, real t);

  real RampTanh(real t);

  void Random(Field* out, real range, bool is_absolute, int seed, real step_size);

  void SinSinSin(Field* out);

  void Uniform(Field* out, real val);

  void Vortex(Field *outx, Field *outy, Field *outz, Field *outp);
  void VortexY(Field *outx, Field *outy, Field *outz, Field *outp);

  void Vortex(Field *outx, Field *outy, Field *outz, Field *outp);
  void VortexY(Field *outx, Field *outy, Field *outz, Field *outp);

  void Zero(Field* field, size_t* arr_idx, size_t arr_idx_size);
};

#endif /* ARTSS_FUNCTIONS_H_ */
