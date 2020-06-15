/// \file         ColoredGaussSeidelDiffuse.h
/// \brief        solves Diffusion equation \f$ \partial_t \phi_2 = \nu \nabla^2 \phi_2 \f$
///               via calculated iterations of CGS step (dependent on residual/ maximal 1000)
/// \date         June 21, 2016
/// \author       lgewuerz
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H
#define ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H

#include "../interfaces/DiffusionI.h"
#include "../Field.h"

class ColoredGaussSeidelDiffuse: public DiffusionI {
public:
    ColoredGaussSeidelDiffuse();

    void diffuse(Field* out, Field* in, const Field* b, const real D, bool sync = true);
    void diffuse(Field* out, Field* in, const Field* b, const real D, const Field* EV, bool sync = true);  // turbulent version
    static void ColoredGaussSeidelStep(Field* out, const Field* b, const real alphaX, const real alphaY, const real alphaZ, const real beta, const real dsign, const real w, bool sync = true);
    static void ColoredGaussSeidelStep(Field* out, const Field* b, const real dsign, const real w, const real D, const Field* EV, const real dt, bool sync = true); // turbulent version
    static void ColoredGaussSeidelStencil(size_t i, size_t j, size_t k, real* out, real* b, const real alphaX, const real alphaY, const real alphaZ, const real dsign, const real beta, const real w, const size_t nx, const size_t ny);

private:
#ifndef PROFILING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    real m_dt;
    real m_dsign;
    real m_w;
    size_t m_max_iter;
    real m_tol_res;
};

#endif /* ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H */
