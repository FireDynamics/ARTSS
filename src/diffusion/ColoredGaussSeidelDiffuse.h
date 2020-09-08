/// \file         ColoredGaussSeidelDiffuse.h
/// \brief        solves Diffusion equation \f$ \partial_t \phi_2 = \nu \nabla^2 \phi_2 \f$
///               via calculated iterations of CGS step (dependent on residual/ maximal 1000)
/// \date         June 21, 2016
/// \author       My Linh Wuerzburger
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H
#define ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H

#include "../interfaces/IDiffusion.h"
#include "../Field.h"
#include "../utility/Utility.h"

class ColoredGaussSeidelDiffuse: public IDiffusion {
 public:
    ColoredGaussSeidelDiffuse();

    void diffuse(Field* out, Field* in, const Field* b, const real D, bool sync = true);
    void diffuse(Field* out, Field* in, const Field* b, const real D, const Field* EV, bool sync = true);  // turbulent version
    static void colored_gauss_seidel_step(Field* out, const Field* b, const real alpha_x, const real alpha_y, const real alpha_z, const real beta, const real dsign, const real w, bool sync = true);
    static void colored_gauss_seidel_step(Field* out, const Field* b, const real dsign, const real w, const real D, const Field* EV, const real dt, bool sync = true); // turbulent version
    static void colored_gauss_seidel_stencil(size_t i, size_t j, size_t k, real* out, real* b, const real alpha_x, const real alpha_y, const real alpha_z, const real dsign, const real beta, const real w, const size_t Nx, const size_t Ny);

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    real m_dt;
    real m_dsign;
    real m_w;
    size_t m_max_iter;
    real m_tol_res;
};

#endif /* ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H */
