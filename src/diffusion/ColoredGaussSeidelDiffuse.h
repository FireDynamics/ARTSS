/// \file         ColoredGaussSeidelDiffuse.h
/// \brief        solves Diffusion equation \f$ \partial_t \phi_2 = \nu \nabla^2 \phi_2 \f$
///               via calculated iterations of CGS step (dependent on residual/ maximal 1000)
/// \date         June 21, 2016
/// \author       My Linh Wuerzburger
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H
#define ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H

#include "../interfaces/IDiffusion.h"
#include "../field/Field.h"
#include "../utility/Utility.h"

class ColoredGaussSeidelDiffuse: public IDiffusion {
 public:
    explicit ColoredGaussSeidelDiffuse(Settings const &settings);

    void diffuse(
            Field &out, const Field &in, const Field &b,
            real D, bool sync = true) override;
    void diffuse(
            Field &out, const Field &in, const Field &b,
            real D, const Field &EV, bool sync = true) override;

    static void colored_gauss_seidel_step(
            Field &out, const Field &b,
            real alpha_x, real alpha_y, real alpha_z,
            real beta, real dsign, real w, bool sync = true);

    static void colored_gauss_seidel_step(
            Field &out, const Field &b,
            real dsign, real w, real D,
            const Field &EV, real dt, bool sync = true);  // turbulent version

    static void colored_gauss_seidel_stencil(
            size_t i, size_t j, size_t k,
            real *out, const real *b,
            real alpha_x, real alpha_y, real alpha_z,
            real dsign, real beta, real w,
            size_t Nx, size_t Ny);

 private:
    Settings const &m_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    real m_dsign;
    real m_w;
    size_t m_max_iter;
    real m_tol_res;
};

#endif /* ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H */
