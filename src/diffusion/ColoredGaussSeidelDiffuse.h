/// \file         ColoredGaussSeidelDiffuse.h
/// \brief        solves Diffusion equation \f$ \partial_t \phi_2 = \nu \nabla^2 \phi_2 \f$
///               via calculated iterations of CGS step (dependent on residual/ maximal 1000)
/// \date         June 21, 2016
/// \author       My Linh Wuerzburger
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H
#define ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H

#include <memory>
#include "../interfaces/IDiffusion.h"
#include "../field/Field.h"
#include "../utility/Utility.h"

class ColoredGaussSeidelDiffuse: public IDiffusion {
 public:
    ColoredGaussSeidelDiffuse();

    void diffuse(Field &out, Field &in, Field const &b,
            const real D, bool sync = true) override;
    void diffuse(Field &out, Field &in, Field const &b,
            const real D, Field const &EV, bool sync = true) override;

    static void colored_gauss_seidel_step(Field &out, Field const &b,
            real alpha_x, real alpha_y, real alpha_z,
            real beta, real dsign, real w, bool sync = true);

    static void colored_gauss_seidel_step(Field &out, Field const &b,
            real const dsign, real const w, real const D,
            Field const &EV, const real dt, bool sync = true);  // turbulent version

    static void colored_gauss_seidel_stencil(size_t i, size_t j, size_t k,
            real *out, real *b,
            real const alpha_x, real const alpha_y, real const alpha_z,
            real const dsign, real const beta, real const w,
            size_t const Nx, size_t const Ny);

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
