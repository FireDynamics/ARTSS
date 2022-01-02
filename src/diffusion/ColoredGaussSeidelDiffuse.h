/// \file         ColoredGaussSeidelDiffuse.h
/// \brief        solves Diffusion equation \f$ \partial_t \phi_2 = \nu \nabla^2 \phi_2 \f$
///               via calculated iterations of CGS step (dependent on residual/ maximal 1000)
/// \date         June 21, 2016
/// \author       My Linh Wuerzburger
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H
#define ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H

#include <vector>

#include "../field/Field.h"
#include "../interfaces/IDiffusion.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"

class ColoredGaussSeidelDiffuse: public IDiffusion {
 public:
    explicit ColoredGaussSeidelDiffuse(const Settings::solver::diffusion_solvers::colored_gauss_seidel &settings);

    void diffuse(
            Field &out, const Field &in, const Field &b,
            real D, bool sync) override;
    void diffuse(
            Field &out, const Field &in, const Field &b,
            real D, const Field &EV, bool sync) override;

    static void colored_gauss_seidel_step(
            Field &out, const Field &b,
            real alpha_x, real alpha_y, real alpha_z,
            real beta, real dsign, real w,
            const std::vector<size_t> &odd,
            const std::vector<size_t> &even,
            bool sync = true);

    static void colored_gauss_seidel_step(
            Field &out, const Field &b,
            real dsign, real w, real D,
            const Field &EV, real dt,
            const std::vector<size_t> &odd,
            const std::vector<size_t> &even,
            bool sync = true);  // turbulent version

    static void colored_gauss_seidel_stencil(
            size_t index,
            real *out, const real *b,
            real alpha_x, real alpha_y, real alpha_z,
            real dsign, real beta, real w,
            size_t Nx, size_t Ny);

    static void create_red_black_lists(size_t level, std::vector<size_t> &odd, std::vector<size_t> &even);
 private:
    const Settings::solver::diffusion_solvers::colored_gauss_seidel &m_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    real m_dsign;

    std::vector<size_t> even_indices;
    std::vector<size_t> odd_indices;
};

#endif /* ARTSS_DIFFUSION_COLOREDGAUSSSEIDEL_H */
