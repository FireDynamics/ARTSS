/// \file       JacobiDiffuse.h
/// \brief      Solves diffusion equation with Jacobian method
/// \details    Solves Diffusion equation \f$ \partial_t \phi_2 = \nu \nabla^2 \phi_2 \f$ via calculated iterations of Jacobi step (dependent on residual/ maximal 1000)
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DIFFUSION_JACOBIDIFFUSE_H_
#define ARTSS_DIFFUSION_JACOBIDIFFUSE_H_

#include <memory>
#include "../interfaces/IDiffusion.h"
#include "../field/Field.h"
#include "../utility/Utility.h"

class JacobiDiffuse : public IDiffusion {
 public:
    JacobiDiffuse();

    void diffuse(Field &out, Field &in, Field const &b,
            real const D, bool sync) override;
    void diffuse(Field &out, Field &in, Field const &b,
            real const D, Field const &EV, bool sync) override;  // turbulent version

    static void JacobiStep(Field &out, Field const &in, Field const &b,
            real const alphaX, real const alphaY, real const alphaZ,
            real const beta, real const dsign, real const w, bool sync = true);
    static void JacobiStep(size_t level, Field &out, Field const &in, Field const &b,
            real const alphaX, real const alphaY, real const alphaZ,
            real const beta, real const dsign, real const w, bool sync = true);  // Multigrid version
    static void JacobiStep(Field &out, Field const &in, Field const &b,
            real const dsign, real const w, real const D,
            Field const &EV, real dt, bool sync = true);  // turbulent version

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

#endif /* ARTSS_DIFFUSION_JACOBIDIFFUSE_H_ */
