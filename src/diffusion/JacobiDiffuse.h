/// \file       JacobiDiffuse.h
/// \brief      Solves diffusion equation with Jacobian method
/// \details    Solves Diffusion equation \f$ \partial_t \phi_2 = \nu \nabla^2 \phi_2 \f$ via calculated iterations of Jacobi step (dependent on residual/ maximal 1000)
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DIFFUSION_JACOBIDIFFUSE_H_
#define ARTSS_DIFFUSION_JACOBIDIFFUSE_H_

#include "../interfaces/IDiffusion.h"
#include "../Field.h"

#include "../utility/Utility.h"

class JacobiDiffuse : public IDiffusion {
public:
    JacobiDiffuse();

    void diffuse(Field *out, Field *in, const Field *b, real D, bool sync) override;
    void diffuse(Field *out, Field *in, const Field *b, real D, const Field *EV, bool sync) override;  // turbulent version

    static void JacobiStep(Field *out, const Field *in, const Field *b, real alphaX, real alphaY, real alphaZ, real beta, real dsign, real w, bool sync = true);
    static void JacobiStep(size_t level, Field *out, const Field *in, const Field *b, real alphaX, real alphaY, real alphaZ, real beta, real dsign, real w, bool sync = true); // Multigrid version
    static void JacobiStep(Field *out, const Field *in, const Field *b, real dsign, real w, real D, const Field *EV, real dt, bool sync = true); // turbulent version

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
