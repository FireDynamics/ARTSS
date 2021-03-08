/// \file       ExplicitDiffuse.h
/// \brief      Solves diffusion equation with an explicit method
/// \date       December 12, 2019
/// \author     Max Boehler
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DIFFUSION_EXPLICITDIFFUSE_H_
#define ARTSS_DIFFUSION_EXPLICITDIFFUSE_H_

#include "../interfaces/IDiffusion.h"

#ifndef BENCHMARKING
#include <memory>
#include "../utility/Utility.h"
#endif

#include "../Domain.h"
#include "../field/Field.h"
#include "../utility/Parameters.h"
#include "../boundary/BoundaryController.h"

class ExplicitDiffuse : public IDiffusion {
 public:
    ExplicitDiffuse();
    // ExplicitDiffuse(real dt, const Domain &domain, const BoundaryController &boundary);

    void diffuse(Field *out, const Field &in, const Field &b,
            const Field &u, const Field &v, const Field &w,
            real D, bool sync) override;
    void diffuse(Field *out, const Field &in, const Field &b,
            const Field &u, const Field &v, const Field &w,
            real D, const Field &EV, bool sync) override;  // turbulent version

    void ExplicitStep(Field *out, const Field &in, real D, bool sync = true);
    void ExplicitStep(Field *out, const Field &in,
        const Field &u, const Field &v, const Field &w,
        const Field &EV, real D, real C, bool sync);

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    real m_dt;
    real m_cs;
    const Domain &m_domain;
    BoundaryController *m_boundary;
};

#endif /* ARTSS_DIFFUSION_EXPLICITDIFFUSE_H_ */
