/// \file       VCycleMG.h
/// \brief      Defines V-cycle of geometric multigrid method
/// \date       Sep 14, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_PRESSURE_VCYCLEMG_H_
#define ARTSS_PRESSURE_VCYCLEMG_H_

#include "../interfaces/IPressure.h"
#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"



class VCycleMG: public IPressure{
 public:
    VCycleMG(Field const &out, Field const &b);
    ~VCycleMG() override;

    void pressure(Field &out, Field const &b, real t, bool sync) override;

 private:
    void VCycleMultigrid(Field &out, bool sync = true);
    void UpdateInput(Field &out, Field const &b, bool sync = true);
    void Smooth(Field &out, Field &tmp, Field const &b, size_t level, bool sync = true);
    void Residuum(Field &out, Field const &in, Field const &b, size_t level, bool sync = true);
    void Restrict(Field &out, Field const &in, size_t level, bool sync = true);
    void Prolongate(Field &out, Field const &in, size_t level, bool sync = true);
    void Solve(Field &out, Field &tmp, Field const &b, size_t level, bool sync = true);

    void call_smooth_colored_gauss_seidel(Field &out, Field &tmp, Field const &b, size_t level, bool sync);
    void call_smooth_jacobi(Field &out, Field &tmp, Field const &b, size_t level, bool sync);
    void call_solve_colored_gauss_seidel(Field &out, Field &tmp, Field const &b, size_t level, bool sync);
    void call_solve_jacobi(Field &out, Field &tmp, Field const &b, size_t level, bool sync);

    const size_t m_levels;
    const int m_n_cycle;
    const int m_n_relax;

    void (VCycleMG::*m_smooth_diffusion_function)(Field &, Field &, Field const &, const size_t, bool);
    void (VCycleMG::*m_solve_diffusion_function)(Field &, Field &, Field const &, const size_t, bool);
    size_t m_diffusion_max_iter;
    real m_diffusion_tol_res;

    real m_dt;
    real m_dsign;
    real m_w;

    Field **m_residuum0;
    /**
     * stores on level 0 rhs/b field
     */
    Field **m_residuum1;
    Field **m_error0;
    /**
     * stores on level 0 out/p field
     */
    Field **m_error1;
    Field **m_mg_temporal_solution; // only as storage? or even sometimes as output field?
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_PRESSURE_VCYCLEMG_H_ */
