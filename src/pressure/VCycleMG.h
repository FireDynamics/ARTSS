/// \file       VCycleMG.h
/// \brief      Defines V-cycle of geometric multigrid method
/// \date       Sep 14, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_PRESSURE_VCYCLEMG_H_
#define ARTSS_PRESSURE_VCYCLEMG_H_

#include "../interfaces/IPressure.h"
#include "../field/Field.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"


class VCycleMG: public IPressure{
 public:
    VCycleMG(Settings::Settings const &settings);
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

    Settings::Settings const &m_settings;
    const size_t m_levels;
    const int m_n_cycle;
    const int m_n_relax;

    /**
     * difference between smooth and solve function:
     * - solve function repeats the chosen algorithm (m_solve_function) until the residuum is
     *   smaller than the given residual number or the number of repetitions exceeds the set number
     *   of iterations
     * - smooth function repeats the chosen algorithm (m_smooth_function) for a set number of times
     *   (m_n_cycles)
     */
    void (VCycleMG::*m_smooth_function)(Field &, Field &, Field const &, const size_t, bool);
    void (VCycleMG::*m_solve_function)(Field &, Field &, Field const &, const size_t, bool);
    size_t m_diffusion_max_iter;
    real m_diffusion_tol_res;

    const real m_dsign = -1;
    const real m_w;

    /**
     * size: 0 .. m_levels
     * level's going coarse: residuum: result of error1 + residuum1
     */
    Field **m_residuum0;
    /**
     * stores on level 0 rhs/b field
     * size: 0 ... m_levels + 1
     */
    Field **m_residuum1;
    /**
     * size: 0 .. m_levels
     * level's going coarse: restrict: result from residuum0
     * level's going fine: prolongate: result from error1
     */
    Field **m_error0;
    /**
     * stores on level 0 out/p field
     * size: 0 ... m_levels + 1
     * level's going coarse: smooth: result of residuum1
     * level's going fine: correction through error0
     * level's going fine: solve/smooth: result from residuum1
     */
    Field **m_error1;
    // storage only
    Field **m_mg_temporal_solution;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_PRESSURE_VCYCLEMG_H_ */
