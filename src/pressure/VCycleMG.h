/// \file       VCycleMG.h
/// \brief      Defines V-cycle of geometric multigrid method
/// \date       Sep 14, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_PRESSURE_VCYCLEMG_H_
#define ARTSS_PRESSURE_VCYCLEMG_H_

#include <vector>
#include <src/interfaces/IDiffusion.h>
#include "../interfaces/IPressure.h"
#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"



class VCycleMG: public IPressure{
public:
    VCycleMG(Field *out, Field *b);
    ~VCycleMG() override;

    void pressure(Field *out, Field *b, real t, bool sync) override;

private:
    void VCycleMultigrid(Field *field_out, bool sync = true);
    void UpdateInput(Field *out, Field *b, bool sync = true);
    void Smooth(Field *out, Field *tmp, Field *b, size_t level, bool sync = true);
    void Residuum(Field *out, Field *in, Field *b, size_t level, bool sync = true);
    void Restrict(Field *out, Field *in, size_t level, bool sync = true);
    static void Prolongate(Field *out, Field *in, size_t level, bool sync = true);
    void Solve( Field *out, Field *tmp, Field *b, size_t level, bool sync = true);

    void call_smooth_colored_gauss_seidel(Field *out, Field *tmp, Field *b, size_t level, bool sync);
    void call_smooth_jacobi(Field *out, Field *tmp, Field *b, size_t level, bool sync);
    void call_solve_colored_gauss_seidel(Field *out, Field *tmp, Field *b, size_t level, bool sync);
    void call_solve_jacobi(Field *out, Field *tmp, Field *b, size_t level, bool sync);
    void call_colored_gauss_seidel(Field *out, Field *tmp, Field *b, size_t level, bool sync);
    void call_jacobi(Field *out, Field *tmp, Field *b, size_t level, bool sync);

    size_t m_levels;
    int m_n_cycle;
    int m_n_relax;

    void (VCycleMG::*m_smooth_diffusion_function)(Field *, Field *, Field *, size_t, bool);
    void (VCycleMG::*m_solve_diffusion_function)(Field *, Field *, Field *, size_t, bool);
    size_t m_diffusion_max_iter;
    real m_diffusion_tol_res;

    real m_dt;
    real m_dsign;
    real m_w;

    /**
     *
     */
    std::vector<Field*> m_residuum0;
    /**
     * stores on level 0 rhs/b field
     */
    std::vector<Field*> m_residuum1;
    std::vector<Field*> m_error0;
    /**
     * stores on level 0 out/p field
     */
    std::vector<Field*> m_error1;
    std::vector<Field*> m_mg_temporal_solution; // only as storage? or even sometimes as output field?
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

};

#endif /* ARTSS_PRESSURE_VCYCLEMG_H_ */
