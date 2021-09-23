/// \file       Solution.h
/// \brief      Calculates analytical solution
/// \details    This class calculates the analytical solution of different test cases
/// \date       Jul 11, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ANALYSIS_SOLUTION_H_
#define ARTSS_ANALYSIS_SOLUTION_H_

#include "../field/Field.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"

class Solution {
 public:
    explicit Solution(const std::string &initial_condition, bool has_analytical_solution);

    void calc_analytical_solution(real t);

    // Getter
    real* GetU_data() const { return m_u_analytical_solution.data; }
    real* GetV_data() const { return m_v_analytical_solution.data; }
    real* GetW_data() const { return m_w_analytical_solution.data; }
    real* GetP_data() const { return m_p_analytical_solution.data; }
    real* GetT_data() const { return m_T_analytical_solution.data; }
    return_ptr GetU() const { return m_u_analytical_solution.data; }
    return_ptr GetV() const { return m_v_analytical_solution.data; }
    return_ptr GetW() const { return m_w_analytical_solution.data; }
    return_ptr GetP() const { return m_p_analytical_solution.data; }
    return_ptr GetT() const { return m_T_analytical_solution.data; }

 private:
    void gauss_bubble(real t);
    void exp_sinus_prod(real t);
    void exp_sinus_sum(real t);
    void hat(real t);
    void sin_sin_sin(real t);
    void mcDermott(real t);
    void vortex(real t);
    void vortex_y(real t);
    void beltrami(real t);
    void buoyancy_mms(real t);
    void zero(real t);

    Field m_u_analytical_solution;
    Field m_v_analytical_solution;
    Field m_w_analytical_solution;
    Field m_p_analytical_solution;
    Field m_T_analytical_solution;
    void (Solution::*m_init_function)(const real);

    real m_current_time_step = -1;
    bool m_has_analytical_solution;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_ANALYSIS_SOLUTION_H_ */
