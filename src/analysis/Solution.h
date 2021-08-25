/// \file       Solution.h
/// \brief      Calculates analytical solution
/// \details    This class calculates the analytical solution of different test cases
/// \date       Jul 11, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ANALYSIS_SOLUTION_H_
#define ARTSS_ANALYSIS_SOLUTION_H_

#include <string>
#include "../Domain.h"
#include "../Functions.h"
#include "../field/Field.h"
#include "../utility/Utility.h"
#include "../utility/Parameters.h"
#include "../utility/GlobalMacrosTypes.h"

class Solution {
 public:
    Solution(const Domain &domain, std::string initial_condition);

    void calc_analytical_solution(real t);

    // Getter
    real* GetU_data() const { return u_a.data; }
    real* GetV_data() const { return v_a.data; }
    real* GetW_data() const { return w_a.data; }
    real* GetP_data() const { return p_a.data; }
    real* GetT_data() const { return T_a.data; }
    return_ptr GetU() const { return u_a.data; }
    return_ptr GetV() const { return v_a.data; }
    return_ptr GetW() const { return w_a.data; }
    return_ptr GetP() const { return p_a.data; }
    return_ptr GetT() const { return T_a.data; }

 private:
    const Domain &m_domain;
    void init(std::string initial_condition);

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


    Field u_a, v_a, w_a;
    Field p_a;
    Field T_a;
    void (Solution::*m_init_function)(const real);

    real m_current_time_step = -1;
    bool m_has_analytical_solution;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_ANALYSIS_SOLUTION_H_ */
