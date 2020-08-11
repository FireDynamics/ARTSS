/// \file       Solution.h
/// \brief      Calculates analytical solution
/// \details    This class calculates the analytical solution of different test cases
/// \date       Jul 11, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ANALYSIS_SOLUTION_H_
#define ARTSS_ANALYSIS_SOLUTION_H_

#include "../Field.h"
#include "../utility/Utility.h"
#include "../utility/GlobalMacrosTypes.h"

class Solution {
public:
    Solution();

    virtual ~Solution();

    void calc_analytical_solution(real t);

    // Getter
    return_ptr GetU() const { return u_a->data; }
    return_ptr GetV() const { return v_a->data; }
    return_ptr GetW() const { return w_a->data; }
    return_ptr GetP() const { return p_a->data; }
    return_ptr GetT() const { return T_a->data; }

private:
    void init();

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


    Field *u_a, *v_a, *w_a;
    Field *p_a;
    Field *T_a;
    void (Solution::*m_init_function)(const real);

    real m_current_time_step = -1;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_ANALYSIS_SOLUTION_H_ */
