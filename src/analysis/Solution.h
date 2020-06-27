/// \file 		Solution.h
/// \brief 		Calculates analytical solution
/// \details	This class calculates the analytical solution of different test cases
/// \date 		Jul 11, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ANALYSIS_SOLUTION_H_
#define ARTSS_ANALYSIS_SOLUTION_H_

#include "../Field.h"
#include "../utility/GlobalMacrosTypes.h"

class Solution {
public:
    Solution();

    virtual ~Solution();

    void CalcAnalyticalSolution(real t);

    // Getter
    return_ptr GetU() const {
        return ua->data;
    }

    return_ptr GetV() const {
        return va->data;
    }

    return_ptr GetW() const {
        return wa->data;
    }

    return_ptr GetP() const {
        return pa->data;
    }

    return_ptr GetT() const {
        return Ta->data;
    }

private:
    void init();

    void gauss_bubble(real t);
    void exp_sinus_prod(real t);
    void exp_sinus_sum(real t);
    void hat(real t);
    void fac_sin_sin_sin(real t);
    void mcDermott(real t);
    void vortex(real t);
    void vortex_y(real t);
    void beltrami(real t);
    void buoyancy_mms(real t);
    void zero(real t);


    Field *ua, *va, *wa;
    Field *pa;
    Field *Ta;
    real m_current_timestep = -1;
    void (Solution::*m_init_function)(const real);
};

#endif /* ARTSS_ANALYSIS_SOLUTION_H_ */
