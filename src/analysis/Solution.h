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
    Field *ua, *va, *wa;
    Field *pa;
    Field *Ta;

    void SetUp();
};

#endif /* ARTSS_ANALYSIS_SOLUTION_H_ */
