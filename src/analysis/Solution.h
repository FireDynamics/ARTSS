/// \file 		Solution.h
/// \brief 		Calculates analytical solution
/// \details	This class calculates the analytical solution of different test cases
/// \date 		Jul 11, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ANALYSIS_SOLUTION_H_
#define ARTSS_ANALYSIS_SOLUTION_H_

#ifndef PROFILING
#include <spdlog/logger.h>
#endif
#include "../Field.h"
#include "../utility/GlobalMacrosTypes.h"

class Solution {
public:
	Solution();
	virtual ~Solution();

	void CalcAnalyticalSolution(const real t);

	// Getter
	return_ptr GetU() const {
		return ua->data;
	}
	return_ptr GetU0() const {
		return u0a->data;
	}
	return_ptr GetV() const {
		return va->data;
	}
	return_ptr GetV0() const {
		return v0a->data;
	}
	return_ptr GetW() const {
		return wa->data;
	}
	return_ptr GetW0() const {
		return w0a->data;
	}
	return_ptr GetP() const {
		return pa->data;
	}
	return_ptr GetP0() const {
		return p0a->data;
	}
	return_ptr GetT() const {
		return Ta->data;
	}
	return_ptr GetT0() const {
		return T0a->data;
	}

private:
	Field* ua, *va, *wa;
	Field* u0a, *v0a, *w0a;
	Field* pa, *p0a;
	Field* Ta, *T0a;

    std::shared_ptr<spdlog::logger> m_logger;
    void SetUp();
};

#endif /* ARTSS_ANALYSIS_SOLUTION_H_ */
