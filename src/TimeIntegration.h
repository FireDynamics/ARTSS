/// \file 		TimeIntegration.h
/// \brief 		Runs the time loop
/// \date 		May 20, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_TIMEINTEGRATION_H_
#define ARTSS_TIMEINTEGRATION_H_

#include "interfaces/ISolver.h"
#include "interfaces/ISource.h"
#include "utility/GlobalMacrosTypes.h"
#include "analysis/Analysis.h"
#include "utility/Utility.h"
#ifndef PROFILING
#include <spdlog/logger.h>
#endif

class TimeIntegration {
public:
	TimeIntegration(ISolver *isolv, const char *fname);

	void run();

private:
	ISolver* m_solver;
	const char *m_fname;
	real m_dt;
	real m_t_end;
	real m_t_cur;
	size_t m_size = 0;
#ifndef PROFILING
	std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_TIMEINTEGRATION_H_ */
