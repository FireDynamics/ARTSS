/// \file 		VCycleMG.h
/// \brief 		Defines V-cycle of geometric multigrid method
/// \date 		Sep 14, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_PRESSURE_VCYCLEMG_H_
#define ARTSS_PRESSURE_VCYCLEMG_H_

#include <vector>
#include "../interfaces/IPressure.h"
#include "../Field.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"



class VCycleMG: public IPressure{
public:
	VCycleMG(Field* out, Field* b);
	~VCycleMG() override;

	void pressure(Field* out, Field* b, real t, bool sync);

private:
	void VCycleMultigrid(Field* out, bool sync = true);
	void UpdateInput(Field* out, Field* b, bool sync = true);
	void Smooth(Field* out, Field* tmp, Field* b, size_t level, bool sync = true);
	void Residuum(Field* out, Field* in, Field* b, size_t level, bool sync = true);
	void Restrict(Field* out, Field* in, size_t level, bool sync = true);
	void Prolongate(Field* out, Field* in, size_t level, bool sync = true);
	void Solve( Field* out, Field* tmp, Field* b, size_t level, bool sync = true);

	size_t levels;
	size_t cycles;
	size_t relaxs;

	real m_dsign;
	real m_w;

	std::vector<Field*> residuum0;
	std::vector<Field*> residuum1;
	std::vector<Field*> err0;
	std::vector<Field*> error1;
	std::vector<Field*> mg_temporal_solution;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_PRESSURE_VCYCLEMG_H_ */
