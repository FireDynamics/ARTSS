/// \file 		ExplicitAdvect.h
/// \brief 		Explicit (BD) solver for advection equation
/// \date 		Aug 22, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ADVECTION_EXPLICITADVECT_H_
#define ARTSS_ADVECTION_EXPLICITADVECT_H_

#include "../Interfaces/AdvectionI.h"
#include "../Field.h"
#include "../Utility/GlobalMacrosTypes.h"

class ExplicitAdvect: public AdvectionI {
public:
	ExplicitAdvect();
    ~ExplicitAdvect() override = default;
	void advect(Field* out, Field* in, const Field* u_vel, const Field* v_vel, const Field* w_vel, bool sync) override;

private:
	real m_dt;
};

#endif /* ARTSS_ADVECTION_EXPLICITADVECT_H_ */
