/// \file 		SLAdvect.h
/// \brief 		Solves advection equation via unconditionally stable Semi-Langrangian approach
/// \date 		Aug 23, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ADVECTION_SLADVECT_H_
#define ARTSS_ADVECTION_SLADVECT_H_

#include "../Interfaces/AdvectionI.h"
#include "../Field.h"
#include "../Utility/GlobalMacrosTypes.h"

class SLAdvect: public AdvectionI {
public:
	SLAdvect();
    ~SLAdvect() override = default;

	void advect(Field* out, Field* in, const Field* u_vel, const Field* v_vel, const Field* w_vel, bool sync) override;

private:
	real m_dt;
};

#endif /* ARTSS_ADVECTION_SLADVECT_H_ */
