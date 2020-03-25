/// \file 		ExplicitSource.h
/// \brief 		Adding source via Explicit Euler
/// \date 		Dec 2, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_SOURCE_EXPLICITEULERSOURCE_H_
#define ARTSS_SOURCE_EXPLICITEULERSOURCE_H_

#include "../Interfaces/SourceI.h"
#include "../Field.h"
#include "../Utility/GlobalMacrosTypes.h"

class ExplicitEulerSource: public SourceI {
public:
	ExplicitEulerSource();

	void addSource(Field* outx, Field* outy, Field* outz, Field* Sx, Field* Sy, Field* Sz, bool sync) override;
	void addSource(Field* out, Field* S, bool sync) override;

private:
	real m_dt;
	std::string m_dir_vel ="";
};

#endif /* ARTSS_SOURCE_EXPLICITEULERSOURCE_H_ */
