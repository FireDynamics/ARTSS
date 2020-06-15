/// \file 		PressureI.h
/// \brief 		Interface for pressure method
/// \date 		Sep 14, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACES_PRESSUREI_H_
#define ARTSS_INTERFACES_PRESSUREI_H_

#include "../Field.h"

class IPressure {
public:
	IPressure() = default;
	virtual ~IPressure() = default;

	virtual void pressure(Field* out, Field* b, real t, bool sync)=0;

	void Divergence(Field* out, const Field* inx, const Field* iny, const Field* inz, bool sync);
	void Project(Field* outu, Field* outv, Field* outw, const Field* inu, const Field* inv, const Field* inw, const Field* inp, bool sync);
};

#endif /* ARTSS_INTERFACES_PRESSUREI_H_ */
