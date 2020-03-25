/// \file 		SourceI.h
/// \brief 		Interface for adding sources
/// \date 		Dec 2, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_INTERFACES_SOURCEI_H_
#define ARTSS_INTERFACES_SOURCEI_H_

#include "../Field.h"
#include "../utility/GlobalMacrosTypes.h"

class SourceI {
public:
	SourceI() = default;
	virtual ~SourceI() = default;

	virtual void addSource(Field* outx, Field* outy, Field* outz, Field* Sx, Field* Sy, Field* Sz, bool sync)=0;
	virtual void addSource(Field* out, Field* S, bool sync)=0;

	void BuoyancyForce(Field* out, const Field* in, const Field* ina, bool sync = true);
	void BuoyancyST_MMS(Field* out, real t, bool sync = true);
	void Gauss(Field* out, real HRR, real cp, real x0, real y0, real z0, real sigmax, real sigmay, real sigmaz, bool sync = true);
	void Dissipate(Field* out, const Field* inu, const Field* inv, const Field* inw, bool sync = true);
};

#endif /* ARTSS_INTERFACES_SOURCEI_H_ */
