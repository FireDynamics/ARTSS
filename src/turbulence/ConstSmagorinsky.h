/// \file         ConstSmagorinsky.h
/// \brief        calculates eddy viscosity based on Constant Smagorinsky-Lilly LES model
/// \date         Aug 18, 2016
/// \author       Suryanarayana Maddu
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_TURBULENCE_CONSTSMAGORINSKY_H_
#define ARTSS_TURBULENCE_CONSTSMAGORINSKY_H_

#include "../interfaces/ITurbulence.h"

class ConstSmagorinsky : public ITurbulence {
public:
    ConstSmagorinsky();
    ~ConstSmagorinsky() override = default;

    void CalcTurbViscosity(Field *ev, Field *in_u, Field *in_v, Field *in_w, bool sync) override;

    void ExplicitFiltering(Field *out, const Field *in, bool sync) override;

private:
    real m_nu;
    real m_dt;
    real m_Cs;
};

#endif /* ARTSS_TURBULENCE_CONSTSMAGORINSKY_H_ */
