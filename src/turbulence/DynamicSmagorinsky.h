/// \file         DynamicSmagorinsky.cpp
/// \brief        calculates eddy viscosity based on Dynamic Smagorinsky-Lilly LES model
/// \date         September 8, 2016
/// \author       Suryanarayana Maddu
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_TURBULENCE_DYNAMICSMAGORINSKY_H_
#define ARTSS_TURBULENCE_DYNAMICSMAGORINSKY_H_

#include "../interfaces/ITurbulence.h"
#include "../Field.h"

class DynamicSmagorinsky : public ITurbulence {
public:
    DynamicSmagorinsky();
    ~DynamicSmagorinsky() override;

    void CalcTurbViscosity(Field *ev, Field *in_u, Field *in_v, Field *in_w, bool sync) override;
    void ExplicitFiltering(Field *out, const Field *in, bool sync) override;

private:
    Field *u_f, *v_f, *w_f;                 // filtered velocities
    Field *uu, *vv, *ww, *uv, *uw, *vw;           // velocity products
    Field *uu_f, *vv_f, *ww_f, *uv_f, *uw_f, *vw_f;     // filters of the velocity products
    Field *uf_uf, *vf_vf, *wf_wf, *uf_vf, *uf_wf, *vf_wf; // product of the filtered velocities
    Field *L11, *L22, *L33, *L12, *L13, *L23;         // Leonard stress
    Field *S11, *S22, *S33, *S12, *S13, *S23;         // strain tensor
    Field *S11_f, *S22_f, *S33_f, *S12_f, *S13_f, *S23_f;   // second filtered strain tensor
    Field *P11, *P22, *P33, *P12, *P13, *P23;           // Product of strain modulus and strain tensor
    Field *P11_f, *P22_f, *P33_f, *P12_f, *P13_f, *P23_f;   // second filter for the above
    Field *M11, *M22, *M33, *M12, *M13, *M23;       // High frequency resolved terms
    Field *S_bar, *S_bar_f;                 // modulus of strain tensor
    Field *Cs;                        // dynamic constant
    real m_nu;                        // viscosity
};

#endif /* ARTSS_TURBULENCE_DYNAMICSMAGORINSKY_H_ */
