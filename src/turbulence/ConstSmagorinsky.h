/// \file         ConstSmagorinsky.h
/// \brief        calculates eddy viscosity based on Constant Smagorinsky-Lilly LES model
/// \date         Aug 18, 2016
/// \author       Suryanarayana Maddu
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_TURBULENCE_CONSTSMAGORINSKY_H_
#define ARTSS_TURBULENCE_CONSTSMAGORINSKY_H_

#include "../interfaces/ITurbulence.h"
#include "../utility/Settings.h"

class ConstSmagorinsky : public ITurbulence {
 public:
    explicit ConstSmagorinsky(Settings const &settings);
    ~ConstSmagorinsky() override = default;

    void calc_turbulent_viscosity(
            Field &ev,
            const Field &in_u, const Field &in_v, const Field &in_w,
            bool sync) override;

    void explicit_filtering(Field &out, const Field &in, bool sync) override;

private:
    real m_Cs;
};

#endif /* ARTSS_TURBULENCE_CONSTSMAGORINSKY_H_ */
