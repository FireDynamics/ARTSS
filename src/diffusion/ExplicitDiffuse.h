/// \file       ExplicitDiffuse.h
/// \brief      Solves diffusion equation with an explicit method
/// \date       December 12, 2019
/// \author     Max Boehler
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DIFFUSION_EXPLICITDIFFUSE_H_
#define ARTSS_DIFFUSION_EXPLICITDIFFUSE_H_

#include "../interfaces/IDiffusion.h"
#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"


class ExplicitDiffuse : public IDiffusion {
public:
    explicit ExplicitDiffuse(Settings::Settings const &settings) : m_settings(settings) {}

    void diffuse(
            Field &out, const Field &in, const Field &b,
            real D, bool sync) override;
    void diffuse(
            Field &out, const Field &in, const Field &b,
            real D, const Field &EV, bool sync) override;  // turbulent version

    void ExplicitStep(Field &out, const Field &in,
                      real D, bool sync = true);
    void ExplicitStep(Field &out, const Field &in,
                      real D, const Field &EV, bool sync = true);

private:
    Settings::Settings const &m_settings;
};

#endif /* ARTSS_DIFFUSION_EXPLICITDIFFUSE_H_ */
