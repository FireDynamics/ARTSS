/// \file       Layers.h
/// \brief      Adaption class for initial condition with layers (layersT)
/// \date       Dec 04, 2018
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ADAPTION_LAYERS_H_
#define ARTSS_ADAPTION_LAYERS_H_

#include "../field/Field.h"
#include "../interfaces/IAdaptionFunction.h"
#include "../solver/SolverController.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/settings/Settings.h"

class Layers : public IAdaptionFunction {
public:
    Layers(const Settings::adaption_classes::layers &settings, FieldController *field_controller);

    bool update(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2) override;
    void apply_changes(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2) override;
    bool has_reduction() override { return false; }

private:
    const Settings::adaption_classes::layers &m_settings;
    void adaptXDirection(real temperature, size_t no_buffer_cell, long *p_shift_x1, long *p_shift_x2);
    void adaptXDirection_serial(real temperature, size_t no_buffer_cell, long *p_shift_x1, long *p_shift_x2);

    void setXValues(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2, bool start);

    size_t getExpansionSize();

    size_t m_minimal;
    size_t m_timecounter;

    Field &m_T, &m_Ta, &m_Nu, &m_kappa, &m_gamma;
    real m_x1, m_x2, m_y1, m_y2, m_z1, m_z2;
    size_t m_nx, m_ny, m_nz;
};

#endif /* ARTSS_ADAPTION_LAYERS_H_ */
