/// \file       Vortex.h
/// \brief      Adaption class for initial condition with vortex
/// \date       Dec 04, 2018
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_ADAPTION_VORTEX_H_
#define ARTSS_ADAPTION_VORTEX_H_

#include "../interfaces/IAdaptionFunction.h"
#include "Adaption.h"

class Vortex : public IAdaptionFunction {
public:
    explicit Vortex(ISolver *solver);

    bool update(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2) override;
    void apply_changes(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2) override;
    bool has_reduction() override;

private:
    void Drift_dynamic(const size_t *arr_idx, size_t arr_idx_size);
    void Zero(size_t *arr_idx, size_t arr_idx_size);

    real m_u_lin;
    real m_v_lin;
    real m_w_lin;
    size_t m_minimal;
    Field *u, *v, *w;
    bool m_reduction;
    real m_threshold;
    size_t m_buffer;
    bool m_x_side = false;
    bool m_y_side = false;
    bool m_z_side = false;
};

#endif /* ARTSS_ADAPTION_VORTEX_H_ */
