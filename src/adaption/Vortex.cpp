/// \file       Vortex.cpp
/// \brief      Adaption class for initial condition with vortex
/// \date       Dec 04, 2018
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Vortex.h"
#include "../DomainData.h"

Vortex::Vortex(Settings::Settings const &settings, FieldController *field_controller) :
        m_settings(settings),
        m_u(field_controller->get_field_u()),
        m_v(field_controller->get_field_v()),
        m_w(field_controller->get_field_w()) {
    auto domain_data = DomainData::getInstance();
    m_u_lin = m_settings.get_real("initial_conditions/u_lin");
    m_v_lin = m_settings.get_real("initial_conditions/v_lin");
    m_w_lin = m_settings.get_real("initial_conditions/w_lin");
    m_minimal = static_cast<size_t> (std::pow(2, domain_data->get_levels()));
    m_reduction = m_settings.get_bool("adaption/class/reduction/enabled");
    if (m_reduction) {
        std::string dir = (m_settings.get("adaption/class/reduction/dir"));
        if (dir.find('x') != std::string::npos) m_x_side = true;
        if (dir.find('y') != std::string::npos) m_y_side = true;
        if (dir.find('z') != std::string::npos) m_z_side = true;
    }
    m_buffer = m_settings.get_int("adaption/class/buffer");
    m_threshold = m_u_lin * m_settings.get_real("adaption/class/threshold");
}

// ==================================== Has reduction ===============================
// ***************************************************************************************
/// \brief  Checks if reduction is enabled
/// \return bool true if yes false if no
// ***************************************************************************************
bool Vortex::has_reduction() {
    return m_reduction;
}

// ==================================== Update ====================================
// ********************************************************************************
/// \brief  Checks for adaption
/// \return  bool if adaption is possible true
// ********************************************************************************
bool Vortex::update(Coordinate<long> *shift_start, Coordinate<long> *shift_end) {
    bool adaption = false;

    shift_start->set_coordinate(0, 0, 0);
    shift_end->set_coordinate(0, 0, 0);

    adaption = Adaption::adapt(m_settings, m_u, m_u_lin, m_buffer, m_threshold, shift_start, shift_end, m_minimal, m_reduction, X) || adaption;
    if (m_y_side) {
        adaption = Adaption::adapt(m_settings, m_v, m_v_lin, m_buffer, m_threshold, shift_start, shift_end, m_minimal, m_reduction, Y) || adaption;
    }

    *shift_start *= m_minimal;
    *shift_end *= m_minimal;
    return adaption;
}

// ==================================== Apply changes =============================
// ********************************************************************************
/// \brief  Set values for new domain_data
// ********************************************************************************
void Vortex::apply_changes(Coordinate<long> *shift_start, Coordinate<long> *shift_end) {
    auto domain_data = DomainData::getInstance();
    for (size_t axis = 0; axis < number_of_axis; axis++) {
        auto other_axes = new CoordinateAxis[2];
        if (axis == CoordinateAxis::X) {
            other_axes[0] = CoordinateAxis::Y;
            other_axes[1] = CoordinateAxis::Z;
        } else if (axis == CoordinateAxis::Y) {
            other_axes[0] = CoordinateAxis::X;
            other_axes[1] = CoordinateAxis::Z;
        } else if (axis == CoordinateAxis::Z) {
            other_axes[0] = CoordinateAxis::X;
            other_axes[1] = CoordinateAxis::Y;
        }
        for (auto shift: {shift_start, shift_end}) {
            if ((*shift)[axis] != 0) {
                auto coordinate_axis = CoordinateAxis(axis);
                if ((*shift)[axis] < 0) {
                    size_t len_e =
                            (domain_data->get_end_index_CD(other_axes[0]) - domain_data->get_start_index_CD(other_axes[0]))
                          * (domain_data->get_end_index_CD(other_axes[1]) - domain_data->get_start_index_CD(other_axes[1]))
                          * static_cast<size_t>(fabs((*shift)[axis]));
                    auto *arr_idxExpansion = new size_t[len_e];
#pragma acc enter data create(arr_idxExpansion[:len_e])
                    Adaption::expand(shift, true, arr_idxExpansion, len_e, coordinate_axis);
                    Vortex::Drift_dynamic(arr_idxExpansion, len_e);
#pragma acc exit data delete(arr_idxExpansion[:len_e])
                    delete[] arr_idxExpansion;
                } else {
                    size_t len_r =
                            (domain_data->get_end_index_CD(other_axes[0]) - domain_data->get_start_index_CD(other_axes[0]))
                          * (domain_data->get_end_index_CD(other_axes[1]) - domain_data->get_start_index_CD(other_axes[1]))
                          * static_cast<size_t>(fabs((*shift)[axis]));
                    auto *arr_idxReduction = new size_t[len_r];
#pragma acc enter data create(arr_idxReduction[:len_r])
                    Adaption::reduce(shift, true, arr_idxReduction, len_r, coordinate_axis);
                    Zero(arr_idxReduction, len_r);
#pragma acc exit data delete(arr_idxReduction[:len_r])
                    delete[] arr_idxReduction;
                }
            }
        }
        delete[] other_axes;
    }
}

// ==================================== Zero =============================
// ********************************************************************************
/// \brief  Set values to zero (parallel version of function zero)
/// \param  arr_idx Index list of cell which should be set to zero
/// \param  arr_idx_size Size of arr_idx
// ********************************************************************************
void Vortex::Zero(size_t *arr_idx, size_t arr_idx_size) {
    auto data_u = m_u.data;
    //  auto data_v = v->data;
    //  auto data_w = w->data;
#pragma acc parallel loop independent present(data_u[:arr_idx_size], arr_idx[:arr_idx_size])
    for (size_t idx = 0; idx < arr_idx_size; idx++) {
        *(data_u + arr_idx[idx]) = 0;
        //      *(data_v+arr_idx[idx])=0;
        //      *(data_w+arr_idx[idx])=0;
    }
}

// ==================================== Drift dynamic =============================
// ********************************************************************************
/// \brief  Set values to background flow (parallel version of function drift)
/// \param  arr_idx Index list of cell which should be set
/// \param  arr_idx_size Size of arr_idx
// ********************************************************************************
void Vortex::Drift_dynamic(const size_t *arr_idx, size_t arr_idx_size) {
    auto data_x = m_u.data;
    auto data_y = m_v.data;
    auto data_z = m_w.data;
    // inner cells

#pragma acc parallel loop independent present(data_x[:arr_idx_size], data_z[:arr_idx_size], data_y[:arr_idx_size], arr_idx[:arr_idx_size])
    for (size_t idx = 0; idx < arr_idx_size; idx++) {
        *(data_x + *(arr_idx + idx)) = m_u_lin;
        *(data_y + *(arr_idx + idx)) = m_v_lin;
        *(data_z + *(arr_idx + idx)) = m_w_lin;
    }
}
