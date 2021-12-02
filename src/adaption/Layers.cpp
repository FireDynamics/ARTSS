/// \file       Layers.cpp
/// \brief      Adaption class for initial condition with layers (layersT)
/// \date       Dec 04, 2018
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Layers.h"
#include <chrono>
#include "../DomainData.h"
#include "Adaption.h"

Layers::Layers(Settings::Settings const &settings, FieldController *field_controller) :
        m_T(field_controller->get_field_T()),
        m_Ta(field_controller->get_field_T_ambient()),
        m_Nu(field_controller->get_field_nu_t()),
        m_kappa(field_controller->get_field_kappa()),
        m_gamma(field_controller->get_field_gamma()) {
    auto domain_data = DomainData::getInstance();
    m_minimal = static_cast<size_t> (std::pow(2, domain_data->get_levels()));
    m_time_counter = 0;

    m_no_buffer_cells = settings.get_size_t("adaption/class/buffer");
    m_check_value = settings.get_real("adaption/class/check_value");
    m_time_step = settings.get_size_t("adaption/class/timestep");
    m_expansion_size = settings.get_size_t("adaption/class/expansion_size");
}


// ==================================== Update ====================================
// ********************************************************************************
/// \brief  Checks for adaption
/// \return  bool if adaption is possible true
// ********************************************************************************
bool Layers::update(Coordinate<long> *shift_start, Coordinate<long> *shift_end) {
    m_time_counter++;
    if (m_time_counter < m_time_step) {
        return false;
    } else {
        m_time_counter = 0;
    }
    auto domain_data = DomainData::getInstance();

    shift_start->set_coordinate(0, 0, 0);
    shift_end->set_coordinate(0, 0, 0);

    m_x1 = domain_data->get_x1();
    m_x2 = domain_data->get_x2();
    m_nx = domain_data->get_nx();

    m_y1 = domain_data->get_y1();
    m_y2 = domain_data->get_y2();
    m_ny = domain_data->get_ny();

    m_z1 = domain_data->get_z1();
    m_z2 = domain_data->get_z2();
    m_nz = domain_data->get_nz();

    adapt(m_check_value, m_no_buffer_cells, shift_start, shift_end, X);
    adapt(m_check_value, m_no_buffer_cells, shift_start, shift_end, Z);

    size_t sum_adaption = 0;
    for (size_t a = 0; a < number_of_axis; a++) {
        CoordinateAxis axis = CoordinateAxis(a);
        size_t adaption = (*shift_start)[axis] + (*shift_end)[axis];
        sum_adaption += adaption;
        if (adaption > 0) {
            (*shift_start)[axis] = (-m_minimal * get_expansion_size());
            (*shift_end)[axis] = (-m_minimal * get_expansion_size());

            // boundary check
            long len_x1 = domain_data->get_start_index_CD(axis) - 1;
            long len_x2 = domain_data->get_end_index_CD(axis) + 1;
            if (len_x1 < -(*shift_start)[axis]) {
                (*shift_end)[axis] = ((*shift_end)[axis] - (*shift_start)[axis] - len_x1);
                (*shift_start)[axis] = (-len_x1);
                if (len_x2 < (*shift_end)[axis]) {
                    (*shift_end)[axis] = len_x2;
                }
            }
            if (len_x2 < (*shift_end)[axis]) {
                (*shift_start)[axis] = (*shift_start)[axis] - ((*shift_end)[axis] - len_x2);
                (*shift_end)[axis] = len_x2;
                if (len_x1 < -(*shift_start)[axis]) {
                    (*shift_start)[axis] = (-len_x1);
                }
            }
        }
    }
    return sum_adaption;
}

// ==================================== Get expansion size ====================================
// ********************************************************************************
/// \brief  In case of dynamic expansion size, the calculation should be done here
/// \return  site_t (calculated) expansion size
// ********************************************************************************
size_t Layers::get_expansion_size() {
    return m_expansion_size;
}

// ==================================== Set x values ====================================
// ********************************************************************************
/// \brief  Set values for new domain_data in x-direction
/// \params start x-values at x1 (start = true) or x-values at x2 (start=false)
// ********************************************************************************
void Layers::set_values(Coordinate<long> *shift_start, Coordinate<long> *shift_end, bool start, CoordinateAxis axis) {
    DomainData *domain_data = DomainData::getInstance();
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

    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    size_t j_start = domain_data->get_start_index_CD(other_axes[0]);
    size_t j_end = domain_data->get_end_index_CD(other_axes[0]);
    size_t k_start = domain_data->get_start_index_CD(other_axes[1]);
    size_t k_end = domain_data->get_start_index_CD(other_axes[1]);

    real *data_nu = m_Nu.data;
    real *data_gamma = m_gamma.data;
    real *data_kappa = m_kappa.data;
    real *data_temp = m_T.data;
    real *data_tempA = m_Ta.data;

    auto tmp = new Coordinate<size_t>();
    if (start) {
        size_t i_begin = domain_data->get_start_index_CD(axis);
        long shift = (*shift_start)[axis];
        size_t index;
        size_t idx;
#pragma acc parallel loop collapse(3) present(m_Nu, m_gamma, m_kappa, m_T, m_Ta)
        for (size_t j = j_start; j <= j_end; j++) {
            for (size_t k = k_start; k <= k_end; k++) {
                (*tmp)[other_axes[0]] = j;
                (*tmp)[other_axes[1]] = k;
                (*tmp)[axis] = i_begin + 1;
                index = tmp->get_index(Nx, Ny);
                for (int i = 0; i >= shift; i--) {
                    (*tmp)[axis] = i_begin + i;
                    idx = tmp->get_index(Nx, Ny);
                    *(data_nu + idx) = *(data_nu + index);
                    *(data_gamma + idx) = *(data_gamma + index);
                    *(data_kappa + idx) = *(data_kappa + index);
                    *(data_temp + idx) = *(data_temp + index);
                    *(data_tempA + idx) = *(data_tempA + index);
                }
            }
        }
    } else {
        size_t i_end = domain_data->get_end_index_CD(axis);
        long shift = (*shift_end)[axis];
        size_t index;
        size_t idx;
#pragma acc parallel loop collapse(3) present(m_Nu, m_gamma, m_kappa, m_T, m_Ta)
        for (size_t j = j_start; j < j_end; j++) {
            for (size_t k = k_start; k < k_end; k++) {
                (*tmp)[other_axes[0]] = j;
                (*tmp)[other_axes[1]] = k;
                (*tmp)[axis] = i_end;
                index = tmp->get_index(Nx, Ny);
                for (int i = 0; i <= shift; i++) {
                    (*tmp)[axis] = i_end + i + 1;
                    idx = tmp->get_index(Nx, Ny);
                    *(data_nu + idx) = *(data_nu + index);
                    *(data_gamma + idx) = *(data_gamma + index);
                    *(data_kappa + idx) = *(data_kappa + index);
                    *(data_temp + idx) = *(data_temp + index);
                    *(data_tempA + idx) = *(data_tempA + index);
                }
            }
        }
    }
    delete tmp;
    delete[] other_axes;
}

// ==================================== Apply changes =============================
// ********************************************************************************
/// \brief  Set values for new domain_data
// ********************************************************************************
void Layers::apply_changes(Coordinate<long> *shift_start, Coordinate<long> *shift_end) {
    for (size_t axis = 0; axis < number_of_axis; axis++) {
        if ((*shift_start)[axis] != 0) {
            Layers::set_values(shift_start, shift_end, true, CoordinateAxis(axis));
        }
        if ((*shift_end)[axis] != 0) {
            Layers::set_values(shift_start, shift_end, false, CoordinateAxis(axis));
        }
    }
}

// ==================================== Adaption x direction parallel ===============================
// ***************************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param  checkValue check value
/// \param  no_buffer_cell Buffersize
// ***************************************************************************************
void Layers::adapt(real checkValue, size_t no_buffer_cell, Coordinate<long> *start, Coordinate<long> *end, CoordinateAxis axis) {
    auto domain_data = DomainData::getInstance();
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

    auto data = m_T.data;
    size_t size __attribute__((unused)) = domain_data->get_size();

    size_t expansion_counter_start = 0;
    size_t expansion_counter_end = 0;

    bool expansion_start = (domain_data->get_x1() != domain_data->get_X1());
    bool expansion_end = (domain_data->get_x2() != domain_data->get_X2());
    auto tmp = new Coordinate<size_t>();

#pragma acc data present(data[:size]) copyout(expansion_counter_end, expansion_counter_start) copyin(tmp)
    {
        size_t Nx = domain_data->get_Nx();
        size_t Ny = domain_data->get_Ny();

        //loop through left and right side of cuboid in x direction
#pragma acc parallel loop collapse(2) present(data[:size]) reduction(+:expansion_counter_end, expansion_counter_start)
        for (size_t j = (*start)[other_axes[0]]; j <= (*end)[other_axes[0]]; j++) {
            for (size_t k = (*start)[other_axes[1]]; k <= (*end)[other_axes[1]]; k++) {
                (*tmp)[other_axes[1]] = k;
                (*tmp)[other_axes[0]] = j;
                (*tmp)[axis] = (*end)[axis] - no_buffer_cell + 1;
                size_t index = tmp->get_index(Nx, Ny);
                // check innermost plane of the buffer zone on the right side
                if (*(data + index) > checkValue) {
                    expansion_counter_end++;
                }

                (*tmp)[axis] = (*start)[axis] + no_buffer_cell - 1;
                index = tmp->get_index(Nx, Ny);
                if (*(data + index) > checkValue) {
                    expansion_counter_start++;
                }
            }
        }
    }
#pragma acc exit data delete(tmp)
    if (expansion_counter_end > 0 && expansion_end) {
        (*end)[axis] = 1;
    }
    if (expansion_counter_start > 0 && expansion_start) {
        (*start)[axis] = 1;
    }
    delete tmp;
    delete[] other_axes;
}
