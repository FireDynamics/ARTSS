/// \file       Layers.cpp
/// \brief      Adaption class for initial condition with layers (layersT)
/// \date       Dec 04, 2018
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Layers.h"
#include <chrono>
#include "../Domain.h"
#include "../utility/Parameters.h"
#include "Adaption.h"

Layers::Layers(ISolver *solver) {
    auto params = Parameters::getInstance();
    auto domain = Domain::getInstance();
    m_minimal = static_cast<size_t> (std::pow(2, domain->get_levels()));
    m_timecounter = 0;

    m_no_buffer_cells = static_cast<size_t> (params->get_int("adaption/class/buffer"));
    m_check_value = params->get_real("adaption/class/check_value");
    m_timestep = static_cast<size_t> (params->get_int("adaption/class/timestep"));
    m_expansion_size = static_cast<size_t> (params->get_int("adaption/class/expansion_size"));

    m_T = solver->T;
    m_Ta = solver->T_ambient;

    m_Nu = solver->nu_t;

    m_kappa = solver->kappa_t;
    m_gamma = solver->gamma_t;

    //if (params->get("adaption/version") == "CPU"){
    //    m_fctP_adapt = adapt_x_direction_serial;
    //}else{
    //  m_fctP_adapt = adapt_x_direction;
    //}
}


// ==================================== Update ====================================
// ********************************************************************************
/// \brief  Checks for adaption
/// \return  bool if adaption is possible true
// ********************************************************************************
bool Layers::update(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2) {
    m_timecounter++;
    if (m_timecounter < m_timestep) {
        return false;
    } else {
        m_timecounter = 0;
    }
    auto domain = Domain::getInstance();

    *p_shift_x1 = 0;
    *p_shift_x2 = 0;
    *p_shift_y1 = 0;
    *p_shift_y2 = 0;
    *p_shift_z1 = 0;
    *p_shift_z2 = 0;

    m_x1 = domain->get_x1();
    m_x2 = domain->get_x2();
    m_nx = domain->get_nx();

    m_y1 = domain->get_y1();
    m_y2 = domain->get_y2();
    m_ny = domain->get_ny();

    m_z1 = domain->get_z1();
    m_z2 = domain->get_z2();
    m_nz = domain->get_nz();

    adaptXDirection(m_check_value, m_no_buffer_cells, p_shift_x1, p_shift_x2);
    //TODO z-direction

    size_t adaption = *p_shift_x1 + *p_shift_x2 + *p_shift_z1 + *p_shift_z2;
    if (adaption > 0) {
        //TODO z-direction
        *p_shift_x1 = (-m_minimal * getExpansionSize());
        *p_shift_x2 = (m_minimal * getExpansionSize());

        // boundary check
        long len_x1 = domain->get_index_x1() - 1;
        long len_x2 = static_cast<long> (std::round((domain->get_X2() - domain->get_x2()) / domain->get_dx()));
        if (len_x1 < -*p_shift_x1 || len_x2 < *p_shift_x2) {
            if (len_x1 < -*p_shift_x1) {
                *p_shift_x2 = (*p_shift_x2 - *p_shift_x1 - len_x1);
                *p_shift_x1 = (-len_x1);
                if (len_x2 < *p_shift_x2) {
                    *p_shift_x2 = len_x2;
                }
            }
            if (len_x2 < *p_shift_x2) {
                *p_shift_x1 = *p_shift_x1 - (*p_shift_x2 - len_x2);
                *p_shift_x2 = len_x2;
                if (len_x1 < -*p_shift_x1) {
                    *p_shift_x1 = (-len_x1);
                }
            }
        }

    }
    return adaption;
}

// ==================================== Get expansion size ====================================
// ********************************************************************************
/// \brief  In case of dynamic expansion size, the calculation should be done here
/// \return  site_t (calculated) expansion size
// ********************************************************************************
size_t Layers::getExpansionSize() {
    return m_expansion_size;
}

// ==================================== Set x values ====================================
// ********************************************************************************
/// \brief  Set values for new domain in x-direction
/// \params start x-values at x1 (start = true) or x-values at x2 (start=false)
// ********************************************************************************
void Layers::setXValues(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2, bool start) {
    Domain *domain = Domain::getInstance();

    size_t Nx = domain->get_Nx();
    size_t Ny = domain->get_Ny();

    size_t j_start = static_cast<size_t> (std::round((m_y1 - domain->get_Y1()) / domain->get_dy()));
    size_t j_end = j_start + m_ny;
    size_t k_start = static_cast<size_t> (std::round((m_z1 - domain->get_Z1()) / domain->get_dz()));
    size_t k_end = k_start + m_nz;

    real *data_nu = m_Nu->data;
    real *data_gamma = m_gamma->data;
    real *data_kappa = m_kappa->data;
    real *data_temp = m_T->data;
    real *data_tempA = m_Ta->data;

    size_t size = domain->get_size();

    if (start) {
        size_t nx_begin = static_cast<size_t> (std::round((m_x1 - domain->get_X1()) / domain->get_dx()));
        long shift = *p_shift_x1;
        size_t index;
        size_t idx;
#pragma acc parallel loop collapse(3) present(data_nu[:size], data_gamma[:size], data_kappa[:size], data_temp[:size], data_tempA[:size])
        for (size_t j = j_start; j < j_end; j++) {
            for (size_t k = k_start; k < k_end; k++) {
                for (int i = 0; i >= shift; i--) {
                    index = IX(nx_begin + 1, j, k, Nx, Ny);
                    idx = IX(nx_begin + i, j, k, Nx, Ny);
                    *(data_nu + idx) = *(data_nu + index);
                    *(data_gamma + idx) = *(data_gamma + index);
                    *(data_kappa + idx) = *(data_kappa + index);
                    *(data_temp + idx) = *(data_temp + index);
                    *(data_tempA + idx) = *(data_tempA + index);
                }
            }
        }
    } else {
        size_t nx_end = static_cast<size_t> (std::round((m_x2 - domain->get_X1()) / domain->get_dx()));
        long shift = *p_shift_x2;
        size_t index;
        size_t idx;
#pragma acc parallel loop collapse(3) present(data_nu[:size], data_gamma[:size], data_kappa[:size], data_temp[:size], data_tempA[:size])
        for (size_t j = j_start; j < j_end; j++) {
            for (size_t k = k_start; k < k_end; k++) {
                for (int i = 0; i <= shift; i++) {
                    index = IX(nx_end, j, k, Nx, Ny);
                    idx = IX(nx_end + i + 1, j, k, Nx, Ny);
                    *(data_nu + idx) = *(data_nu + index);
                    *(data_gamma + idx) = *(data_gamma + index);
                    *(data_kappa + idx) = *(data_kappa + index);
                    *(data_temp + idx) = *(data_temp + index);
                    *(data_tempA + idx) = *(data_tempA + index);
                }
            }
        }
    }
}

// ==================================== Apply changes =============================
// ********************************************************************************
/// \brief  Set values for new domain
// ********************************************************************************
void Layers::apply_changes(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2) {
    if (*p_shift_x1) {
        Layers::setXValues(p_shift_x1, p_shift_x2, p_shift_y1, p_shift_y2, p_shift_z1, p_shift_z2, true);
    }
    if (*p_shift_x2) {
        Layers::setXValues(p_shift_x1, p_shift_x2, p_shift_y1, p_shift_y2, p_shift_z1, p_shift_z2, false);
    }
}

// ==================================== Adaption x direction parallel ===============================
// ***************************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param  checkValue check value
/// \param  no_buffer_cell Buffersize
// ***************************************************************************************
void Layers::adaptXDirection(real checkValue, size_t no_buffer_cell, long *p_shift_x1, long *p_shift_x2) {
    auto domain = Domain::getInstance();

    auto data = m_T->data;
    size_t size = domain->get_size();

    size_t expansion_counter_start = 0;
    size_t expansion_counter_end = 0;

    bool expansion_start = (domain->get_x1() != domain->get_X1());
    bool expansion_end = (domain->get_x2() != domain->get_X2());

    //copy(expansion_counter_start,expansion_counter_end)
#pragma acc data present(data[:size]) copyout(expansion_counter_end, expansion_counter_start)
    {
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        size_t nx_begin = domain->get_index_x1();//((domain->get_x1() - domain->get_X1()) / domain->get_dx());

        size_t j_start = domain->get_index_y1();//((domain->get_y1() - domain->get_Y1()) / domain->get_dy());
        size_t j_end = j_start + domain->get_ny() - 1;
        size_t k_start = domain->get_index_z1();//((domain->get_z1() - domain->get_Z1()) / domain->get_dz());
        size_t k_end = k_start + domain->get_nz() - 1;

        size_t i_start = nx_begin;
        size_t i_end = nx_begin + domain->get_nx() - 1;

        //loop through left and right side of cuboid in x direction
#pragma acc parallel loop collapse(2) present(data[:size]) reduction(+:expansion_counter_end)
        for (size_t j = j_start; j < j_end; j++) {
            for (size_t k = k_start; k < k_end; k++) {
                // check innermost plane of the buffer zone on the right side
                if ((*(data + IX(i_end - no_buffer_cell + 1, j, k, Nx, Ny)) > checkValue)) {
                    expansion_counter_end++;
                }
            }
        }
//#pragma acc parallel loop collapse(2) present(data[:size]) reduction(max:expansion_counter_start)
#pragma acc parallel loop collapse(2) present(data[:size]) reduction(+:expansion_counter_start)
        for (size_t j = j_start; j < j_end; j++) {
            for (size_t k = k_start; k < k_end; k++) {
                // check innermost plane of the buffer zone on the left side
                if ((*(data + IX(i_start + no_buffer_cell - 1, j, k, Nx, Ny)) > checkValue)) {
                    expansion_counter_start++;
                }
            }
        }
    }
    if (expansion_counter_end > 0 && expansion_end) {
        *p_shift_x2 = 1;
    }
    if (expansion_counter_start > 0 && expansion_start) {
        *p_shift_x1 = 1;
    }
}

// ==================================== Adaption x direction serial ===============================
// ***************************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param  checkValue check value
/// \param  no_buffer_cell Buffersize
// ***************************************************************************************
void Layers::adaptXDirection_serial(real checkValue, size_t no_buffer_cell, long *p_shift_x1, long *p_shift_x2) {
    auto domain = Domain::getInstance();
    size_t Nx = domain->get_Nx();
    size_t Ny = domain->get_Ny();

    size_t nx_begin = domain->get_index_x1();//((domain->get_x1() - domain->get_X1()) / domain->get_dx());

    size_t j_start = domain->get_index_y1();//((domain->get_y1() - domain->get_Y1()) / domain->get_dy());
    size_t j_end = j_start + domain->get_ny() - 1;
    size_t k_start = domain->get_index_z1();//((domain->get_z1() - domain->get_Z1()) / domain->get_dz());
    size_t k_end = k_start + domain->get_nz() - 1;

    size_t i_start = nx_begin;
    size_t i_end = nx_begin + domain->get_nx() - 1;

    //expansion - expand if there is at least one cell in the buffer area fulfills the condition
    ADTypes expansion_start = ADTypes::UNKNOWN;
    if (domain->get_x1() == domain->get_X1()) {
        expansion_start = ADTypes::NO;
    }
    ADTypes expansion_end = ADTypes::UNKNOWN;
    if (domain->get_x2() == domain->get_X2()) {
        expansion_end = ADTypes::NO;
    }
    ADTypes expansion = ADTypes::UNKNOWN;
    if (expansion_end == expansion_start) {
        expansion = expansion_end;
    }

    auto data = m_T->data;

    //loop through left and right side of cuboid in x direction
    for (size_t j = j_start; j < j_end && expansion == ADTypes::UNKNOWN; j++) {
        for (size_t k = k_start; k < k_end && expansion == ADTypes::UNKNOWN; k++) {
            // check innermost plane of the buffer zone on the left side
            //size_t idx_s1 = IX(i_start + no_buffer_cell - 1, j, k, Nx, Ny);
            if (expansion_start == ADTypes::UNKNOWN &&
                (*(data + IX(i_start + no_buffer_cell - 1, j, k, Nx, Ny)) > checkValue)) {
                expansion_start = ADTypes::YES;
            }
            // check innermost plane of the buffer zone on the right side
//            size_t idx_s2 = IX(i_end - no_buffer_cell + 1, j, k, Nx, Ny);
            if (expansion_end == ADTypes::UNKNOWN &&
                (*(data + IX(i_end - no_buffer_cell + 1, j, k, Nx, Ny)) > checkValue)) {
                expansion_end = ADTypes::YES;
            }
            if (expansion_end == expansion_start) {
                expansion = expansion_end;
            }
        }
    }
    if (expansion_end == ADTypes::YES) {
        *p_shift_x2 = 1;
    }
    if (expansion_start == ADTypes::YES) {
        *p_shift_x1 = 1;
    }
}
