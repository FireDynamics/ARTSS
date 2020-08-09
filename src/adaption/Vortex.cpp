/// \file       Vortex.cpp
/// \brief      Adaption class for initial condition with vortex
/// \date       Dec 04, 2018
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Vortex.h"
#include "../utility/Parameters.h"
#include "../Domain.h"

Vortex::Vortex(ISolver *solver) {
    auto domain = Domain::getInstance();
    auto params = Parameters::getInstance();
    m_u_lin = params->get_real("initial_conditions/u_lin");
    m_v_lin = params->get_real("initial_conditions/v_lin");
    m_w_lin = params->get_real("initial_conditions/w_lin");
    m_minimal = static_cast<size_t> (std::pow(2, domain->get_levels()));
    m_reduction = (params->get("adaption/class/reduction/enabled") == "Yes");
    if (m_reduction) {
        std::string dir = (params->get("adaption/class/reduction/dir"));
        if (dir.find('x') != std::string::npos) m_x_side = true;
        if (dir.find('y') != std::string::npos) m_y_side = true;
        if (dir.find('z') != std::string::npos) m_z_side = true;
    }
    m_buffer = params->get_int("adaption/class/buffer");
    m_threshold = m_u_lin * params->get_real("adaption/class/threshold");

    u = solver->u;
    v = solver->v;
    w = solver->w;
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
bool Vortex::update(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2) {
    auto d_u = u->data;
    auto d_v = v->data;
    bool adaption = false;

    *p_shift_x1 = 0;
    *p_shift_x2 = 0;
    *p_shift_y1 = 0;
    *p_shift_y2 = 0;
    *p_shift_z1 = 0;
    *p_shift_z2 = 0;

    adaption = Adaption::adapt_x_direction(d_u, m_u_lin, m_buffer, m_threshold, p_shift_x1, p_shift_x2, m_minimal, m_reduction) || adaption;
    if (m_y_side)
        adaption = Adaption::adapt_y_direction(d_v, m_v_lin, m_buffer, m_threshold, p_shift_y1, p_shift_y2, m_minimal, m_reduction) || adaption;

    *p_shift_x1 *= m_minimal;
    *p_shift_x2 *= m_minimal;
    *p_shift_y1 *= m_minimal;
    *p_shift_y2 *= m_minimal;
    *p_shift_z1 *= m_minimal;
    *p_shift_z2 *= m_minimal;

    return adaption;
}

// ==================================== Apply changes =============================
// ********************************************************************************
/// \brief  Set values for new domain
// ********************************************************************************
void Vortex::apply_changes(long *p_shift_x1, long *p_shift_x2, long *p_shift_y1, long *p_shift_y2, long *p_shift_z1, long *p_shift_z2) {
    auto domain = Domain::getInstance();

    size_t i_start = domain->get_index_x1();//(x1 - X1) / dx;
    size_t i_end = domain->get_index_x2() + 2;//(x2 - X1) / dx + 2;
    size_t j_start = domain->get_index_y1();//(y1 - Y1) / dy;
    size_t j_end = domain->get_index_y2() + 2;//(y2 - Y1) / dy + 2;
    size_t k_start = domain->get_index_z1();//(z1 - Z1) / dz;
    size_t k_end = domain->get_index_z2() + 2;//(z2 - Z1) / dz + 2;

    if (*p_shift_x1 != 0) {
        if (*p_shift_x1 < 0) {
            size_t len_e = (j_end - j_start) * (k_end - k_start) * static_cast<size_t>(fabs(*p_shift_x1));
            auto *arr_idxExpansion = new size_t[len_e];
#pragma acc enter data create(arr_idxExpansion[:len_e])
            Adaption::expand_x_direction(*p_shift_x1, true, arr_idxExpansion, len_e);
            Vortex::Drift_dynamic(arr_idxExpansion, len_e);
#pragma acc exit data delete(arr_idxExpansion[:len_e])
            delete[] arr_idxExpansion;
        }
#ifndef BENCHMARKING
        else {
            size_t len_r = (j_end - j_start) * (k_end - k_start) * static_cast<size_t>(fabs(*p_shift_x1));
            auto *arr_idxReduction = new size_t[len_r];
#pragma acc enter data create(arr_idxReduction[:len_r])
            Adaption::reduce_x_direction(*p_shift_x1, true, arr_idxReduction, len_r);
            Zero(arr_idxReduction, len_r);
#pragma acc exit data delete(arr_idxReduction[:len_r])
            delete[] arr_idxReduction;
        }
#endif
    }

    if (*p_shift_x2 != 0) {
        if (*p_shift_x2 > 0) {
            size_t len_e = (j_end - j_start) * (k_end - k_start) * static_cast<size_t>(fabs(*p_shift_x2));
            auto *arr_idxExpansion = new size_t[len_e];
#pragma acc enter data create(arr_idxExpansion[:len_e])
            Adaption::expand_x_direction(*p_shift_x2, false, arr_idxExpansion, len_e);
            Vortex::Drift_dynamic(arr_idxExpansion, len_e);
#pragma acc exit data delete(arr_idxExpansion[:len_e])
            delete[] arr_idxExpansion;
        }
#ifndef BENCHMARKING
        else {
            size_t len_r = (j_end - j_start) * (k_end - k_start) * static_cast<size_t>(fabs(*p_shift_x2));
            auto *arr_idxReduction = new size_t[len_r];
#pragma acc enter data create(arr_idxReduction[:len_r])
            Adaption::reduce_x_direction(*p_shift_x2, false, arr_idxReduction, len_r);
            Zero(arr_idxReduction, len_r);
#pragma acc exit data delete(arr_idxReduction[:len_r])
            delete[] arr_idxReduction;
        }
#endif
    }

    if (*p_shift_y1 != 0) {
        if (*p_shift_y1 < 0) {
            size_t len_e = (i_end - i_start) * (k_end - k_start) * static_cast<size_t>(fabs(*p_shift_y1));
            auto *arr_idxExpansion = new size_t[len_e];
#pragma acc enter data create(arr_idxExpansion[:len_e])
            Adaption::expand_y_direction(*p_shift_y1, true, arr_idxExpansion, len_e);
            Vortex::Drift_dynamic(arr_idxExpansion, len_e);
#pragma acc exit data delete(arr_idxExpansion[:len_e])
            delete[] arr_idxExpansion;
        }
#ifndef BENCHMARKING
        else {
            size_t len_r = (i_end - i_start) * (k_end - k_start) * static_cast<size_t> (fabs(*p_shift_y1));
            auto *arr_idxReduction = new size_t[len_r];
#pragma acc enter data create(arr_idxReduction[:len_r])
            Adaption::reduce_y_Direction(*p_shift_y1, true, arr_idxReduction, len_r);
            Zero(arr_idxReduction, len_r);
#pragma acc exit data delete(arr_idxReduction[:len_r])
            delete[] arr_idxReduction;
        }
#endif
    }

    if (*p_shift_y2 != 0) {
        if (*p_shift_y2 > 0) {
            size_t len_e = (i_end - i_start) * (k_end - k_start) * static_cast<size_t> (fabs(*p_shift_y2));
            auto *arr_idxExpansion = new size_t[len_e];
#pragma acc enter data create(arr_idxExpansion[:len_e])
            Adaption::expand_y_direction(*p_shift_y2, false, arr_idxExpansion, len_e);
            Vortex::Drift_dynamic(arr_idxExpansion, len_e);
#pragma acc exit data delete(arr_idxExpansion[:len_e])
            delete[] arr_idxExpansion;
        }
#ifndef BENCHMARKING
        else {
            size_t len_r = (i_end - i_start) * (k_end - k_start) * static_cast<size_t> (fabs(*p_shift_y2));
            auto *arr_idxReduction = new size_t[len_r];
#pragma acc enter data create(arr_idxReduction[:len_r])
            Adaption::reduce_y_Direction(*p_shift_y2, false, arr_idxReduction, len_r);
            Zero(arr_idxReduction, len_r);
#pragma acc exit data delete(arr_idxReduction[:len_r])
            delete[] arr_idxReduction;
        }
#endif
    }
}

// ==================================== Zero =============================
// ********************************************************************************
/// \brief  Set values to zero (parallel version of function zero)
/// \param  arr_idx Index list of cell which should be set to zero
/// \param  arr_idx_size Size of arr_idx
// ********************************************************************************
void Vortex::Zero(size_t *arr_idx, size_t arr_idx_size) {
    auto data_u = u->data;
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
    auto data_x = u->data;
    auto data_y = v->data;
    auto data_z = w->data;
    // inner cells

#pragma acc parallel loop independent present(data_x[:arr_idx_size], data_z[:arr_idx_size], data_y[:arr_idx_size], arr_idx[:arr_idx_size])
    for (size_t idx = 0; idx < arr_idx_size; idx++) {
        *(data_x + *(arr_idx + idx)) = m_u_lin;
        *(data_y + *(arr_idx + idx)) = m_v_lin;
        *(data_z + *(arr_idx + idx)) = m_w_lin;
    }
}
