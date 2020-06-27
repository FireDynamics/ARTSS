/// \file 		Layers.h
/// \brief 		Adaption class for initial condition with layers (layersT)
/// \date 		Dec 04, 2018
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Layers.h"
#include "Adaption.h"
#include <chrono>
#include "../Domain.h"
#include "../utility/Parameters.h"

Layers::Layers(Adaption *pAdpation, Field ** fields) {
    m_pAdaption = pAdpation;
    auto params = Parameters::getInstance();
    auto domain = Domain::getInstance();
    m_minimal = static_cast<size_t> (std::pow(2, domain->GetLevels()));
    m_timecounter = 0;

    m_noBufferCells = static_cast<size_t> (params->get_int("adaption/class/buffer"));//2;
    m_checkValue = params->get_real("adaption/class/check_value");//temperature;
    m_timestep = static_cast<size_t> (params->get_int("adaption/class/timestep"));//1;
    m_expansion_size = static_cast<size_t> (params->get_int("adaption/class/expansion_size"));

    m_T = fields[VectorFieldsTypes::TEMPERATURE];
    m_Ta = fields[VectorFieldsTypes::TEMPERATURE_A];

    m_u = fields[VectorFieldsTypes::VEL_U];
    m_v = fields[VectorFieldsTypes::VEL_V];
    m_w = fields[VectorFieldsTypes::VEL_W];

    m_Nu = fields[VectorFieldsTypes::NU_T];

    m_P = fields[VectorFieldsTypes::PRESSURE];
    m_kappa = fields[VectorFieldsTypes::KAPPA_T];
    m_gamma = fields[VectorFieldsTypes::GAMMA_T];

    //if (params->get("adaption/version") == "CPU"){
    //    m_fctP_adapt = adaptXDirection_serial;
    //}else{
    //  m_fctP_adapt = adaptXDirection;
    //}
}


 // ==================================== Update ====================================
 // ********************************************************************************
 /// \brief  Checks for adaption
 /// \return  bool if adaption is possible true
 // ********************************************************************************
bool Layers::update() {
    m_timecounter++;
    if (m_timecounter < m_timestep) {
        return false;
    } else {
        m_timecounter = 0;
    }
    auto domain = Domain::getInstance();

    m_pAdaption->m_shift_x1=0;
    m_pAdaption->m_shift_x2=0;
    m_pAdaption->m_shift_y1=0;
    m_pAdaption->m_shift_y2=0;
    m_pAdaption->m_shift_z1=0;
    m_pAdaption->m_shift_z2=0;

    m_x1 = domain->Getx1();
    m_x2 = domain->Getx2();
    m_nx = domain->Getnx();

    m_y1 = domain->Gety1();
    m_y2 = domain->Gety2();
    m_ny = domain->Getny();

    m_z1 = domain->Getz1();
    m_z2 = domain->Getz2();
    m_nz = domain->Getnz();

    adaptXDirection(m_checkValue, m_noBufferCells);
    //TODO z-direction

    size_t adaption =  m_pAdaption->m_shift_x1 + m_pAdaption->m_shift_x2 + m_pAdaption->m_shift_z1 +m_pAdaption->m_shift_z2;
    if (adaption > 0) {
        //TODO z-direction
        m_pAdaption->m_shift_x1 = (-m_minimal * getExpansionSize());
        m_pAdaption->m_shift_x2 = (m_minimal * getExpansionSize());

        // boundary check
        long len_x1 = domain->GetIndexx1()-1;
        long len_x2 = static_cast<long> (std::round((domain->GetX2() - domain->Getx2())/domain->Getdx()));
        if (len_x1 < -m_pAdaption->m_shift_x1 || len_x2 < m_pAdaption->m_shift_x2) {
            if (len_x1 < -m_pAdaption->m_shift_x1) {
                m_pAdaption->m_shift_x2 = (m_pAdaption->m_shift_x2 - m_pAdaption->m_shift_x1 - len_x1);
                m_pAdaption->m_shift_x1 = (-len_x1);
                if (len_x2 < m_pAdaption->m_shift_x2) {
                    m_pAdaption->m_shift_x2 = len_x2;
                }
            }
            if (len_x2 < m_pAdaption->m_shift_x2) {
                m_pAdaption->m_shift_x1 = m_pAdaption->m_shift_x1 - (m_pAdaption->m_shift_x2 - len_x2);
                m_pAdaption->m_shift_x2 = len_x2;
                if (len_x1 < -m_pAdaption->m_shift_x1) {
                    m_pAdaption->m_shift_x1 = (-len_x1);
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
void Layers::setXValues(bool start) {
    Domain *domain = Domain::getInstance();

    size_t Nx = domain->GetNx();
    size_t Ny = domain->GetNy();

    size_t j_start = static_cast<size_t> (std::round((m_y1 - domain->GetY1()) / domain->Getdy()));
    size_t j_end = j_start + m_ny;
    size_t k_start = static_cast<size_t> (std::round((m_z1 - domain->GetZ1()) / domain->Getdz()));
    size_t k_end = k_start + m_nz;

    real *data_nu = m_Nu->data;
    real *data_gamma = m_gamma->data;
    real *data_kappa = m_kappa->data;
    real *data_temp = m_T->data;
    real *data_tempA = m_Ta->data;

    size_t size = domain->GetSize();

    if (start) {
        size_t nx_begin = static_cast<size_t> (std::round((m_x1 - domain->GetX1()) / domain->Getdx()));
        long shift = m_pAdaption->m_shift_x1;
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
        size_t nx_end = static_cast<size_t> (std::round((m_x2 - domain->GetX1()) / domain->Getdx()));
        long shift = m_pAdaption->m_shift_x2;
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
void Layers::applyChanges() {
    if (m_pAdaption->m_shift_x1) {
        Layers::setXValues(true);
    }
    if (m_pAdaption->m_shift_x2) {
        Layers::setXValues(false);
    }
}

 // ==================================== Adaption x direction parallel ===============================
 // ***************************************************************************************
 /// \brief  Checks if adaption is possible and allowed
 /// \param  checkValue check value
 /// \param  noBufferCell Buffersize
 // ***************************************************************************************
void Layers::adaptXDirection(real checkValue, size_t noBufferCell) {
    auto domain = Domain::getInstance();

    auto data = m_T->data;
    size_t size = domain->GetSize();

    size_t expansion_counter_start = 0;
    size_t expansion_counter_end = 0;

    bool expansion_start = (domain->Getx1() != domain->GetX1());
    bool expansion_end = (domain->Getx2() != domain->GetX2());

    //copy(expansion_counter_start,expansion_counter_end)
#pragma acc data present(data[:size]) copyout(expansion_counter_end, expansion_counter_start)
    {
    size_t Nx = domain->GetNx();
    size_t Ny = domain->GetNy();

    size_t nx_begin = domain->GetIndexx1();//((domain->Getx1() - domain->GetX1()) / domain->Getdx());

    size_t j_start = domain->GetIndexy1();//((domain->Gety1() - domain->GetY1()) / domain->Getdy());
    size_t j_end = j_start + domain->Getny() - 1;
    size_t k_start = domain->GetIndexz1();//((domain->Getz1() - domain->GetZ1()) / domain->Getdz());
    size_t k_end = k_start + domain->Getnz() - 1;

    size_t i_start = nx_begin;
    size_t i_end = nx_begin + domain->Getnx() - 1;

    //loop through left and right side of cuboid in x direction
#pragma acc parallel loop collapse(2) present(data[:size]) reduction(+:expansion_counter_end)
    for (size_t j = j_start; j < j_end; j++) {
        for (size_t k = k_start; k < k_end; k++) {
            // check innermost plane of the buffer zone on the right side
            if ((*(data + IX(i_end - noBufferCell + 1, j, k, Nx, Ny)) > checkValue)) {
                expansion_counter_end++;
            }
        }
    }
//#pragma acc parallel loop collapse(2) present(data[:size]) reduction(max:expansion_counter_start)
#pragma acc parallel loop collapse(2) present(data[:size]) reduction(+:expansion_counter_start)
    for (size_t j = j_start; j < j_end; j++) {
        for (size_t k = k_start; k < k_end; k++) {
            // check innermost plane of the buffer zone on the left side
            if ((*(data + IX(i_start + noBufferCell - 1, j, k, Nx, Ny)) > checkValue)) {
                expansion_counter_start++;
            }
        }
    }
    }
    if (expansion_counter_end > 0 && expansion_end) {
        m_pAdaption->m_shift_x2 = 1;
    }
    if (expansion_counter_start > 0 && expansion_start) {
        m_pAdaption->m_shift_x1 = 1;
    }
}

// ==================================== Adaption x direction serial ===============================
// ***************************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param  checkValue check value
/// \param  noBufferCell Buffersize
// ***************************************************************************************
void Layers::adaptXDirection_serial(real checkValue, size_t noBufferCell) {
    auto domain = Domain::getInstance();
    size_t Nx = domain->GetNx();
    size_t Ny = domain->GetNy();

    size_t nx_begin = domain->GetIndexx1();//((domain->Getx1() - domain->GetX1()) / domain->Getdx());

    size_t j_start = domain->GetIndexy1();//((domain->Gety1() - domain->GetY1()) / domain->Getdy());
    size_t j_end = j_start + domain->Getny() - 1;
    size_t k_start = domain->GetIndexz1();//((domain->Getz1() - domain->GetZ1()) / domain->Getdz());
    size_t k_end = k_start + domain->Getnz() - 1;

    size_t i_start = nx_begin;
    size_t i_end = nx_begin + domain->Getnx() - 1;

    //expansion - expand if there is at least one cell in the buffer area fulfills the condition
    ADTypes expansion_start = ADTypes::UNKNOWN;
    if (domain->Getx1() == domain->GetX1()) {
        expansion_start = ADTypes::NO;
    }
    ADTypes expansion_end = ADTypes::UNKNOWN;
    if (domain->Getx2() == domain->GetX2()) {
        expansion_end = ADTypes::NO;
    }
    ADTypes expansion = ADTypes::UNKNOWN;
    if (expansion_end == expansion_start) {
        expansion = expansion_end;
    }

    auto data = m_T->data;
    size_t size = domain->GetSize();

    //loop through left and right side of cuboid in x direction
    for (size_t j = j_start; j < j_end && expansion == ADTypes::UNKNOWN; j++) {
        for (size_t k = k_start; k < k_end && expansion == ADTypes::UNKNOWN; k++) {
            // check innermost plane of the buffer zone on the left side
            //size_t idx_s1 = IX(i_start + noBufferCell - 1, j, k, Nx, Ny);
            if (expansion_start == ADTypes::UNKNOWN &&
                (*(data + IX(i_start + noBufferCell - 1, j, k, Nx, Ny)) > checkValue)) {
                expansion_start = ADTypes::YES;
            }
            // check innermost plane of the buffer zone on the right side
//            size_t idx_s2 = IX(i_end - noBufferCell + 1, j, k, Nx, Ny);
            if (expansion_end == ADTypes::UNKNOWN &&
                (*(data + IX(i_end - noBufferCell + 1, j, k, Nx, Ny)) > checkValue)) {
                expansion_end = ADTypes::YES;
            }
            if (expansion_end == expansion_start) {
                expansion = expansion_end;
            }
        }
    }
    if (expansion_end == ADTypes::YES) {
        m_pAdaption->m_shift_x2 = 1;
    }
    if (expansion_start == ADTypes::YES) {
        m_pAdaption->m_shift_x1 = 1;
    }
}
