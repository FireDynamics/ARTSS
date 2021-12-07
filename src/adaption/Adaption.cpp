/// \file       Adaption.cpp
/// \brief      Controller class for adaption
/// \date       Nov 29, 2018
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.


#ifdef _OPENACC
#include <accelmath.h>
#endif

#include <chrono>
#include <fstream>

#include "Adaption.h"
#include "Layers.h"
#include "Vortex.h"
#include "../DomainData.h"
#include "../boundary/BoundaryController.h"

Adaption::Adaption(Settings::Settings const &settings, FieldController *field_controller) {
    m_field_controller = field_controller;
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    auto domain_data = DomainData::getInstance();
    m_dynamic = settings.get_bool("adaption/dynamic");
    m_filename = settings.get_filename();
    m_filename.resize(m_filename.size() - 4);//remove .xml from filename
    m_has_data_extraction = settings.get_bool("adaption/data_extraction");
    if (m_has_data_extraction) {
        m_has_data_extraction_before = settings.get_bool("adaption/data_extraction/before/enabled");
        m_has_data_extraction_after = settings.get_bool("adaption/data_extraction/after/enabled");
        m_has_data_extraction_endresult = settings.get_bool("adaption/data_extraction/endresult/enabled");
        m_has_time_measuring = settings.get_bool("adaption/data_extraction/time_measuring/enabled");
        //m_has_write_runtime = (settings.get("adaption/data_extraction/runtime/enabled") == "Yes");
    }
    if (m_dynamic) {
        std::string init = settings.get("adaption/class/name");
        if (init == "Layers") {
            func = new Layers(settings, m_field_controller);
        } else if (init == "Vortex" || init == "VortexY") {
            func = new Vortex(settings, m_field_controller);
        } else {
#ifndef BENCHMARKING
            m_logger->critical("Type {} is not defined", init);
#endif
            std::exit(1);
            //TODO Error Handling
        }

        m_dynamic_end = false;
        m_minimal = static_cast<size_t> (std::pow(2, domain_data->get_levels()));
        m_shift_x1 = 0;
        m_shift_x2 = 0;
        m_shift_y1 = 0;
        m_shift_y2 = 0;
        m_shift_z1 = 0;
        m_shift_z2 = 0;
    }
}

// ==================================== Run ====================================
// ***************************************************************************************
/// \brief  starts adaption process
/// \param  t_cur current timestep
// ***************************************************************************************
void Adaption::run(real) {
    if (m_dynamic && isUpdateNecessary()) {
        applyChanges();
#ifndef BENCHMARKING
        std::ofstream file;
        file.open(get_time_measuring_name(), std::ios::app);
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
#endif
        BoundaryController::getInstance()->update_lists();
#ifndef BENCHMARKING
        end = std::chrono::system_clock::now();
        if (m_has_time_measuring) {
            long ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            file << "refresh: " << ms << " microsec\n";
        }
        file.close();
#endif
    }
}

// ==================================== Extract data ====================================
// ***************************************************************************************
/// \brief  takes 2D slice of 3D mesh at position y = height and writes slice into file
/// \param  filename  filename
/// \param  height y-value
/// \param  time timestep
// ***************************************************************************************
void Adaption::extractData(const std::string &filename, real height, real time) {
    auto domain_data = DomainData::getInstance();
    size_t x_start = 0;
    size_t x_end = domain_data->get_Nx();
    auto y = static_cast<size_t>(std::round((height - domain_data->get_Y1()) / domain_data->get_dy()));
    auto z = static_cast<size_t >(std::round(domain_data->get_Nz() / 2));
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real *data_temp = m_field_controller->get_field_T().data;
    real *data_tempA = m_field_controller->get_field_T_ambient().data;
    real *data_u = m_field_controller->get_field_u().data;
    real *data_v = m_field_controller->get_field_v().data;
    real *data_w = m_field_controller->get_field_w().data;
    real *data_nu = m_field_controller->get_field_nu_t().data;
    real *data_kappa = m_field_controller->get_field_kappa().data;

    std::ofstream file;
    file.open(filename, std::ios::app);
    file << time << ";";
    for (size_t x = x_start; x < x_end; x++) {
        size_t idx = IX(x, y, z, Nx, Ny);
        file << x << "|" << *(data_temp + idx) << "|" << *(data_tempA + idx) << "|" << *(data_u + idx) << "|"
             << *(data_v + idx) << "|" << *(data_w + idx) << "|" << *(data_nu + idx) << "|" << *(data_kappa + idx)
             << ";";
    }
    file << std::endl;
    file.close();
}

// ==================================== Extract data ====================================
// ***************************************************************************************
/// \brief  takes all cells and writes into a file
/// \param  filename  filename
// ***************************************************************************************
void Adaption::extractData(const std::string &filename) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();
    size_t x_end = Nx;
    size_t y_end = Ny;
    size_t z_end = domain_data->get_Nz();

    real *data_temp = m_field_controller->get_field_T().data;
    real *data_u = m_field_controller->get_field_u().data;
    real *data_v = m_field_controller->get_field_v().data;
    real *data_w = m_field_controller->get_field_w().data;

    std::ofstream file;
    file.open(filename, std::ios::app);
    for (size_t z = 1; z < z_end - 1; z++) {
        for (size_t y = 1; y < y_end - 1; y++) {
            for (size_t x = 1; x < x_end - 1; x++) {
                size_t idx = IX(x, y, z, Nx, Ny);
                file << x << "|" << y << "|" << z << "|" << *(data_temp + idx) << "|" << *(data_u + idx) << "|"
                     << *(data_v + idx) << "|" << *(data_w + idx) << ";";
            }
            file << std::endl;
        }
    }
    file << std::endl;
    file.close();
}


// ==================================== Apply changes ====================================
// ***************************************************************************************
/// \brief  Applies domain adaption
// ***************************************************************************************
void Adaption::applyChanges() {
#ifndef BENCHMARKING
    std::ofstream file;
    file.open(get_time_measuring_name(), std::ios::app);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
#endif
    auto domain = DomainData::getInstance();
    auto shift_start = new Coordinate<long>(m_shift_x1, m_shift_y1, m_shift_z1);
    auto shift_end = new Coordinate<long>(m_shift_x2, m_shift_y2, m_shift_z2);
    if (domain->resize(*shift_start, *shift_end)) {
        if (domain->get_X1() == domain->get_x1() &&
            domain->get_X2() == domain->get_x2() &&
            domain->get_Y1() == domain->get_y1() &&
            domain->get_Y2() == domain->get_y2() &&
            domain->get_Z1() == domain->get_z1() &&
            domain->get_Z2() == domain->get_z2()) {
            m_dynamic_end = true;
        }
        func->apply_changes(&m_shift_x1, &m_shift_x2, &m_shift_y1, &m_shift_y2, &m_shift_z1, &m_shift_z2);
    }
#ifndef BENCHMARKING
    end = std::chrono::system_clock::now();
    long ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    if (m_has_time_measuring) {
        file << "apply: " << ms << " microsec\n";
    }
    file.close();
#endif
}

// ==================================== Is update necessary ==============================
// ***************************************************************************************
/// \brief  Checks if adaption should be done
/// \return bool true if yes false if no
// ***************************************************************************************
bool Adaption::isUpdateNecessary() {
#ifndef BENCHMARKING
    std::ofstream file;
    file.open(get_time_measuring_name(), std::ios::app);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    if (m_has_time_measuring) {
        file.open(get_time_measuring_name(), std::ios::app);
        start = std::chrono::system_clock::now();
    }
#endif
    bool update = false;
    if (!m_dynamic_end) {
        update = func->update(&m_shift_x1, &m_shift_x2, &m_shift_y1, &m_shift_y2, &m_shift_z1, &m_shift_z2);
    }
#ifndef BENCHMARKING
    end = std::chrono::system_clock::now();
    long ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    if (m_has_time_measuring) {
        file << "update: " << ms << " microsec\n";
        file.close();
    }
#endif
    return update;
}

// ==================================== Expand x direction ===============================
// ***************************************************************************************
/// \brief  Basic implementation for expansion in x direction (parallelized)
/// \param  shift expansion size value in x direction
/// \param  start indicates whether the expansion is at the beginning or the end of the computational domain
/// \param  arr_idx_expansion  Index list of cells to be newly added
/// \param len_e  size of arr_idx_expansion
// ***************************************************************************************
void Adaption::expand_x_direction(long shift, bool start, size_t *arr_idx_expansion, size_t len_e __attribute__((unused))) {
    auto domain = DomainData::getInstance();
#pragma acc data present(arr_idx_expansion[:len_e])
    {
        size_t j_start = domain->get_index_y1();
        size_t j_end = domain->get_index_y2() + 1;

        size_t k_start = domain->get_index_z1();
        size_t k_end = domain->get_index_z2() + 1;

        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        if (start) { // shift = negative
            size_t i1 = domain->get_index_x1();
            shift *= (-1);
#pragma acc parallel loop collapse(3) present(arr_idx_expansion[:len_e])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t j = j_start; j < j_end; j++) {
                    for (long ii = 0; ii < shift; ii++) {
                        size_t idx = IX(i1 + ii, j, k, Nx, Ny);
                        size_t counter = (k - k_start) * (shift) * (j_end - j_start) + (j - j_start) * (shift) + ii;
                        *(arr_idx_expansion + counter) = idx;
                    }
                }
            }
        } else { // shift = positive
            size_t i2 = domain->get_index_x2();
#pragma acc parallel loop collapse(3) present(arr_idx_expansion[:len_e])
            for (size_t j = j_start; j < j_end; j++) {
                for (size_t k = k_start; k < k_end; k++) {
                    for (long ii = 0; ii < shift; ii++) {
                        size_t idx = IX(i2 - ii, j, k, Nx, Ny);
                        size_t counter = (k - k_start) * (shift) * (j_end - j_start) + (j - j_start) * (shift) + ii;
                        *(arr_idx_expansion + counter) = idx;
                    }
                }
            }
        }
    }
}

// ==================================== Expand y direction ===============================
// ***************************************************************************************
/// \brief  Basic implementation for expansion in y direction (parallelized)
/// \param  shift expansion size value in y direction
/// \param  start indicates whether the expansion is at the beginning or the end of the computational domain
/// \param  arr_idx_expansion  Index list of cells to be newly added
/// \param len_e  size of arr_idx_expansion
// ***************************************************************************************
void Adaption::expand_y_direction(long shift, bool start, size_t *arr_idx_expansion, size_t len_e __attribute__((unused))) {
    auto domain = DomainData::getInstance();

#pragma acc data present(arr_idx_expansion[:len_e])
    {
        size_t i_start = domain->get_index_x1();
        size_t i_end = domain->get_index_x2() + 1;
        size_t k_start = domain->get_index_z1();
        size_t k_end = domain->get_index_z2() + 1;
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        if (start) { // shift = negative
            size_t j1 = domain->get_index_y1();
            shift *= (-1);
#pragma acc parallel loop collapse(3) present(arr_idx_expansion[:len_e])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t i = i_start; i < i_end; i++) {
                    for (long jj = 0; jj < shift; jj++) {
                        size_t idx = IX(i, j1 + jj, k, Nx, Ny);
                        size_t counter = (k - k_start) * (shift) * (i_end - i_start) + (i - i_start) * (shift) + jj;
                        *(arr_idx_expansion + counter) = idx;
                    }
                }
            }
        } else { // shift = positive
            size_t j2 = domain->get_index_y2();
#pragma acc parallel loop collapse(3) present(arr_idx_expansion[:len_e])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t i = i_start; i < i_end; i++) {
                    for (long jj = 0; jj < shift; jj++) {
                        size_t idx = IX(i, j2 - jj, k, Nx, Ny);
                        size_t counter = (k - k_start) * (shift) * (i_end - i_start) + (i - i_start) * (shift) + jj;
                        *(arr_idx_expansion + counter) = idx;
                    }
                }
            }
        }
    }
}

// ==================================== Reduction x direction ===============================
// ***************************************************************************************
/// \brief  Basic implementation for expansion in x direction (parallelized)
/// \param  shift reduction size value in x direction
/// \param  start indicates whether the reduction is at the beginning or the end of the computational domain
/// \param  arr_idx_reduction  Index list of cells to be newly added
/// \param len_e  size of arr_idx_reduction
// ***************************************************************************************
void Adaption::reduce_x_direction(long shift, bool start, size_t *arr_idx_reduction, size_t len_r __attribute__((unused))) {
    auto domain = DomainData::getInstance();
#pragma acc data present(arr_idx_reduction[:len_r])
    {
        size_t j_start = domain->get_index_y1();
        size_t j_end = domain->get_index_y2() + 1;
        size_t k_start = domain->get_index_z1();
        size_t k_end = domain->get_index_z2() + 1;
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        if (start) { // shift = positive
            size_t i1 = domain->get_index_x1();
#pragma acc parallel loop collapse(3) present(arr_idx_reduction[:len_r])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t j = j_start; j < j_end; j++) {
                    for (long ii = 0; ii < shift; ii++) {
                        size_t counter = (k - k_start) * (shift) * (j_end - j_start) + (j - j_start) * (shift) + ii;
                        auto idx = IX(i1 - ii - 1, j, k, Nx, Ny);
                        *(arr_idx_reduction + counter) = idx;
                    }
                }
            }
        } else { // shift = negative
            shift *= (-1);
            size_t i2 = domain->get_index_x2();
#pragma acc parallel loop collapse(3) present(arr_idx_reduction[:len_r])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t j = j_start; j < j_end; j++) {
                    for (long ii = 0; ii < shift; ii++) {
                        size_t counter = (k - k_start) * (shift) * (j_end - j_start) + (j - j_start) * (shift) + ii;
                        auto idx = IX(i2 + ii + 2, j, k, Nx, Ny);
                        *(arr_idx_reduction + counter) = idx;
                    }
                }
            }
        }
    }
}

// ==================================== Reduction y direction ===============================
// ***************************************************************************************
/// \brief  Basic implementation for expansion in y direction (parallelized)
/// \param  shift reduction size value in y direction
/// \param  start indicates whether the reduction is at the beginning or the end of the computational domain
/// \param  arr_idx_reduction  Index list of cells to be newly added
/// \param len_e  size of arr_idx_reduction
// ***************************************************************************************
void Adaption::reduce_y_Direction(long shift, bool start, size_t *arr_idx_reduction, size_t len_r __attribute__((unused))) {
    auto domain = DomainData::getInstance();
#pragma acc data present(arr_idx_reduction[:len_r])
    {
        size_t i_start = domain->get_index_x1();
        size_t i_end = domain->get_index_x2() + 2;
        size_t k_start = domain->get_index_z1();
        size_t k_end = domain->get_index_z2() + 2;

        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        if (start) { // shift = positive
            size_t j1 = domain->get_index_y1();
#pragma acc parallel loop collapse(3) present(arr_idx_reduction[:len_r])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t i = i_start; i < i_end; i++) {
                    for (long jj = 0; jj < shift; jj++) {
                        size_t counter = (k - k_start) * (shift) * (i_end - i_start) + (i - i_start) * (shift) + jj;
                        auto idx = IX(i, j1 - jj - 1, k, Nx, Ny);
                        *(arr_idx_reduction + counter) = idx;
                    }
                }
            }
        } else { // shift = negative
            shift *= (-1);
            size_t j2 = domain->get_index_y2();
#pragma acc parallel loop collapse(3) present(arr_idx_reduction[:len_r])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t i = i_start; i < i_end; i++) {
                    for (long jj = 0; jj < shift; jj++) {
                        size_t counter = (k - k_start) * (shift) * (i_end - i_start) + (i - i_start) * (shift) + jj;
                        auto idx = IX(i, j2 + jj + 2, k, Nx, Ny);
                        *(arr_idx_reduction + counter) = idx;
                    }
                }
            }
        }
    }
}

// ==================================== Adaption x direction serial ===============================
// ***************************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param  f field
/// \param  check_value check value
/// \param  no_buffer_cell Buffersize
/// \param  threshold precision of comparison
// ***************************************************************************************
bool Adaption::adapt_x_direction_serial(const real *f, real check_value, size_t no_buffer_cell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce) {
    auto domain = DomainData::getInstance();
    size_t Nx = domain->get_Nx();
    size_t Ny = domain->get_Ny();

    size_t i_start = domain->get_index_x1();
    size_t i_end = domain->get_index_x2() + 1;
    size_t j_start = domain->get_index_y1();
    size_t j_end = domain->get_index_y2() + 1;
    size_t k_start = domain->get_index_z1();
    size_t k_end = domain->get_index_z2() + 1;

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

    //reduction - reduce if all cells do not fulfil the condition any longer
    ADTypes reduction_start = ADTypes::NO;
    ADTypes reduction_end = ADTypes::NO;
    ADTypes reduction = ADTypes::NO;
    if (reduce && domain->get_nx() > (minimal + no_buffer_cell) * 2) {
        reduction_start = ADTypes::UNKNOWN;
        reduction_end = ADTypes::UNKNOWN;
        reduction = ADTypes::UNKNOWN;
    }
    bool check = !(expansion == ADTypes::YES || (reduction == ADTypes::NO && expansion == ADTypes::NO));
    //loop through left and right side of cuboid in x direction
    for (size_t j = j_start; j < j_end && check; j++) {
        for (size_t k = k_start; k < k_end && check; k++) {
            // check innermost plane of the buffer zone on the left side
            size_t idx_s1 = IX(i_start + no_buffer_cell - 1, j, k, Nx, Ny);
            if (expansion_start == ADTypes::UNKNOWN && std::fabs(*(f + idx_s1) - check_value) > threshold) {
                *p_shift_x1 = -1;
                expansion_start = ADTypes::YES;
            } else {
                // check innermost plane of the minimal zone to reduce on the left side
                size_t idx = IX(i_start + minimal - 1, j, k, Nx, Ny);
                size_t idx2 = IX(i_start + minimal - 1 + no_buffer_cell, j, k, Nx, Ny);
                if (reduction_start == ADTypes::UNKNOWN && (std::fabs(*(f + idx2) - check_value) > threshold ||
                                                            std::fabs(*(f + idx) - check_value) > threshold)) {
                    reduction_start = ADTypes::NO;
                }
            }
            // check innermost plane of the buffer zone on the right side
            size_t idx_s2 = IX(i_end - no_buffer_cell + 1, j, k, Nx, Ny);
            if (expansion_end == ADTypes::UNKNOWN && std::fabs(*(f + idx_s2) - check_value) > threshold) {
                *p_shift_x2 = 1;
                expansion_end = ADTypes::YES;
            } else {
                // check innermost plane of the minimal zone to reduce on the right side
                size_t idx = IX(i_end - minimal + 1, j, k, Nx, Ny);
                size_t idx2 = IX(i_end - minimal + 1 - no_buffer_cell, j, k, Nx, Ny);
                if (reduction_end == ADTypes::UNKNOWN && (std::fabs(*(f + idx2) - check_value) > threshold ||
                                                          std::fabs(*(f + idx) - check_value) > threshold)) {
                    reduction_end = ADTypes::NO;
                }
            }
            if (expansion_end == expansion_start) {
                expansion = expansion_end;
            }
            if (reduction_end == reduction_start) {
                reduction = reduction_end;
            }
            check = !(expansion == ADTypes::YES || (reduction == ADTypes::NO && expansion == ADTypes::NO));
        }
    }
    if ((expansion_start == reduction_start && expansion_start == ADTypes::YES) ||
        (expansion_end == reduction_end && expansion_end == ADTypes::YES)) {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(Adaption).name());
        m_logger->error("Exception in x-Adaption: {} {} {} {}",
                        size_t(expansion_start),
                        size_t(reduction_start),
                        size_t(expansion_end),
                        size_t(reduction_end));
#endif
        // TODO Error handling
        throw std::exception();
    }
    if (reduction_start == ADTypes::UNKNOWN && expansion_start != ADTypes::YES) {
        reduction_start = ADTypes::YES;
        *p_shift_x1 = 1;
    }
    if (reduction_end == ADTypes::UNKNOWN && expansion_end != ADTypes::YES) {
        reduction_end = ADTypes::YES;
        *p_shift_x2 = -1;
    }
    return expansion_end == ADTypes::YES || expansion_start == ADTypes::YES || reduction_start == ADTypes::YES ||
           reduction_end == ADTypes::YES;
}

// ========================= Adaption x direction parallel ====================
// *****************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param  f field
/// \param  check_value check value
/// \param  no_buffer_cell Buffersize
/// \param  threshold precision of comparison
// ***************************************************************************************
bool Adaption::adapt_x_direction(const real *f, real check_value, size_t no_buffer_cell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce) {
    auto domain = DomainData::getInstance();
    size_t expansion_counter_start = 0;
    size_t expansion_counter_end = 0;
    size_t reduction_counter_start = 0;
    size_t reduction_counter_end = 0;

    bool reduction_start = reduce;
    bool reduction_end = reduce;
    if (reduce && domain->get_nx() > (minimal + no_buffer_cell) * 2) {
        reduction_start = false;
        reduction_end = false;
    }
#pragma acc data present(f) copy(expansion_counter_start, expansion_counter_end, reduction_counter_start, reduction_counter_end)
    {
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        size_t i_start = domain->get_index_x1();
        size_t i_end = domain->get_index_x2() + 1;
        size_t j_start = domain->get_index_y1();
        size_t j_end = domain->get_index_y2() + 1;
        size_t k_start = domain->get_index_z1();
        size_t k_end = domain->get_index_z2() + 1;

        //expansion - expand if there is at least one cell in the buffer area fulfills the condition
        //reduction - reduce if all cells do not fulfil the condition any longer
        //loop through left side of cuboid in x direction
#pragma acc parallel loop collapse(2) present(f) reduction(+:expansion_counter_start, reduction_counter_start, expansion_counter_end, reduction_counter_end)
        for (size_t j = j_start; j < j_end; j++) {
            for (size_t k = k_start; k < k_end; k++) {
                // check innermost plane of the buffer zone on the left side
                size_t idx_s1 = IX(i_start + no_buffer_cell - 1, j, k, Nx, Ny);
                if (std::fabs(*(f + idx_s1) - check_value) > threshold) {
                    expansion_counter_start++;
                } else {
                    size_t idx = IX(i_start + minimal - 1, j, k, Nx, Ny);
                    size_t idx2 = IX(i_start + minimal - 1 + no_buffer_cell, j, k, Nx, Ny);
                    if ((std::fabs(*(f + idx2) - check_value) > threshold ||
                         std::fabs(*(f + idx) - check_value) > threshold)) {
                        reduction_counter_start++;
                    }
                }
                size_t idx_s2 = IX(i_end - no_buffer_cell + 1, j, k, Nx, Ny);
                if (std::fabs(*(f + idx_s2) - check_value) > threshold) {
                    expansion_counter_end++;
                } else {
                    size_t idx = IX(i_end - minimal + 1, j, k, Nx, Ny);
                    size_t idx2 = IX(i_end - minimal + 1 - no_buffer_cell, j, k, Nx, Ny);
                    if ((std::fabs(*(f + idx2) - check_value) > threshold ||
                         std::fabs(*(f + idx) - check_value) > threshold)) {
                        reduction_counter_end++;
                    }
                }
            }
        }
    }
    if ((expansion_counter_start > 0 && reduction_counter_start == 0 && reduction_start) ||
        (expansion_counter_end > 0 && reduction_counter_end == 0 && reduction_end)) {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(Adaption).name());
        m_logger->error("Trying to reduce and expand at the same time (x): {},{} | {},{}",
                expansion_counter_start,
                reduction_counter_start,
                expansion_counter_end,
                reduction_counter_end);
#endif
        //TODO Error Handling
        //throw std::exception();
    }
    if (expansion_counter_start > 0) {
        *p_shift_x1 = -1;
    } else {
        if (reduction_counter_start == 0 && reduce) {
            *p_shift_x1 = 1;
        }
    }
    if (expansion_counter_end > 0) {
        *p_shift_x2 = 1;
    } else {
        if (reduction_counter_end == 0 && reduce) {
            *p_shift_x2 = -1;
        }
    }
    return (expansion_counter_end + expansion_counter_start) > 0 || reduction_counter_start == 0 ||
           reduction_counter_end == 0;
}

// ==================================== Adaption y direction serial ===============================
// ***************************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param  f field
/// \param  check_value check value
/// \param  no_buffer_cell Buffersize
/// \param  threshold precision of comparison
// ***************************************************************************************
bool Adaption::adapt_y_direction_serial(const real *f, real check_value, size_t no_buffer_cell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce) {
    auto domain = DomainData::getInstance();
    size_t Nx = domain->get_Nx();
    size_t Ny = domain->get_Ny();

    size_t i_start = domain->get_index_x1();
    size_t i_end = domain->get_index_x2() + 1;
    size_t j_start = domain->get_index_y1();
    size_t j_end = domain->get_index_y2() + 1;
    size_t k_start = domain->get_index_z1();
    size_t k_end = domain->get_index_z2() + 1;

    //expansion - expand if there is at least one cell in the buffer area fulfills the condition
    ADTypes expansion_start = ADTypes::UNKNOWN;
    if (domain->get_y1() == domain->get_Y1()) {
        expansion_start = ADTypes::NO;
    }
    ADTypes expansion_end = ADTypes::UNKNOWN;
    if (domain->get_y2() == domain->get_Y2()) {
        expansion_end = ADTypes::NO;
    }
    ADTypes expansion = ADTypes::UNKNOWN;
    if (expansion_end == expansion_start) {
        expansion = expansion_end;
    }

    //reduction - reduce if all cells do not fulfil the condition any longer
    ADTypes reduction_start = ADTypes::NO;
    ADTypes reduction_end = ADTypes::NO;
    ADTypes reduction = ADTypes::NO;
    if (reduce && domain->get_nx() > (minimal + no_buffer_cell) * 2) {
        reduction_start = ADTypes::UNKNOWN;
        reduction_end = ADTypes::UNKNOWN;
        reduction = ADTypes::UNKNOWN;
    }
    bool check = !(expansion == ADTypes::YES || (reduction == ADTypes::NO && expansion == ADTypes::NO));
    //loop through left and right side of cuboid in x direction
    for (size_t i = i_start; i < i_end && check; i++) {
        for (size_t k = k_start; k < k_end && check; k++) {
            // check innermost plane of the buffer zone on the lower side
            size_t idx_s1 = IX(i, j_start + no_buffer_cell - 1, k, Nx, Ny);
            if (expansion_start == ADTypes::UNKNOWN && std::fabs(*(f + idx_s1) - check_value) > threshold) {
                *p_shift_x1 = -1;
                expansion_start = ADTypes::YES;
            } else {
                // check innermost plane of the minimal zone to reduce on the lower side
                size_t idx = IX(i, j_start + minimal - 1, k, Nx, Ny);
                size_t idx2 = IX(i, j_start + minimal - 1 + no_buffer_cell, k, Nx, Ny);
                if (reduction_start == ADTypes::UNKNOWN && (std::fabs(*(f + idx2) - check_value) > threshold ||
                                                            std::fabs(*(f + idx) - check_value) > threshold)) {
                    reduction_start = ADTypes::NO;
                }
            }
            // check innermost plane of the buffer zone on the upper side
            size_t idx_s2 = IX(i, j_end - no_buffer_cell + 1, k, Nx, Ny);
            if (expansion_end == ADTypes::UNKNOWN && std::fabs(*(f + idx_s2) - check_value) > threshold) {
                *p_shift_x2 = 1;
                expansion_end = ADTypes::YES;
            } else {
                // check innermost plane of the minimal zone to reduce on the upper side
                size_t idx = IX(i, j_end - minimal + 1, k, Nx, Ny);
                size_t idx2 = IX(i, j_end - minimal + 1 - no_buffer_cell, k, Nx, Ny);
                if (reduction_end == ADTypes::UNKNOWN && (std::fabs(*(f + idx2) - check_value) > threshold ||
                                                          std::fabs(*(f + idx) - check_value) > threshold)) {
                    reduction_end = ADTypes::NO;
                }
            }
            if (expansion_end == expansion_start) {
                expansion = expansion_end;
            }
            if (reduction_end == reduction_start) {
                reduction = reduction_end;
            }
            check = !(expansion == ADTypes::YES || (reduction == ADTypes::NO && expansion == ADTypes::NO));
        }
    }
    if ((expansion_start == reduction_start && expansion_start == ADTypes::YES) ||
        (expansion_end == reduction_end && expansion_end == ADTypes::YES)) {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(Adaption).name());
        m_logger->error("Exception in y-Adaption: {} {} {} {}",
                        size_t(expansion_start),
                        size_t(reduction_start),
                        size_t(expansion_end),
                        size_t(reduction_end));
#endif
        // TODO Error handling
        throw std::exception();
    }
    if (reduction_start == ADTypes::UNKNOWN && expansion_start != ADTypes::YES) {
        reduction_start = ADTypes::YES;
        *p_shift_x1 = 1;
    }
    if (reduction_end == ADTypes::UNKNOWN && expansion_end != ADTypes::YES) {
        reduction_end = ADTypes::YES;
        *p_shift_x2 = -1;
    }
    return expansion_end == ADTypes::YES || expansion_start == ADTypes::YES ||
           reduction_start == ADTypes::YES ||
           reduction_end == ADTypes::YES;
}

// ========================= Adaption x direction parallel ====================
// *****************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param  f field
/// \param  check_value check value
/// \param  no_buffer_cell Buffersize
/// \param  threshold precision of comparison
// ***************************************************************************************
bool Adaption::adapt_y_direction(const real *f, real check_value, size_t no_buffer_cell, real threshold, long *p_shift_x1, long *p_shift_x2, size_t minimal, bool reduce) {
    auto domain = DomainData::getInstance();

    size_t expansion_counter_start = 0;
    size_t expansion_counter_end = 0;
    size_t reduction_counter_start = 0;
    size_t reduction_counter_end = 0;
    //reduction - reduce if all cells do not fulfil the condition any longer
    bool reduction_start = reduce;
    bool reduction_end = reduce;
    if (reduce && domain->get_ny() > (minimal + no_buffer_cell) * 2) {
        reduction_start = false;
        reduction_end = false;
    }
#pragma acc data present(f) copy(expansion_counter_start, expansion_counter_end, reduction_counter_start, reduction_counter_end)
    {
        size_t Nx = domain->get_Nx();
        size_t Ny = domain->get_Ny();

        size_t i_start = domain->get_index_x1();
        size_t i_end = domain->get_index_x2() + 1;
        size_t j_start = domain->get_index_y1();
        size_t j_end = domain->get_index_y2() + 1;
        size_t k_start = domain->get_index_z1();
        size_t k_end = domain->get_index_z2() + 1;

        //expansion - expand if there is at least one cell in the buffer area fulfills the condition
        //loop through lower side of cuboid in y direction
#pragma acc parallel loop collapse(2) present(f) reduction(+:expansion_counter_start, expansion_counter_end, reduction_counter_start, reduction_counter_end)
        for (size_t i = i_start; i < i_end; i++) {
            for (size_t k = k_start; k < k_end; k++) {
                // check innermost plane of the buffer zone on the lower side
                size_t idx_s1 = IX(i, j_start + no_buffer_cell - 1, k, Nx, Ny);
                if (std::fabs(*(f + idx_s1) - check_value) > threshold) {
                    expansion_counter_start++;
                } else {
                    size_t idx = IX(i, j_start + minimal - 1, k, Nx, Ny);
                    size_t idx2 = IX(i, j_start + minimal - 1 + no_buffer_cell, k, Nx, Ny);
                    if ((std::fabs(*(f + idx2) - check_value) > threshold ||
                         std::fabs(*(f + idx) - check_value) > threshold)) {
                        reduction_counter_start++;
                    }
                }
                size_t idx_s2 = IX(i, j_end - no_buffer_cell + 1, k, Nx, Ny);
                if (std::fabs(*(f + idx_s2) - check_value) > threshold) {
                    expansion_counter_end++;
                } else {
                    size_t idx = IX(i, j_end - minimal + 1, k, Nx, Ny);
                    size_t idx2 = IX(i, j_end - minimal + 1 - no_buffer_cell, k, Nx, Ny);
                    if ((std::fabs(*(f + idx2) - check_value) > threshold ||
                         std::fabs(*(f + idx) - check_value) > threshold)) {
                        reduction_counter_end++;
                    }
                }
            }
        }
    }
    if ((expansion_counter_start > 0 && reduction_counter_start == 0 && reduction_start) ||
        (expansion_counter_end > 0 && reduction_counter_end == 0 && reduction_end)) {
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(Adaption).name());
        m_logger->error(
            "Trying to reduce and expand at the same time (y): {}, {} | {}, {}",
            expansion_counter_start,
            reduction_counter_start,
            expansion_counter_end,
            reduction_counter_end);
#endif
        throw std::exception();
    }
    if (expansion_counter_start > 0) {
        *p_shift_x1 = -1;
    } else if (reduction_counter_start == 0 && reduce) {
        *p_shift_x1 = 1;
    }
    if (expansion_counter_end > 0) {
        *p_shift_x2 = 1;
    } else if (reduction_counter_end == 0 && reduce) {
        *p_shift_x2 = -1;
    }
    return (expansion_counter_end + expansion_counter_start) > 0 || reduction_counter_start == 0 ||
           reduction_counter_end == 0;
}
