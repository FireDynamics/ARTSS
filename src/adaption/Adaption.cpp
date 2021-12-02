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

Adaption::Adaption(Settings::Settings const &settings, FieldController *field_controller) :
        m_shift_start(), m_shift_end(), m_field_controller(field_controller) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(settings, typeid(this).name());
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
        m_shift_start.set_coordinate(0, 0, 0);
        m_shift_end.set_coordinate(0, 0, 0);
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
/// \brief  Applies domain_data adaption
// ***************************************************************************************
void Adaption::applyChanges() {
#ifndef BENCHMARKING
    std::ofstream file;
    file.open(get_time_measuring_name(), std::ios::app);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
#endif
    auto domain_data = DomainData::getInstance();
    if (domain_data->resize(m_shift_start, m_shift_end)) {
        if (domain_data->get_X1() == domain_data->get_x1() &&
            domain_data->get_X2() == domain_data->get_x2() &&
            domain_data->get_Y1() == domain_data->get_y1() &&
            domain_data->get_Y2() == domain_data->get_y2() &&
            domain_data->get_Z1() == domain_data->get_z1() &&
            domain_data->get_Z2() == domain_data->get_z2()) {
            m_dynamic_end = true;
        }
        func->apply_changes(&m_shift_start, &m_shift_end);
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
        update = func->update(&m_shift_start, &m_shift_end);
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
/// \param  shift_value expansion size value in x direction
/// \param  start indicates whether the expansion is at the beginning or the end of the computational domain_data
/// \param  arr_idx_expansion  Index list of cells to be newly added
/// \param len_e  size of arr_idx_expansion
// ***************************************************************************************
void Adaption::expand(Coordinate<long> *shift, bool start, size_t *arr_idx_expansion, size_t len_e, CoordinateAxis axis) {
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
    auto tmp = new Coordinate<size_t>();
#pragma acc data present(arr_idx_expansion[:len_e]) copyin(shift, tmp, other_axes)
    {
        long shift_value = (*shift)[axis];
        size_t other_axes0_start = domain_data->get_start_index_CD(other_axes[0]);
        size_t other_axes0_end = domain_data->get_end_index_CD(other_axes[0]);
        size_t other_axes1_start = domain_data->get_start_index_CD(other_axes[1]);
        size_t other_axes1_end = domain_data->get_end_index_CD(other_axes[1]);

        size_t Nx = domain_data->get_Nx();
        size_t Ny = domain_data->get_Ny();

        size_t counter = 0;
        if (start) { // shift = negative
            size_t i1 = domain_data->get_start_index_CD(axis);
            shift_value *= (-1);
#pragma acc parallel loop collapse(3) present(arr_idx_expansion[:len_e])
            for (size_t k = other_axes1_start; k <= other_axes1_end; k++) {
                for (size_t j = other_axes0_start; j <= other_axes0_end; j++) {
                    for (long ii = 0; ii < shift_value; ii++) {
                        (*tmp)[other_axes[1]] = k;
                        (*tmp)[other_axes[0]] = j;
                        (*tmp)[axis] = i1 + ii;
                        size_t idx = tmp->get_index(Nx, Ny);
                        *(arr_idx_expansion + counter++) = idx;
                        delete tmp;
                    }
                }
            }
        } else { // shift = positive
            size_t i2 = domain_data->get_end_index_CD(axis);
#pragma acc parallel loop collapse(3) present(arr_idx_expansion[:len_e])
            for (size_t j = other_axes0_start; j <= other_axes0_end; j++) {
                for (size_t k = other_axes1_start; k <= other_axes1_end; k++) {
                    for (long ii = 0; ii < shift_value; ii++) {
                        auto tmp = new Coordinate<size_t>();
                        (*tmp)[other_axes[1]] = k;
                        (*tmp)[other_axes[0]] = j;
                        (*tmp)[axis] = i2 - ii;
                        size_t idx = tmp->get_index(Nx, Ny);
                        *(arr_idx_expansion + counter++) = idx;
                        delete tmp;
                    }
                }
            }
        }
    }
    delete[] other_axes;
}

// ==================================== Reduction x direction ===============================
// ***************************************************************************************
/// \brief  Basic implementation for expansion in x direction (parallelized)
/// \param  shift_value reduction size value in x direction
/// \param  start indicates whether the reduction is at the beginning or the end of the computational domain_data
/// \param  arr_idx_reduction  Index list of cells to be newly added
/// \param len_e  size of arr_idx_reduction
// ***************************************************************************************
void Adaption::reduce(Coordinate<long> *shift, bool start, size_t *arr_idx_reduction, size_t len_r, CoordinateAxis axis) {
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
    auto tmp = new Coordinate<size_t>();
#pragma acc data present(arr_idx_reduction[:len_r]) copyin(other_axes, tmp, shift)
    {
        long shift_value = (*shift)[axis];
        size_t other_axes0_start = domain_data->get_start_index_CD(other_axes[0]);
        size_t other_axes0_end = domain_data->get_end_index_CD(other_axes[0]);
        size_t other_axes1_start = domain_data->get_start_index_CD(other_axes[1]);
        size_t other_axes1_end = domain_data->get_end_index_CD(other_axes[1]);
        size_t Nx = domain_data->get_Nx();
        size_t Ny = domain_data->get_Ny();

        size_t counter = 0;
        if (start) { // shift = positive
            size_t i1 = domain_data->get_start_index_CD(axis);
#pragma acc parallel loop collapse(3) present(arr_idx_reduction[:len_r])
            for (size_t k = other_axes1_start; k <= other_axes1_end; k++) {
                for (size_t j = other_axes0_start; j <= other_axes0_end; j++) {
                    for (long ii = 0; ii < shift_value; ii++) {
                        (*tmp)[other_axes[1]] = k;
                        (*tmp)[other_axes[0]] = j;
                        (*tmp)[axis] = i1 - ii - 1;
                        size_t idx = tmp->get_index(Nx, Ny);
                        *(arr_idx_reduction + counter++) = idx;
                    }
                }
            }
        } else { // shift = negative
            shift_value *= (-1);
            size_t i2 = domain_data->get_end_index_CD(axis);
#pragma acc parallel loop collapse(3) present(arr_idx_reduction[:len_r])
            for (size_t k = other_axes1_start; k <= other_axes1_end; k++) {
                for (size_t j = other_axes0_start; j <= other_axes0_end; j++) {
                    for (long ii = 0; ii < shift_value; ii++) {
                        (*tmp)[other_axes[1]] = k;
                        (*tmp)[other_axes[0]] = j;
                        (*tmp)[axis] = i2 + ii + 2;
                        size_t idx = tmp->get_index(Nx, Ny);
                        *(arr_idx_reduction + counter++) = idx;
                    }
                }
            }
        }
    }
    delete tmp;
    delete[] other_axes;
}

// ========================= Adaption x direction parallel ====================
// *****************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param  field field
/// \param  check_value check value
/// \param  no_buffer_cell Buffersize
/// \param  threshold precision of comparison
// ***************************************************************************************
bool Adaption::adapt(Settings::Settings const &settings, const Field &field,
                     real check_value,
                     size_t no_buffer_cell,
                     real threshold,
                     Coordinate<long> *shift_start, Coordinate<long> *shift_end,
                     size_t minimal,
                     bool reduce,
                     CoordinateAxis axis) {
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

    size_t expansion_counter_start = 0;
    size_t expansion_counter_end = 0;
    size_t reduction_counter_start = 0;
    size_t reduction_counter_end = 0;

    bool reduction_start = reduce;
    bool reduction_end = reduce;
    if (reduce && domain_data->get_number_of_inner_cells(axis) > (minimal + no_buffer_cell) * 2) {
        reduction_start = false;
        reduction_end = false;
    }

    auto start = new Coordinate<size_t>(
            domain_data->get_start_index_CD(X),
            domain_data->get_start_index_CD(Y),
            domain_data->get_start_index_CD(Z));
    auto end = new Coordinate<size_t>(
            domain_data->get_end_index_CD(X),
            domain_data->get_end_index_CD(Y),
            domain_data->get_end_index_CD(Z));
    auto tmp = new Coordinate<size_t>();
#pragma acc data present(field) copy(expansion_counter_start, expansion_counter_end, reduction_counter_start, reduction_counter_end) copyin(shift_start, shift_end, tmp, other_axes)
    {
        size_t Nx = domain_data->get_Nx();
        size_t Ny = domain_data->get_Ny();

        // comments refer to the x-direction to make it easier to understand
        // expansion - expand if there is at least one cell in the buffer area fulfills the condition
        // reduction - reduce if all cells do not fulfil the condition any longer
        // loop through left side of cuboid in x direction
#pragma acc parallel loop collapse(2) present(field) reduction(+:expansion_counter_start, reduction_counter_start, expansion_counter_end, reduction_counter_end)
        for (size_t j = (*start)[other_axes[0]]; j <= (*end)[other_axes[0]]; j++) {
            for (size_t k = (*start)[other_axes[1]]; k <= (*end)[other_axes[1]]; k++) {
                // check innermost plane of the buffer zone on the left side
                (*tmp)[other_axes[1]] = k;
                (*tmp)[other_axes[0]] = j;
                (*tmp)[axis] = (*start)[axis] + no_buffer_cell - 1;

                size_t idx_s1 = tmp->get_index(Nx, Ny);
                if (std::fabs(field[idx_s1] - check_value) > threshold) {
                    expansion_counter_start++;
                } else {
                    (*tmp)[axis] = (*start)[axis] + minimal - 1;
                    size_t idx = tmp->get_index(Nx, Ny);
                    (*tmp)[axis] += no_buffer_cell;
                    size_t idx2 = tmp->get_index(Nx, Ny);
                    if ((std::fabs(field[idx2] - check_value) > threshold ||
                         std::fabs(field[idx] - check_value) > threshold)) {
                        reduction_counter_start++;
                    }
                }

                (*tmp)[axis] = (*end)[axis] - no_buffer_cell + 2;
                size_t idx_s2 = tmp->get_index(Nx, Ny);
                if (std::fabs(field[idx_s2] - check_value) > threshold) {
                    expansion_counter_end++;
                } else {
                    (*tmp)[axis] = (*end)[axis] - minimal + 2;
                    size_t idx = tmp->get_index(Nx, Ny);
                    (*tmp)[axis] -= no_buffer_cell;
                    size_t idx2 = tmp->get_index(Nx, Ny);
                    if ((std::fabs(field[idx2] - check_value) > threshold ||
                         std::fabs(field[idx] - check_value) > threshold)) {
                        reduction_counter_end++;
                    }
                }
            }
        }
    }
    if ((expansion_counter_start > 0 && reduction_counter_start == 0 && reduction_start) ||
        (expansion_counter_end > 0 && reduction_counter_end == 0 && reduction_end)) {
#ifndef BENCHMARKING
        auto logger = Utility::create_logger(settings, typeid(Adaption).name());
        logger->error("Trying to reduce and expand at the same time (x): {},{} | {},{}",
                      expansion_counter_start,
                      reduction_counter_start,
                      expansion_counter_end,
                      reduction_counter_end);
#endif
        //TODO Error Handling
        //throw std::exception();
    }
    if (expansion_counter_start > 0) {
        (*shift_start)[axis] = -1;
    } else {
        if (reduction_counter_start == 0 && reduce) {
            (*shift_start)[axis] = 1;
        }
    }
    if (expansion_counter_end > 0) {
        (*shift_end)[axis] = 1;
    } else {
        if (reduction_counter_end == 0 && reduce) {
            (*shift_end)[axis] = -1;
        }
    }
    delete tmp;
    delete[] other_axes;
    return (expansion_counter_end + expansion_counter_start) > 0 || reduction_counter_start == 0 ||
           reduction_counter_end == 0;
}
