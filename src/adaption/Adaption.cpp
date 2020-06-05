/// \file 		Adaption.h
/// \brief 		Controll class for adaption
/// \date 		Nov 29, 2018
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifdef _OPENACC
#include <accelmath.h>
#else
#include <cmath>
#endif

#include <chrono>
#include <fstream>
#include <iostream>

#include "Adaption.h"
#include "Layers.h"
#include "Vortex.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"
#include "../utility/Utility.h"

class Vortex;

class Layers;

Adaption::Adaption(Field **fields) {
    m_logger = Utility::createLogger(typeid(this).name());
    auto params = Parameters::getInstance();
    auto domain = Domain::getInstance();
    m_dynamic = (params->get("adaption/dynamic") == "Yes");
    std::string tmp = params->get("adaption/dynamic");
    m_filename = Parameters::getInstance()->get("xml_filename");
    m_filename.resize(m_filename.size()-4);//remove .xml from filename
    m_hasDataExtraction = (params->get("adaption/data_extraction") == "Yes");
    if (m_hasDataExtraction) {
        m_hasDataExtractionBefore = (params->get("adaption/data_extraction/before/enabled") == "Yes");
        m_hasDataExtractionAfter = (params->get("adaption/data_extraction/after/enabled") == "Yes");
        m_hasDataExtractionEndresult = (params->get("adaption/data_extraction/endresult/enabled") == "Yes");
        m_hasTimeMeasuring = (params->get("adaption/data_extraction/time_measuring/enabled") == "Yes");

        m_hasWriteField = (params->get("adaption/data_extraction/write_field/enabled") == "Yes");
        //m_hasWriteRuntime = (params->get("adaption/data_extraction/runtime/enabled") == "Yes");
    }
    this->fields = fields;
    if (m_dynamic) {
        std::string init = params->get("adaption/class/name");
        if (init == "Layers") {
            func = new Layers(this, fields);
        } else if (init == "Vortex" || init == "VortexY") {
            func = new Vortex(this, fields);
        } else {
            m_logger->critical("Type {} is not defined", init);
            std::exit(1);
            ///TODO Error Handling
        }

        m_reduction = func->hasReduction();
        m_dynamic_end = false;
        m_minimal = static_cast<size_t> (std::pow(2, domain->GetLevels()));
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
/// \param	t_cur	current timestep
// ***************************************************************************************
void Adaption::run(real t_cur) {

    if (m_dynamic && isUpdateNecessary()) {
        applyChanges();
#ifndef PROFILING
        std::ofstream file;
        file.open(getTimeMeasuringName(), std::ios::app);
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
#endif
        BoundaryController::getInstance()->updateLists();
#ifndef PROFILING
            end = std::chrono::system_clock::now();
        if(m_hasTimeMeasuring) {
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
/// \param	filename	filename
/// \param	height y-value
/// \param  time timestep
// ***************************************************************************************
void Adaption::extractData(const std::string filename, real height, real time) {
    auto domain = Domain::getInstance();
    size_t x_start = 0;
    size_t x_end = domain->GetNx();
    auto y = static_cast<size_t>(std::round((height - domain->GetY1()) / domain->Getdy()));
    auto z = static_cast<size_t >(std::round(domain->GetNz() / 2));
    size_t Nx = domain->GetNx();
    size_t Ny = domain->GetNy();

    real *data_temp = fields[VectorFieldsTypes::TEMPERATURE]->data;
    real *data_tempA = fields[VectorFieldsTypes::TEMPERATURE_A]->data;
    real *data_u = fields[VectorFieldsTypes::VEL_U]->data;
    real *data_v = fields[VectorFieldsTypes::VEL_V]->data;
    real *data_w = fields[VectorFieldsTypes::VEL_W]->data;
    real *data_nu = fields[VectorFieldsTypes::NU_T]->data;
    real *data_kappa = fields[VectorFieldsTypes::KAPPA_T]->data;

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
/// \param	filename	filename
// ***************************************************************************************
void Adaption::extractData(const std::string filename) {
    auto domain = Domain::getInstance();
    size_t Nx = domain->GetNx();
    size_t Ny = domain->GetNy();
    size_t x_end = Nx;
    size_t y_end = Ny;
    size_t z_end = domain->GetNz();

    real *data_temp = fields[VectorFieldsTypes::TEMPERATURE]->data;
    real *data_u = fields[VectorFieldsTypes::VEL_U]->data;
    real *data_v = fields[VectorFieldsTypes::VEL_V]->data;
    real *data_w = fields[VectorFieldsTypes::VEL_W]->data;

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
#ifndef PROFILING
  std::ofstream file;
    file.open(getTimeMeasuringName(), std::ios::app);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
#endif
    auto domain = Domain::getInstance();
    if (domain->Resize(m_shift_x1, m_shift_x2, m_shift_y1, m_shift_y2, m_shift_z1, m_shift_z2)) {
        if (domain->GetX1() == domain->Getx1() &&
            domain->GetX2() == domain->Getx2() &&
            domain->GetY1() == domain->Gety1() &&
            domain->GetY2() == domain->Gety2() &&
            domain->GetZ1() == domain->Getz1() &&
            domain->GetZ2() == domain->Getz2()) {
            m_dynamic_end = true;
        }
        func->applyChanges();
    }
#ifndef PROFILING
    end = std::chrono::system_clock::now();
    long ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    if(m_hasTimeMeasuring) {
        file << "apply: " << ms << " microsec\n";
    }
    file.close();
#endif
}

// ==================================== Is update necessary ==============================
// ***************************************************************************************
/// \brief  Checks if adaption should be done
/// \return	bool true if yes false if no
// ***************************************************************************************
bool Adaption::isUpdateNecessary() {
#ifndef PROFILING
    std::ofstream file;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    if (m_hasTimeMeasuring) {
        file.open(getTimeMeasuringName(), std::ios::app);
        start = std::chrono::system_clock::now();
    }
#endif
    bool update = false;
    if (!m_dynamic_end) {
        update = func->update();
    }
#ifndef PROFILING
    if(m_hasTimeMeasuring) {
        end = std::chrono::system_clock::now();
        long ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        file << "update: " << ms << " microsec\n";
        file.close();
    }
#endif
    return update;
}

// ==================================== Expand x direction ===============================
// ***************************************************************************************
/// \brief  Basic implementation for expansion in x direction (parallelized)
/// \param	shift expansion size value in x direction
/// \param  start indicates whether the expansion is at the beginning or the end of the computational domain
/// \param  arr_idxExpansion  Index list of cells to be newly added
/// \param len_e  size of arr_idxExpansion
// ***************************************************************************************
void Adaption::expandXDirection(long shift, bool start, size_t *arr_idxExpansion, size_t len_e) {
    auto domain = Domain::getInstance();
#pragma acc data present(arr_idxExpansion[:len_e])
    {
        size_t j_start = domain->GetIndexy1();// (y1 - Y1) / dy;
        size_t j_end = domain->GetIndexy2() + 2;//(y2 - Y1) / dy + 2;

        size_t k_start = domain->GetIndexz1();//(z1 - Z1) / dz;
        size_t k_end = domain->GetIndexz2() + 2;//(z2 - Z1) / dz + 2;

        size_t Nx = domain->GetNx();
        size_t Ny = domain->GetNy();

        if (start) { // shift = negative
            size_t i1 = domain->GetIndexx1();//(x1 - X1) / dx;
            shift *= (-1);
#pragma acc parallel loop collapse(3) present(arr_idxExpansion[:len_e])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t j = j_start; j < j_end; j++) {
                    for (long ii = 0; ii < shift; ii++) {
                        size_t idx = IX(i1 + ii, j, k, Nx, Ny);
                        size_t counter = (k - k_start) * (shift) * (j_end - j_start) + (j - j_start) * (shift) + ii;
                        *(arr_idxExpansion + counter) = idx;
                    }
                }
            }
        } else { // shift = positive
            size_t i2 = domain->GetIndexx2();//(x2 - X1) / dx + 1;
#pragma acc parallel loop collapse(3) present(arr_idxExpansion[:len_e])
            for (size_t j = j_start; j < j_end; j++) {
                for (size_t k = k_start; k < k_end; k++) {
                    for (long ii = 0; ii < shift; ii++) {
                        size_t idx = IX(i2 - ii, j, k, Nx, Ny);
                        size_t counter = (k - k_start) * (shift) * (j_end - j_start) + (j - j_start) * (shift) + ii;
                        *(arr_idxExpansion + counter) = idx;
                    }
                }
            }
        }
    }
}

// ==================================== Expand y direction ===============================
// ***************************************************************************************
/// \brief  Basic implementation for expansion in y direction (parallelized)
/// \param	shift expansion size value in y direction
/// \param  start indicates whether the expansion is at the beginning or the end of the computational domain
/// \param  arr_idxExpansion  Index list of cells to be newly added
/// \param len_e  size of arr_idxExpansion
// ***************************************************************************************
void Adaption::expandYDirection(long shift, bool start, size_t *arr_idxExpansion, size_t len_e) {
    auto domain = Domain::getInstance();

#pragma acc data present(arr_idxExpansion[:len_e])
    {
        size_t i_start = domain->GetIndexx1();//(x1 - X1) / dx;
        size_t i_end = domain->GetIndexx2() + 2;//(x2 - X1) / dx + 2;
        size_t k_start = domain->GetIndexz1();//(z1 - Z1) / dz;
        size_t k_end = domain->GetIndexz2() + 2;//(z2 - Z1) / dz + 2;
        size_t Nx = domain->GetNx();
        size_t Ny = domain->GetNy();

        if (start) { // shift = negative
            size_t j1 = domain->GetIndexy1();//(y1 - Y1) / dy;
            shift *= (-1);
#pragma acc parallel loop collapse(3) present(arr_idxExpansion[:len_e])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t i = i_start; i < i_end; i++) {
                    for (long jj = 0; jj < shift; jj++) {
                        size_t idx = IX(i, j1 + jj, k, Nx, Ny);
                        size_t counter = (k - k_start) * (shift) * (i_end - i_start) + (i - i_start) * (shift) + jj;
                        *(arr_idxExpansion + counter) = idx;
                    }
                }
            }
        } else { // shift = positive
            size_t j2 = domain->GetIndexy2() + 1;//(y2 - Y1) / dy + 1;
#pragma acc parallel loop collapse(3) present(arr_idxExpansion[:len_e])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t i = i_start; i < i_end; i++) {
                    for (long jj = 0; jj < shift; jj++) {
                        size_t idx = IX(i, j2 - jj, k, Nx, Ny);
                        size_t counter = (k - k_start) * (shift) * (i_end - i_start) + (i - i_start) * (shift) + jj;
                        *(arr_idxExpansion + counter) = idx;
                    }
                }
            }
        }
    }
}

// ==================================== Reduction x direction ===============================
// ***************************************************************************************
/// \brief  Basic implementation for expansion in x direction (parallelized)
/// \param	shift reduction size value in x direction
/// \param  start indicates whether the reduction is at the beginning or the end of the computational domain
/// \param  arr_idxReduction  Index list of cells to be newly added
/// \param len_e  size of arr_idxReduction
// ***************************************************************************************
void Adaption::reduceXDirection(long shift_inp, bool start, size_t *arr_idxReduction, size_t len_r) {

    auto domain = Domain::getInstance();
    unsigned long shift = shift_inp;
#pragma acc data present(arr_idxReduction[:len_r])
    {
        size_t j_start = domain->GetIndexy1();//(y1 - Y1) / dy;
        size_t j_end = domain->GetIndexy2() + 2;//(y2 - Y1) / dy + 2;
        size_t k_start = domain->GetIndexz1();//(z1 - Z1) / dz;
        size_t k_end = domain->GetIndexz2() + 2;//(z2 - Z1) / dz + 2;
        size_t Nx = domain->GetNx();
        size_t Ny = domain->GetNy();

        if (start) { // shift = positive
            size_t i1 = domain->GetIndexx1();//(domain->Getx1() - X1) / dx;
#pragma acc parallel loop collapse(3) present(arr_idxReduction[:len_r])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t j = j_start; j < j_end; j++) {
                    for (size_t ii = 0; ii < shift; ii++) {
                        size_t counter = (k - k_start) * (shift) * (j_end - j_start) + (j - j_start) * (shift) + ii;
                        auto idx = IX(i1 - ii - 1, j, k, Nx, Ny);
                        *(arr_idxReduction + counter) = idx;
                    }
                }
            }
        } else { // shift = negative
            shift *= (-1);
            size_t i2 = domain->GetIndexx2();//(domain->Getx2() - X1) / dx;
#pragma acc parallel loop collapse(3) present(arr_idxReduction[:len_r])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t j = j_start; j < j_end; j++) {
                    for (size_t ii = 0; ii < shift; ii++) {
                        size_t counter = (k - k_start) * (shift) * (j_end - j_start) + (j - j_start) * (shift) + ii;
                        auto idx = IX(i2 + ii + 2, j, k, Nx, Ny);
                        *(arr_idxReduction + counter) = idx;
                    }
                }
            }
        }
    }
}

// ==================================== Reduction y direction ===============================
// ***************************************************************************************
/// \brief  Basic implementation for expansion in y direction (parallelized)
/// \param	shift reduction size value in y direction
/// \param  start indicates whether the reduction is at the beginning or the end of the computational domain
/// \param  arr_idxReduction  Index list of cells to be newly added
/// \param len_e  size of arr_idxReduction
// ***************************************************************************************
void Adaption::reduceYDirection(long shift_inp, bool start, size_t *arr_idxReduction, size_t len_r) {
    // std::cout << "reduceYDirection" << std::endl;

    auto domain = Domain::getInstance();
    unsigned long shift = shift_inp;
#pragma acc data present(arr_idxReduction[:len_r])
    {
        size_t i_start = domain->GetIndexx1();//(x1 - X1) / dx;
        size_t i_end = domain->GetIndexx2() + 2;//(x2 - X1) / dx + 2;
        size_t k_start = domain->GetIndexz1();//(z1 - Z1) / dz;
        size_t k_end = domain->GetIndexz2() + 2;//(z2 - Z1) / dz + 2;

        size_t Nx = domain->GetNx();
        size_t Ny = domain->GetNy();

        if (start) { // shift = positive
            size_t j1 = domain->GetIndexy1();//(domain->Gety1() - Y1) / dy;
#pragma acc parallel loop collapse(3) present(arr_idxReduction[:len_r])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t i = i_start; i < i_end; i++) {
                    for (size_t jj = 0; jj < shift; jj++) {
                        size_t counter = (k - k_start) * (shift) * (i_end - i_start) + (i - i_start) * (shift) + jj;
                        auto idx = IX(i, j1 - jj - 1, k, Nx, Ny);
                        *(arr_idxReduction + counter) = idx;
                    }
                }
            }
        } else { // shift = negative
            shift *= (-1);
            size_t j2 = domain->GetIndexy2();//(domain->Gety2() - Y1) / dy;
#pragma acc parallel loop collapse(3) present(arr_idxReduction[:len_r])
            for (size_t k = k_start; k < k_end; k++) {
                for (size_t i = i_start; i < i_end; i++) {
                    for (size_t jj = 0; jj < shift; jj++) {
                        size_t counter = (k - k_start) * (shift) * (i_end - i_start) + (i - i_start) * (shift) + jj;
                        auto idx = IX(i, j2 + jj + 2, k, Nx, Ny);
                        *(arr_idxReduction + counter) = idx;
                    }
                }
            }
        }
    }
}

// ==================================== Adaption x direction serial ===============================
// ***************************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param	f field
/// \param  checkValue check value
/// \param  noBufferCell Buffersize
/// \param  threshold precision of comparison
// ***************************************************************************************
bool Adaption::adaptXDirection_serial(const real *f, real checkValue, size_t noBufferCell, real threshold) {
    auto domain = Domain::getInstance();
    size_t Nx = domain->GetNx();
    size_t Ny = domain->GetNy();

    size_t nx_begin = domain->GetIndexx1();// ((domain->Getx1() - domain->GetX1()) / domain->Getdx());

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

    //reduction - reduce if all cells do not fulfil the condition any longer
    ADTypes reduction_start = ADTypes::NO;
    ADTypes reduction_end = ADTypes::NO;
    ADTypes reduction = ADTypes::NO;
    if (m_reduction && domain->Getnx() > (m_minimal + noBufferCell) * 2) {
        reduction_start = ADTypes::UNKNOWN;
        reduction_end = ADTypes::UNKNOWN;
        reduction = ADTypes::UNKNOWN;
    }
    bool check = !(expansion == ADTypes::YES || (reduction == ADTypes::NO && expansion == ADTypes::NO));
    //loop through left and right side of cuboid in x direction
    for (size_t j = j_start; j < j_end && check; j++) {
        for (size_t k = k_start; k < k_end && check; k++) {
            // check innermost plane of the buffer zone on the left side
            size_t idx_s1 = IX(i_start + noBufferCell - 1, j, k, Nx, Ny);
            if (expansion_start == ADTypes::UNKNOWN && std::fabs(*(f + idx_s1) - checkValue) > threshold) {
                m_shift_x1 = -1;
                expansion_start = ADTypes::YES;
            } else {
                // check innermost plane of the minimal zone to reduce on the left side
                size_t idx = IX(i_start + m_minimal - 1, j, k, Nx, Ny);
                size_t idx2 = IX(i_start + m_minimal - 1 + noBufferCell, j, k, Nx, Ny);
                if (reduction_start == ADTypes::UNKNOWN && (std::fabs(*(f + idx2) - checkValue) > threshold ||
                                                            std::fabs(*(f + idx) - checkValue) > threshold)) {
                    reduction_start = ADTypes::NO;
                }
            }
            // check innermost plane of the buffer zone on the right side
            size_t idx_s2 = IX(i_end - noBufferCell + 1, j, k, Nx, Ny);
            if (expansion_end == ADTypes::UNKNOWN && std::fabs(*(f + idx_s2) - checkValue) > threshold) {
                m_shift_x2 = 1;
                expansion_end = ADTypes::YES;
            } else {
                // check innermost plane of the minimal zone to reduce on the right side
                size_t idx = IX(i_end - m_minimal + 1, j, k, Nx, Ny);
                size_t idx2 = IX(i_end - m_minimal + 1 - noBufferCell, j, k, Nx, Ny);
                if (reduction_end == ADTypes::UNKNOWN && (std::fabs(*(f + idx2) - checkValue) > threshold ||
                                                          std::fabs(*(f + idx) - checkValue) > threshold)) {
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
        m_logger->error("Exception in x-Adaption: {} {} {} {}", size_t(expansion_start),
                                                              size_t(reduction_start),
                                                              size_t(expansion_end),
                                                              size_t(reduction_end));
        //TODO Error handling
        //throw std::exception();
    }
    if (reduction_start == ADTypes::UNKNOWN && expansion_start != ADTypes::YES) {
        reduction_start = ADTypes::YES;
        m_shift_x1 = 1;
    }
    if (reduction_end == ADTypes::UNKNOWN && expansion_end != ADTypes::YES) {
        reduction_end = ADTypes::YES;
        m_shift_x2 = -1;
    }
    return expansion_end == ADTypes::YES || expansion_start == ADTypes::YES || reduction_start == ADTypes::YES ||
           reduction_end == ADTypes::YES;
}

// ==================================== Adaption x direction parallel ===============================
// ***************************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param	f field
/// \param  checkValue check value
/// \param  noBufferCell Buffersize
/// \param  threshold precision of comparison
// ***************************************************************************************
bool Adaption::adaptXDirection(const real *f, real checkValue, size_t noBufferCell, real threshold) {
    auto domain = Domain::getInstance();
    size_t expansion_counter_start = 0;
    size_t expansion_counter_end = 0;
    size_t reduction_counter_start = 0;
    size_t reduction_counter_end = 0;

    bool reduction_start = m_reduction;
    bool reduction_end = m_reduction;
    if (m_reduction && domain->Getnx() > (m_minimal + noBufferCell) * 2) {
        reduction_start = false;
        reduction_end = false;
    }
#pragma acc data present(f) copy(expansion_counter_start, expansion_counter_end, reduction_counter_start, reduction_counter_end)
    {
        size_t Nx = domain->GetNx();
        size_t Ny = domain->GetNy();

        size_t nx_begin = domain->GetIndexx1();// ((domain->Getx1() - domain->GetX1()) / domain->Getdx());

        size_t j_start = domain->GetIndexy1();//((domain->Gety1() - domain->GetY1()) / domain->Getdy());
        size_t j_end = j_start + domain->Getny() - 1;
        size_t k_start = domain->GetIndexz1();//((domain->Getz1() - domain->GetZ1()) / domain->Getdz());
        size_t k_end = k_start + domain->Getnz() - 1;

        size_t i_start = nx_begin;
        size_t i_end = nx_begin + domain->Getnx() - 1;

        size_t minimal = m_minimal;

        //expansion - expand if there is at least one cell in the buffer area fulfills the condition
        //reduction - reduce if all cells do not fulfil the condition any longer
        //loop through left side of cuboid in x direction
#pragma acc parallel loop collapse(2) present(f) reduction(+:expansion_counter_start, reduction_counter_start, expansion_counter_end, reduction_counter_end)
        for (size_t j = j_start; j < j_end; j++) {
            for (size_t k = k_start; k < k_end; k++) {
                // check innermost plane of the buffer zone on the left side
                size_t idx_s1 = IX(i_start + noBufferCell - 1, j, k, Nx, Ny);
                if (std::fabs(*(f + idx_s1) - checkValue) > threshold) {
                    expansion_counter_start++;
                } else {
                    size_t idx = IX(i_start + minimal - 1, j, k, Nx, Ny);
                    size_t idx2 = IX(i_start + minimal - 1 + noBufferCell, j, k, Nx, Ny);
                    if ((std::fabs(*(f + idx2) - checkValue) > threshold ||
                         std::fabs(*(f + idx) - checkValue) > threshold)) {
                        reduction_counter_start++;
                    }
                }
                size_t idx_s2 = IX(i_end - noBufferCell + 1, j, k, Nx, Ny);
                if (std::fabs(*(f + idx_s2) - checkValue) > threshold) {
                    expansion_counter_end++;
                } else {
                    size_t idx = IX(i_end - minimal + 1, j, k, Nx, Ny);
                    size_t idx2 = IX(i_end - minimal + 1 - noBufferCell, j, k, Nx, Ny);
                    if ((std::fabs(*(f + idx2) - checkValue) > threshold ||
                         std::fabs(*(f + idx) - checkValue) > threshold)) {
                        reduction_counter_end++;
                    }
                }
            }
        }
    }
    if ((expansion_counter_start > 0 && reduction_counter_start == 0 && reduction_start) ||
        (expansion_counter_end > 0 && reduction_counter_end == 0 && reduction_end)) {
        m_logger->error("Trying to reduce and expand at the same time (x): {},{} | {},{}",
                expansion_counter_start,
                reduction_counter_start,
                expansion_counter_end,
                reduction_counter_end);
        //TODO Error Handling
        //throw std::exception();
    }
    if (expansion_counter_start > 0) {
        m_shift_x1 = -1;
    } else {
        if (reduction_counter_start == 0 && m_reduction) {
            m_shift_x1 = 1;
        }
    }
    if (expansion_counter_end > 0) {
        m_shift_x2 = 1;
    } else {
        if (reduction_counter_end == 0 && m_reduction) {
            m_shift_x2 = -1;
        }
    }
    return (expansion_counter_end + expansion_counter_start) > 0 || reduction_counter_start == 0 ||
           reduction_counter_end == 0;
}

// ==================================== Adaption y direction serial ===============================
// ***************************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param	f field
/// \param  checkValue check value
/// \param  noBufferCell Buffersize
/// \param  threshold precision of comparison
// ***************************************************************************************
bool Adaption::adaptYDirection_serial(const real *f, real checkValue, size_t noBufferCell, real threshold) {
    auto domain = Domain::getInstance();
    size_t Nx = domain->GetNx();
    size_t Ny = domain->GetNy();


    size_t i_start = domain->GetIndexx1();//((domain->Getx1() - domain->GetX1()) / domain->Getdx());
    size_t i_end = i_start + domain->Getnx() - 1;
    size_t k_start = domain->GetIndexz1();//((domain->Getz1() - domain->GetZ1()) / domain->Getdz());
    size_t k_end = k_start + domain->Getnz() - 1;

    size_t ny_begin = domain->GetIndexy1();//((domain->Gety1() - domain->GetY1()) / domain->Getdy());
    size_t j_start = ny_begin;
    size_t j_end = ny_begin + domain->Getny() - 1;

    //expansion - expand if there is at least one cell in the buffer area fulfills the condition
    ADTypes expansion_start = ADTypes::UNKNOWN;
    if (domain->Gety1() == domain->GetY1()) {
        expansion_start = ADTypes::NO;
    }
    ADTypes expansion_end = ADTypes::UNKNOWN;
    if (domain->Gety2() == domain->GetY2()) {
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
    if (m_reduction && domain->Getnx() > (m_minimal + noBufferCell) * 2) {
        reduction_start = ADTypes::UNKNOWN;
        reduction_end = ADTypes::UNKNOWN;
        reduction = ADTypes::UNKNOWN;
    }
    bool check = !(expansion == ADTypes::YES || (reduction == ADTypes::NO && expansion == ADTypes::NO));
    //loop through left and right side of cuboid in x direction
    for (size_t i = i_start; i < i_end && check; i++) {
        for (size_t k = k_start; k < k_end && check; k++) {
            // check innermost plane of the buffer zone on the lower side
            size_t idx_s1 = IX(i, j_start + noBufferCell - 1, k, Nx, Ny);
            if (expansion_start == ADTypes::UNKNOWN && std::fabs(*(f + idx_s1) - checkValue) > threshold) {
                m_shift_y1 = -1;
                expansion_start = ADTypes::YES;
            } else {
                // check innermost plane of the minimal zone to reduce on the lower side
                size_t idx = IX(i, j_start + m_minimal - 1, k, Nx, Ny);
                size_t idx2 = IX(i, j_start + m_minimal - 1 + noBufferCell, k, Nx, Ny);
                if (reduction_start == ADTypes::UNKNOWN && (std::fabs(*(f + idx2) - checkValue) > threshold ||
                                                            std::fabs(*(f + idx) - checkValue) > threshold)) {
                    reduction_start = ADTypes::NO;
                }
            }
            // check innermost plane of the buffer zone on the upper side
            size_t idx_s2 = IX(i, j_end - noBufferCell + 1, k, Nx, Ny);
            if (expansion_end == ADTypes::UNKNOWN && std::fabs(*(f + idx_s2) - checkValue) > threshold) {
                m_shift_y2 = 1;
                expansion_end = ADTypes::YES;
            } else {
                // check innermost plane of the minimal zone to reduce on the upper side
                size_t idx = IX(i, j_end - m_minimal + 1, k, Nx, Ny);
                size_t idx2 = IX(i, j_end - m_minimal + 1 - noBufferCell, k, Nx, Ny);
                if (reduction_end == ADTypes::UNKNOWN && (std::fabs(*(f + idx2) - checkValue) > threshold ||
                                                          std::fabs(*(f + idx) - checkValue) > threshold)) {
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
        m_logger->error("Exception in y-Adaption: {} {} {} {}", size_t(expansion_start),
                                                              size_t(reduction_start),
                                                              size_t(expansion_end),
                                                              size_t(reduction_end));
        //TODO Error handling
        //throw std::exception();
    }
    if (reduction_start == ADTypes::UNKNOWN && expansion_start != ADTypes::YES) {
        reduction_start = ADTypes::YES;
        m_shift_y1 = 1;
    }
    if (reduction_end == ADTypes::UNKNOWN && expansion_end != ADTypes::YES) {
        reduction_end = ADTypes::YES;
        m_shift_y2 = -1;
    }
    return expansion_end == ADTypes::YES || expansion_start == ADTypes::YES ||
           reduction_start == ADTypes::YES ||
           reduction_end == ADTypes::YES;
}

// ==================================== Adaption x direction parallel ===============================
// ***************************************************************************************
/// \brief  Checks if adaption is possible and allowed
/// \param	f field
/// \param  checkValue check value
/// \param  noBufferCell Buffersize
/// \param  threshold precision of comparison
// ***************************************************************************************
bool Adaption::adaptYDirection(const real *f, real checkValue, size_t noBufferCell, real threshold) {
    auto domain = Domain::getInstance();

    size_t expansion_counter_start = 0;
    size_t expansion_counter_end = 0;
    size_t reduction_counter_start = 0;
    size_t reduction_counter_end = 0;
    //reduction - reduce if all cells do not fulfil the condition any longer
    bool reduction_start = m_reduction;
    bool reduction_end = m_reduction;
    if (m_reduction && domain->Getny() > (m_minimal + noBufferCell) * 2) {
        reduction_start = false;
        reduction_end = false;
    }
#pragma acc data present(f) copy(expansion_counter_start, expansion_counter_end, reduction_counter_start, reduction_counter_end)
    {
        size_t Nx = domain->GetNx();
        size_t Ny = domain->GetNy();

        size_t ny_begin = domain->GetIndexy1();//((domain->Gety1() - domain->GetY1()) / domain->Getdy());

        size_t i_start = domain->GetIndexx1();//((domain->Getx1() - domain->GetX1()) / domain->Getdx());
        size_t i_end = i_start + domain->Getnx() - 1;
        size_t k_start = domain->GetIndexz1();//((domain->Getz1() - domain->GetZ1()) / domain->Getdz());
        size_t k_end = k_start + domain->Getnz() - 1;

        size_t j_start = ny_begin;
        size_t j_end = ny_begin + domain->Getny() - 1;

        size_t minimal = m_minimal;

        //expansion - expand if there is at least one cell in the buffer area fulfills the condition
        //loop through lower side of cuboid in y direction
#pragma acc parallel loop collapse(2) present(f) reduction(+:expansion_counter_start, expansion_counter_end, reduction_counter_start, reduction_counter_end)
        for (size_t i = i_start; i < i_end; i++) {
            for (size_t k = k_start; k < k_end; k++) {
                // check innermost plane of the buffer zone on the lower side
                size_t idx_s1 = IX(i, j_start + noBufferCell - 1, k, Nx, Ny);
                if (std::fabs(*(f + idx_s1) - checkValue) > threshold) {
                    expansion_counter_start++;
                } else {
                    size_t idx = IX(i, j_start + minimal - 1, k, Nx, Ny);
                    size_t idx2 = IX(i, j_start + minimal - 1 + noBufferCell, k, Nx, Ny);
                    if ((std::fabs(*(f + idx2) - checkValue) > threshold ||
                         std::fabs(*(f + idx) - checkValue) > threshold)) {
                        reduction_counter_start++;
                    }
                }
                size_t idx_s2 = IX(i, j_end - noBufferCell + 1, k, Nx, Ny);
                if (std::fabs(*(f + idx_s2) - checkValue) > threshold) {
                    expansion_counter_end++;
                } else {
                    size_t idx = IX(i, j_end - minimal + 1, k, Nx, Ny);
                    size_t idx2 = IX(i, j_end - minimal + 1 - noBufferCell, k, Nx, Ny);
                    if ((std::fabs(*(f + idx2) - checkValue) > threshold ||
                         std::fabs(*(f + idx) - checkValue) > threshold)) {
                        reduction_counter_end++;
                    }
                }
            }
        }
    }
    if ((expansion_counter_start > 0 && reduction_counter_start == 0 && reduction_start) ||
        (expansion_counter_end > 0 && reduction_counter_end == 0 && reduction_end)) {
        m_logger->error("Trying to reduce and expand at the same time (y): {}, {} | {}, {}",
            expansion_counter_start,
            reduction_counter_start,
            expansion_counter_end,
            reduction_counter_end);
        //TODO Error handling
        //throw std::exception();
    }
    if (expansion_counter_start > 0) {
        m_shift_y1 = -1;
    } else if (reduction_counter_start == 0 && m_reduction) {
        m_shift_y1 = 1;
    }
    if (expansion_counter_end > 0) {
        m_shift_y2 = 1;
    } else if (reduction_counter_end == 0 && m_reduction) {
        m_shift_y2 = -1;
    }
    return (expansion_counter_end + expansion_counter_start) > 0 || reduction_counter_start == 0 ||
           reduction_counter_end == 0;
}
