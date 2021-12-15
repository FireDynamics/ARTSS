/// \file       SingleJoinedList.cpp
/// \brief      stores indices lists for GPU, only differs between multilevel for multigrid.
/// \details    Meant to be for Domain object, as there is only one for each level
/// \date       Oct 22, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include "SingleJoinedList.h"

SingleJoinedList::SingleJoinedList(size_t multigrid_level) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_index_list = new size_t[multigrid_level + 2];
    std::fill(m_index_list, m_index_list + multigrid_level + 2, 0);

    m_size_list = new size_t[multigrid_level + 1];
    std::fill(m_size_list, m_size_list + multigrid_level + 1, 0);
}

/// set size of GPU array
/// @param size size of GPU array
void SingleJoinedList::set_size(const size_t size) {
    m_size = size;
    m_data = new size_t[m_size];
}

/// get size of GPU array for a specific level
/// @param level multigrid level
/// @return size of GPU array for a specific level
size_t SingleJoinedList::get_slice_size(size_t level) const {
    return m_size_list[level];
}

/// get starting index of GPU array for a specific level
/// @param level multigrid level
/// @return starting index of GPU array for the specified level
size_t SingleJoinedList::get_first_index(size_t level) const {
    return m_index_list[level];
}

/// get ending index of GPU array for a specific level. if size of level equals 0
/// the first index will be returned instead. Only possible if there aren't any
/// objects, which always results in a return of 0.
/// @param level multigrid level
/// @return ending index of GPU array for the specified level
size_t SingleJoinedList::get_last_index(size_t level) const {
    size_t size = get_slice_size(level);
    if (size == 0) {
        return get_first_index(level);
    }
    return get_first_index(level) + size - 1;
}

/// add data (index list) for one multigrid level. data will be copied
/// @param data index list to be added
/// @param level multigrid level of index list
/// @param size size of data
void SingleJoinedList::add_data(size_t level, size_t size, const size_t *data) {
    size_t index = m_index_list[level];
    m_index_list[level + 1] = index + size;
    m_size_list[level] = size;
    for (size_t i = 0; i < size; i++) {
        m_data[index + i] = data[i];
    }
#ifndef BENCHMARKING
    m_logger->debug("SJL: add data for level {} with size {}|{} first index: {} last index: {}",
                    level, size, m_size_list[level], m_data[index], m_data[index + size - 1]);
#endif
}

SingleJoinedList::~SingleJoinedList() {
#pragma acc exit data delete(m_data[:m_size])
    delete[] m_index_list;
    if (m_size > 0) {
        delete[] m_data;
    }
    delete[] m_size_list;
}

/// get a copy of the index list of a specific level
/// @param level multigrid level
/// @return copy of a part of index list. Transfers the ownership to the caller.
///         Caller needs to free the memory
size_t *SingleJoinedList::get_slice(size_t level) {
    size_t slice_size = get_slice_size(level);
    auto slice = new size_t[slice_size];
    if (slice_size > 0) {
        //TODO(cvm) legit? want to return a copy of a part of 'data'
        std::copy(m_data + get_first_index(level), m_data + get_last_index(level) + 1, slice);
    }
    return slice;
}
