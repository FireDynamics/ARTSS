/// \file       JoinedList.cpp
/// \brief      
/// \date       Oct 22, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include "SingleJoinedList.h"

SingleJoinedList::SingleJoinedList(Settings::Settings const &settings, size_t multigrid_level) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(settings, typeid(this).name());
#endif
#ifdef GPU_DEBUG
    m_gpu_logger = Utility::create_gpu_logger(typeid(this).name());
#endif
    m_index_list = new size_t[multigrid_level + 2];
    std::fill(m_index_list, m_index_list + multigrid_level + 2, 0);

    m_size_list = new size_t[multigrid_level + 1];
    std::fill(m_size_list, m_size_list + multigrid_level + 1, 0);
}

void SingleJoinedList::set_size(const size_t size) {
    m_size = size;
    m_data = new size_t[m_size];
}

size_t SingleJoinedList::get_slice_size(size_t level) const {
    return m_size_list[level];
}

size_t SingleJoinedList::get_first_index(size_t level) const {
    return m_index_list[level];
}

/**
 * for data: get last index of slice of level
 * @param level  level of slice
 * @return
 */
size_t SingleJoinedList::get_last_index(size_t level) const {
    size_t size = get_slice_size(level);
    if (size == 0) {
        return get_first_index(level);
    }
    return get_first_index(level) + size - 1;
}

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
    delete[] m_data;
    delete[] m_size_list;
}

size_t *SingleJoinedList::get_slice(size_t level) {
    size_t slice_size = get_slice_size(level);
    size_t *slice = new size_t[slice_size];
    if (slice_size > 0) {
        //TODO(cvm) legit? want to return a copy of a part of 'data'
        std::copy(m_data + get_first_index(level), m_data + get_last_index(level) + 1, slice);
    }
    return slice;
}
