/// \file       JoinedListObstacle.cpp
/// \brief      
/// \date       Oct 23, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include "MultipleJoinedList.h"

MultipleJoinedList::MultipleJoinedList(Settings::Settings const &settings, size_t multigrid_level, size_t number_of_objects) : m_number_of_objects(number_of_objects) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(settings, typeid(this).name());
#endif
#ifdef GPU_DEBUG
    m_gpu_logger = Utility::create_gpu_logger(settings, typeid(this).name());
#endif
    m_index_list = new size_t[(multigrid_level + 1) * number_of_objects + 1];
    std::fill(m_index_list, m_index_list + (multigrid_level + 1) * number_of_objects + 1, 0);

    m_size_list = new size_t[(multigrid_level + 1) * number_of_objects];
    std::fill(m_size_list, m_size_list + (multigrid_level + 1) * number_of_objects, 0);
}

void MultipleJoinedList::set_size(const size_t size) {
    m_size = size;
    m_data = new size_t[m_size];
}

size_t MultipleJoinedList::get_slice_size(size_t level) const {
    if (m_number_of_objects > 0) {
        return get_last_index(level, m_number_of_objects - 1) - get_first_index(level, 0) + 1;
    }
    return 0;
}

size_t MultipleJoinedList::get_slice_size(size_t level, size_t obstacle_id) const {
    return m_size_list[level * m_number_of_objects + obstacle_id];
}

size_t MultipleJoinedList::get_first_index(size_t level, size_t obstacle_id) const {
    return m_index_list[level * m_number_of_objects + obstacle_id];
}

/**
 * for data: get last index of slice of level
 * @param level  level of slice
 * @return
 */
size_t MultipleJoinedList::get_last_index(size_t level, size_t obstacle_id) const {
    size_t size = get_slice_size(level, obstacle_id);
    if (size == 0) {
        return get_first_index(level, obstacle_id);
    }
    return get_first_index(level, obstacle_id) + size - 1;
}

void MultipleJoinedList::add_data(size_t level, size_t obstacle_id, size_t size, const size_t *data) {
#ifndef BENCHMARKING
    m_logger->debug("MJL: add data for obstacle id={} with level {} and size {}", obstacle_id, level, size);
#endif
    size_t index = m_index_list[level * m_number_of_objects + obstacle_id];
    m_index_list[level * m_number_of_objects + obstacle_id + 1] = index + size;
    m_size_list[level * m_number_of_objects + obstacle_id] = size;
    for (size_t i = 0; i < size; i++) {
        m_data[index + i] = data[i];
    }
}

MultipleJoinedList::~MultipleJoinedList() {
#pragma acc exit data delete(m_data[:m_size])
    delete[] m_index_list;
    delete[] m_data;
    delete[] m_size_list;
}

size_t *MultipleJoinedList::get_slice(size_t level) {
    size_t slice_size = get_slice_size(level);
    size_t *slice = new size_t[slice_size];
    if (slice_size > 0) {
        size_t first_index = get_first_index(level, 0);
        //TODO(cvm) legit? want to return a copy of a part of 'data'
        std::copy(m_data + first_index, m_data + first_index + slice_size, slice);
    }
    return slice;
}
