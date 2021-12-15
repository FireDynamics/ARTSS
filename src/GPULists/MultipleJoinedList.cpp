/// \file       MultipleJoinedList.cpp
/// \brief      stores indices lists for GPU, differs between multilevel for multigrid and multiple objects
/// \details    Meant to be for surface and obstacle objects, as there may be multiple objects for each multigrid level
/// \date       Oct 23, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include "MultipleJoinedList.h"

MultipleJoinedList::MultipleJoinedList(size_t multigrid_level, size_t number_of_objects) : m_number_of_objects(number_of_objects) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_index_list = new size_t[(multigrid_level + 1) * number_of_objects + 1];
    std::fill(m_index_list, m_index_list + (multigrid_level + 1) * number_of_objects + 1, 0);

    m_size_list = new size_t[(multigrid_level + 1) * number_of_objects];
    std::fill(m_size_list, m_size_list + (multigrid_level + 1) * number_of_objects, 0);
}

/// set size of GPU array
/// @param size size of GPU array
void MultipleJoinedList::set_size(const size_t size) {
    m_size = size;
    m_data = new size_t[m_size];
}

/// get size of GPU array for a specific level
/// @param level multigrid level
/// @return size of GPU array for a specific level
size_t MultipleJoinedList::get_slice_size(size_t level) const {
    if (m_number_of_objects > 0) {
        return get_last_index(level, m_number_of_objects - 1) - get_first_index(level, 0) + 1;
    }
    return 0;
}

/// get size of GPU array for a specific level with all objects of that level
/// @param level multigrid level
/// @return size of GPU array for a specific level
size_t MultipleJoinedList::get_slice_size(size_t level, size_t obstacle_id) const {
    return m_size_list[level * m_number_of_objects + obstacle_id];
}

/// get starting index of GPU array for a specific level
/// @param level multigrid level
/// @return starting index of GPU array for the specified level
size_t MultipleJoinedList::get_first_index(size_t level, size_t obstacle_id) const {
    return m_index_list[level * m_number_of_objects + obstacle_id];
}

/// get ending index of GPU array for a specific level. if size of level equals 0
/// the first index will be returned instead. Only possible if there aren't any
/// objects, which always results in a return of 0.
/// @param level multigrid level
/// @return ending index of GPU array for the specified level
size_t MultipleJoinedList::get_last_index(size_t level, size_t obstacle_id) const {
    size_t size = get_slice_size(level, obstacle_id);
    if (size == 0) {
        return get_first_index(level, obstacle_id);
    }
    return get_first_index(level, obstacle_id) + size - 1;
}

/// add data (index list) for one multigrid level for one object. data will be copied
/// @param level multigrid level of index list
/// @param object_id id of object, necessary to identify individual objects
/// @param size size of data
/// @param data index list to be added
void MultipleJoinedList::add_data(size_t level, size_t object_id, size_t size, const size_t *data) {
#ifndef BENCHMARKING
    m_logger->debug("MJL: add data for obstacle id={} with level {} and size {}", object_id, level, size);
#endif
    size_t index = m_index_list[level * m_number_of_objects + object_id];
    m_index_list[level * m_number_of_objects + object_id + 1] = index + size;
    m_size_list[level * m_number_of_objects + object_id] = size;
    for (size_t i = 0; i < size; i++) {
        m_data[index + i] = data[i];
    }
}

MultipleJoinedList::~MultipleJoinedList() {
#pragma acc exit data delete(m_data[:m_size])
    delete[] m_index_list;
    if (m_size > 0) {
        delete[] m_data;
    }
    delete[] m_size_list;
}

/// get a copy of the index list of all objects of a specific level
/// @param level multigrid level
/// @return copy of a part of index list. Transfers the ownership to the caller.
///         Caller needs to free the memory
size_t *MultipleJoinedList::get_slice(size_t level) {
    size_t slice_size = get_slice_size(level);
    auto slice = new size_t[slice_size];
    if (slice_size > 0) {
        size_t first_index = get_first_index(level, 0);
        //TODO(cvm) legit? want to return a copy of a part of 'data'
        std::copy(m_data + first_index, m_data + first_index + slice_size, slice);
    }
    return slice;
}
