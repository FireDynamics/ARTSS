/// \file       JoinedList.cpp
/// \brief      
/// \date       Oct 22, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include "SimpleJoinedList.h"

SimpleJoinedList::SimpleJoinedList(size_t multigrid_level) {
    m_index_list = new size_t[multigrid_level + 2];
    m_index_list[0] = 0;
}

void SimpleJoinedList::set_size(const size_t size) {
    m_size = size;
    m_data = new size_t[m_size];
}

size_t SimpleJoinedList::get_slice_size(size_t level) const {
    return m_size_list[level];
}

size_t SimpleJoinedList::get_first_index(size_t level) const {
    return m_index_list[level];
}

/**
 * for data: get last index of slice of level
 * @param level  level of slice
 * @return
 */
size_t SimpleJoinedList::get_last_index(size_t level) const {
    size_t size = get_slice_size(level);
    if (size == 0) {
        return get_first_index(level);
    }
    return get_first_index(level) + size - 1;
}

void SimpleJoinedList::add_data(size_t level, size_t size, const size_t *data) {
    size_t index = m_index_list[level];
    m_index_list[level + 1] = index + size;
    m_size_list[level] = size;
    for (size_t i = 0; i < size; i++) {
        m_data[index + i] = data[i];
    }
}

void SimpleJoinedList::control() {

}
