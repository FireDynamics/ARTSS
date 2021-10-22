/// \file       JoinedList.cpp
/// \brief      
/// \date       Oct 22, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include "JoinedList.h"

JoinedList::JoinedList(size_t multigrid_level) {
    m_index_list = new size_t[multigrid_level + 2];
    m_index_list[0] = 0;
}

void JoinedList::set_size(const size_t size) {
    m_size = size;
    m_data = new size_t[m_size];
}

size_t JoinedList::get_slice_size(size_t level) const {
    return get_last_index(level) - get_first_index(level);
}

size_t JoinedList::get_first_index(size_t level) const {
    return m_index_list[level];
}

/**
 * for data: get last index of slice of level
 * @param level  level of slice
 * @return
 */
size_t JoinedList::get_last_index(size_t level) const {
    return m_index_list[level + 1] - 1;
}

void JoinedList::add_data(size_t level, size_t size, size_t *data) {

}

/**
 * add size of slice of level
 * @param level level of slice
 * @param size  size of slice
 */
void JoinedList::add_size(size_t level, size_t size) {
    m_size_list[level] = size;
}

void control() {

}
