/// \file       JoinedListObstacle.cpp
/// \brief      
/// \date       Oct 23, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#include "ObstacleJoinedList.h"

ObstacleJoinedList::ObstacleJoinedList(size_t multigrid_level, size_t number_of_obstacles) : m_number_of_obstacles(number_of_obstacles) {
    m_index_list = new size_t[(multigrid_level + 1) * number_of_obstacles + 1];
    for (size_t id = 0; id < number_of_obstacles; id++) {
        m_index_list[id] = 0;
    }
    m_size_list = new size_t[multigrid_level  * number_of_obstacles];
}

void ObstacleJoinedList::set_size(const size_t size) {
    m_size = size;
    m_data = new size_t[m_size];
}

size_t ObstacleJoinedList::get_slice_size(size_t level) const {
    return get_last_index(level, m_number_of_obstacles - 1) - get_first_index(level, 0) + 1;
}

size_t ObstacleJoinedList::get_slice_size(size_t level, size_t obstacle_id) const {
    return get_last_index(level, obstacle_id) - get_first_index(level, obstacle_id) + 1;
}

size_t ObstacleJoinedList::get_first_index(size_t level, size_t obstacle_id) const {
    return m_index_list[level * m_number_of_obstacles + obstacle_id];
}

/**
 * for data: get last index of slice of level
 * @param level  level of slice
 * @return
 */
size_t ObstacleJoinedList::get_last_index(size_t level, size_t obstacle_id) const {
    size_t size = get_slice_size(level, obstacle_id);
    if (size == 0) {
        return get_first_index(level, obstacle_id);
    }
    return get_first_index(level, obstacle_id) + size - 1;
}

void ObstacleJoinedList::add_data(size_t level, size_t obstacle_id, size_t size, const size_t *data) {
    size_t index = m_index_list[level * m_number_of_obstacles];
    m_index_list[level * m_number_of_obstacles + 1] = index + size;
    m_size_list[level * m_number_of_obstacles + obstacle_id] = size;
    for (size_t i = 0; i < size; i++) {
        m_data[index + i] = data[i];
    }
}

void ObstacleJoinedList::control() {

}

ObstacleJoinedList::~ObstacleJoinedList() {
#pragma acc exit data delete(m_data[:m_size])
    delete[] m_index_list;
    delete[] m_data;
    delete[] m_size_list;
}
