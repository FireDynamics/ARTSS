/// \file       JoinedListObstacle.h
/// \brief      
/// \date       Oct 23, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_JOINEDLISTS_JOINEDLISTOBSTACLE_H_
#define ARTSS_JOINEDLISTS_JOINEDLISTOBSTACLE_H_

#include "../utility/Utility.h"

class ObstacleJoinedList {
  public:
    ObstacleJoinedList(size_t multigrid_level, size_t m_number_of_obstacles);
    inline size_t &operator[](size_t i) const { return m_data[i]; }
    size_t get_first_index(size_t level, size_t obstacle_id) const;
    size_t get_last_index(size_t level, size_t obstacle_id) const;
    size_t get_slice_size(size_t level) const;
    size_t get_slice_size(size_t level, size_t obstacle_id) const;
    size_t get_size() const { return m_size; }
    void set_size(size_t size);
    void add_data(size_t level, size_t obstacle_id, size_t size, const size_t *data);
    size_t * get_data() const { return m_data; }
  private:
    void control();
    size_t *m_data;
    size_t m_size;
    size_t *m_size_list;  // size of each level
    size_t *m_index_list;  // starting index of each level
    size_t m_number_of_obstacles;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
#ifdef GPU_DEBUG
    std::shared_ptr<spdlog::logger> m_gpu_logger;
#endif
};


#endif /* ARTSS_JOINEDLISTS_JOINEDLISTOBSTACLE_H_ */
