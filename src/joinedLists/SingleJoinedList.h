/// \file       JoinedList.h
/// \brief      
/// \date       Oct 22, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich All rights reserved.
//
#ifndef ARTSS_JOINEDLISTS_JOINEDLISTDOMAIN_H_
#define ARTSS_JOINEDLISTS_JOINEDLISTDOMAIN_H_

#include "../utility/Utility.h"

class SingleJoinedList {
  public:
    explicit SingleJoinedList(size_t multigrid_level);
    ~SingleJoinedList();

    void copyin() {
#pragma acc enter data copyin(m_data[:m_size])
#ifdef GPU_DEBUG
        m_gpu_logger->debug("copyin gpu index list with data pointer: {} and size: {}", static_cast<void *>(m_data), m_size);
#endif
    }

    inline size_t &operator[](size_t i) const { return m_data[i]; }

    size_t get_first_index(size_t level) const;
    size_t get_last_index(size_t level) const;
    size_t get_slice_size(size_t level) const;
    size_t get_size() const { return m_size; }
    size_t * get_data() const { return m_data; }
    size_t * get_slice(size_t level);

    void set_size(size_t size);
    void add_data(size_t level, size_t size, const size_t *data);
  private:
    size_t *m_data;  // array for GPU
    size_t m_size = 0;
    size_t *m_size_list;  // size of each level
    size_t *m_index_list;  // starting index of each level
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
#ifdef GPU_DEBUG
    std::shared_ptr<spdlog::logger> m_gpu_logger;
#endif
};


#endif /* ARTSS_JOINEDLISTS_JOINEDLISTDOMAIN_H_ */
