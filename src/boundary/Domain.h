/// \file       Boundary.h
/// \brief      Data class of boundary object
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_DOMAIN_H_
#define ARTSS_BOUNDARY_DOMAIN_H_


#include "BoundaryData.h"
#include "Obstacle.h"
#include "../utility/Utility.h"
#include "PatchObject.h"

class Domain {
 public:
    ~Domain();
    Domain(
            Obstacle** obstacle_list,
            size_t number_of_obstacles,
            size_t size_obstacles,
            size_t multigrid_level = 0);
    explicit Domain(size_t multigrid_level = 0);
    void init(size_t size_obstacles);

    size_t * get_domain_list() const { return m_domain_list; }
    size_t ** get_boundary_list() const { return m_boundary; }
    //size_t* get_boundary_front() const { return m_boundary[FRONT]; }
    //size_t* get_boundary_back() const { return m_boundary[BACK]; }
    //size_t* get_boundary_top() const { return m_boundary[TOP]; }
    //size_t* get_boundary_bottom() const { return m_boundary[BOTTOM]; }
    //size_t* get_boundary_left() const { return m_boundary[LEFT]; }
    //size_t* get_boundary_right() const { return m_boundary[RIGHT]; }

    size_t * get_inner_list() const { return m_inner_list; }
    size_t get_size_inner_list() const { return m_size_inner_list; }

    size_t get_size_domain_list() const { return m_size_domain_list; }
    PatchObject & get_size_boundary_list() { return m_size_boundary; }
    //size_t get_size_boundary_front() const { return m_size_boundary_front; }
    //size_t get_size_boundary_back() const { return m_size_boundary_back; }
    //size_t get_size_boundary_top() const { return m_size_boundary_top; }
    //size_t get_size_boundary_bottom() const { return m_size_boundary_bottom; }
    //size_t get_size_boundary_left() const { return m_size_boundary_left; }
    //size_t get_size_boundary_right() const { return m_size_boundary_right; }

    void update_lists(Obstacle **obstacle_list, size_t number_of_obstacles, size_t size_obstacles);
    void update_lists();
    void control(size_t size_obstacles);

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif

    size_t m_multigrid_level;

    size_t *m_domain_list;
    size_t **m_boundary;

    PatchObject m_size_boundary;
    size_t m_size_domain_list;

    size_t *m_inner_list;
    size_t m_size_inner_list;

    void boundary_cells();
    void inner_cells(Obstacle **obstacle_list, size_t number_of_obstacles);
    void inner_cells();
    void print(size_t size_obstacles);
    void clear_lists();
};

#endif /* ARTSS_BOUNDARY_DOMAIN_H_ */
