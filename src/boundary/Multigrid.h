/// \file       Multigrid.h
/// \brief      Creates all lists needed for multigrid
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_MULTIGRID_H_
#define ARTSS_BOUNDARY_MULTIGRID_H_

#include "Domain.h"
#include "BoundaryDataController.h"
#include "Obstacle.h"
#include "Surface.h"
#include "../DomainData.h"
#include "../field/Field.h"
#include "../utility/Utility.h"
#include "../joinedLists/SimpleJoinedList.h"
#include "../joinedLists/ObstacleJoinedList.h"

class Multigrid {
 public:
    Multigrid(
            size_t number_of_surfaces, Surface** surface_list,
            size_t number_of_obstacles, Obstacle** obstacle_list,
            BoundaryDataController* bdc_boundary,
            BoundaryDataController **bdc_obstacles,
            size_t multigrid_level);
    explicit Multigrid(BoundaryDataController *bdc_boundary, size_t multigrid_level);
    ~Multigrid();

    size_t get_size_domain_inner_cells(size_t level = 0) const;
    size_t get_size_domain_boundary_cells(size_t level = 0) const;
    size_t get_size_obstacle_list() const;
    size_t *get_obstacle_boundary_list() const;

    size_t* get_inner_list_level_joined() const { return m_jl_domain_inner_list.get_data(); }
    size_t get_size_inner_list_level_joined() const { return *(m_size_MG_inner_index_list_level + m_multigrid_levels + 1); }
    size_t get_inner_list_level_joined_start(size_t level) const;
    size_t get_inner_list_level_joined_end(size_t level) const;

    size_t* get_boundary_list_level_joined() const { return m_data_MG_boundary_list_level_joined; }
    size_t get_size_boundary_list_level_joined() const { return *(m_size_MG_boundary_index_list_level + m_multigrid_levels + 1); }
    size_t get_boundary_list_level_joined_start(size_t level) const;
    size_t get_boundary_list_level_joined_end(size_t level) const;

    void update_lists();

    void apply_boundary_condition(Field &field, bool sync = false);

    size_t get_obstacle_stride_x(size_t id, size_t level) const;
    size_t get_obstacle_stride_y(size_t id, size_t level) const;
    size_t get_obstacle_stride_z(size_t id, size_t level) const;

    bool is_obstacle_cell(size_t level, size_t index);
    bool is_obstacle_cell(size_t level, size_t i, size_t j, size_t k);

private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
#ifdef GPU_DEBUG
    std::shared_ptr<spdlog::logger> m_gpu_logger;
#endif
    size_t m_multigrid_levels;
    // all surfaces divided by level
    Surface*** m_MG_surface_object_list;  // m_MG_surface_object_list[level][surfaceID]
    size_t m_number_of_surface_objects = 0;
    // all obstacles divided by level
    Obstacle*** m_MG_obstacle_object_list;  // m_MG_obstacle_object_list[level][obstacleID]
    size_t m_number_of_obstacle_objects = 0;
    // boundary for each level
    Domain** m_MG_domain_object_list;  // m_MG_boundary_object_list[level]

    // obstacle indices divided by level
    size_t** m_MG_obstacle_index_list;
    // surface indices divided by level
    size_t** m_MG_surface_index_list;

    // start index of each level in level joined list
    size_t* m_size_MG_inner_index_list_level;
    size_t* m_size_MG_boundary_index_list_level;
    size_t* m_size_MG_obstacle_index_list_level;
    size_t* m_size_MG_surface_index_list_level;

    // start index of each boundary object in level joined list
    // (slice z = Front/Back, slice y = Bottom/Top, slice x = Left/Right)
    size_t* m_size_MG_boundary_slice_z_level;
    size_t* m_size_MG_boundary_slice_y_level;
    size_t* m_size_MG_boundary_slice_x_level;

    // start index of each obstacle object in level joined list
    size_t* m_size_MG_obstacle_front_level;
    size_t* m_size_MG_obstacle_back_level;
    size_t* m_size_MG_obstacle_bottom_level;
    size_t* m_size_MG_obstacle_top_level;
    size_t* m_size_MG_obstacle_left_level;
    size_t* m_size_MG_obstacle_right_level;

    //---- all level joined / arrays for GPU -----
    SimpleJoinedList m_jl_domain_inner_list;
    SimpleJoinedList **m_jl_domain_boundary_list;

    SimpleJoinedList m_jl_obstacle_list;
    ObstacleJoinedList **m_jl_obstacle_boundary_list;

    size_t* m_data_MG_boundary_list_level_joined;
    size_t* m_data_MG_surface_list_level_joined;

    size_t get_last_index_of_obstacle_front(size_t level, size_t id) const;
    size_t get_last_index_of_obstacle_back(size_t level, size_t id) const;
    size_t get_last_index_of_obstacle_bottom(size_t level, size_t id) const;
    size_t get_last_index_of_obstacle_top(size_t level, size_t id) const;
    size_t get_last_index_of_obstacle_left(size_t level, size_t id) const;
    size_t get_last_index_of_obstacle_right(size_t level, size_t id) const;
    size_t get_first_index_of_obstacle_front(size_t level, size_t id) const;
    size_t get_first_index_of_obstacle_back(size_t level, size_t id) const;
    size_t get_first_index_of_obstacle_bottom(size_t level, size_t id) const;
    size_t get_first_index_of_obstacle_top(size_t level, size_t id) const;
    size_t get_first_index_of_obstacle_left(size_t level, size_t id) const;
    size_t get_first_index_of_obstacle_right(size_t level, size_t id) const;
    size_t get_length_of_obstacle_front(size_t level) const;
    size_t get_length_of_obstacle_back(size_t level) const;
    size_t get_length_of_obstacle_bottom(size_t level) const;
    size_t get_length_of_obstacle_top(size_t level) const;
    size_t get_length_of_obstacle_left(size_t level) const;
    size_t get_length_of_obstacle_right(size_t level) const;

    // get length of bSlice from/for joined array
    size_t get_length_of_boundary_slice_z_joined() const;
    size_t get_length_of_boundary_slice_y_joined() const;
    size_t get_length_of_boundary_slice_x_joined() const;

    size_t get_first_index_of_boundary_slice_z(size_t level) const;
    size_t get_first_index_of_boundary_slice_x(size_t level) const;
    size_t get_first_index_of_boundary_slice_y(size_t level) const;
    size_t get_last_index_of_boundary_slice_z(size_t level) const;
    size_t get_last_index_of_boundary_slice_x(size_t level) const;
    size_t get_last_index_of_boundary_slice_y(size_t level) const;

    size_t get_length_of_boundary_index_list_joined() const;
    size_t get_length_of_surface_index_list_joined() const;
    size_t get_first_index_of_boundary_index_list(size_t level) const;
    size_t get_first_index_of_surface_index_list(size_t level) const;
    size_t get_last_index_of_inner_index_list(size_t level) const;
    size_t get_last_index_boundary_index_list(size_t level) const;
    size_t get_last_index_of_surface_index_list(size_t level) const;

    void init();
    void add_MG_lists();
    size_t calc_obstacles(Obstacle** obstacle_object_list, PatchObject *patch_sum);
    void calc_surfaces(Surface** surface_list);
    void send_lists_to_GPU();
    void send_domain_lists_to_GPU();
    void send_surface_lists_to_GPU();
    void send_obstacle_lists_to_GPU();

    void surface_dominant_restriction(size_t level);
    size_t obstacle_dominant_restriction(size_t level, PatchObject *sum_patches);

    void control();
    void print();

    // size_t **m_data_surfaces_patches_joined;
    size_t **m_data_obstacles_patches_joined;
    BoundaryDataController *m_bdc_boundary;
    BoundaryDataController **m_bdc_obstacle;

    void remove_boundary_lists_from_GPU();

};

#endif /* ARTSS_BOUNDARY_MULTIGRID_H_*/

