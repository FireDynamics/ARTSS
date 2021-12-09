/// \file       Multigrid.h
/// \brief      Creates all lists needed for multigrid
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_MULTIGRID_H_
#define ARTSS_BOUNDARY_MULTIGRID_H_

#include "Boundary.h"
#include "BoundaryDataController.h"
#include "Obstacle.h"
#include "Surface.h"
#include "../Domain.h"
#include "../field/Field.h"
#include "../utility/Utility.h"

class Multigrid {
 public:
    Multigrid(
            Settings::Settings const &settings,
            size_t number_of_surfaces, Surface** surface_list,
            size_t number_of_obstacles, Obstacle** obstacle_list,
            BoundaryDataController* bdc_boundary,
            BoundaryDataController **bdc_obstacles);
    Multigrid(Settings::Settings const &settings, BoundaryDataController *bdc_boundary);
    ~Multigrid();

    size_t get_size_inner_list(size_t level = 0) const;
    size_t get_size_boundary_list(size_t level = 0) const;
    size_t get_size_obstacle_list() const;
    size_t *get_obstacle_list() const;

    size_t* get_inner_list_level_joined() const { return m_data_MG_inner_list_level_joined; }
    size_t get_size_inner_list_level_joined() const { return *(m_size_MG_inner_index_list_level + m_levels + 1); }
    size_t get_inner_list_level_joined_start(size_t level) const;
    size_t get_inner_list_level_joined_end(size_t level) const;

    size_t* get_boundary_list_level_joined() const { return m_data_MG_boundary_list_level_joined; }
    size_t get_size_boundary_list_level_joined() const { return *(m_size_MG_boundary_index_list_level + m_levels + 1); }
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
    Settings::Settings const &m_settings;
    size_t m_levels;
    // all surfaces divided by level
    Surface*** m_MG_surface_object_list;  // m_MG_surface_object_list[level][surfaceID]
    size_t m_number_of_surface_objects = 0;
    // all obstacles divided by level
    Obstacle*** m_MG_obstacle_object_list;  // m_MG_obstacle_object_list[level][obstacleID]
    size_t m_number_of_obstacle_objects = 0;
    // boundary for each level
    Boundary** m_MG_boundary_object_list;  // m_MG_boundary_object_list[level]

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
    size_t* m_data_MG_inner_list_level_joined;
    size_t* m_data_MG_boundary_list_level_joined;
    size_t* m_data_MG_surface_list_level_joined;

    size_t* m_data_MG_boundary_front_level_joined;
    size_t* m_data_MG_boundary_back_level_joined;
    size_t* m_data_MG_boundary_bottom_level_joined;
    size_t* m_data_MG_boundary_top_level_joined;
    size_t* m_data_MG_boundary_left_level_joined;
    size_t* m_data_MG_boundary_right_level_joined;

    size_t* m_data_MG_obstacle_front_level_joined;
    size_t* m_data_MG_obstacle_back_level_joined;
    size_t* m_data_MG_obstacle_bottom_level_joined;
    size_t* m_data_MG_obstacle_top_level_joined;
    size_t* m_data_MG_obstacle_left_level_joined;
    size_t* m_data_MG_obstacle_right_level_joined;
    size_t* m_data_MG_obstacle_list_zero_joined;

    size_t get_size_obstacle_index_list(size_t level) const;
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

    // get length of from/for joined array
    size_t get_length_of_obstacle_front_joined() const;
    size_t get_length_of_obstacle_back_joined() const;
    size_t get_length_of_obstacle_bottom_joined() const;
    size_t get_length_of_obstacle_top_joined() const;
    size_t get_length_of_obstacle_left_joined() const;
    size_t get_length_of_obstacle_right_joined() const;

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

    size_t get_length_of_inner_index_list_joined() const;
    size_t get_length_of_boundary_index_list_joined() const;
    size_t get_length_of_surface_index_list_joined() const;
    size_t get_first_index_of_inner_index_list(size_t level) const;
    size_t get_first_index_of_boundary_index_list(size_t level) const;
    size_t get_first_index_of_surface_index_list(size_t level) const;
    size_t get_last_index_of_inner_index_list(size_t level) const;
    size_t get_last_index_boundary_index_list(size_t level) const;
    size_t get_last_index_of_surface_index_list(size_t level) const;

    void init();
    void add_MG_lists();
    void calc_obstacles(Obstacle** obstacle_object_list);
    void calc_surfaces(Surface** surface_list);
    void send_lists_to_GPU();
    void send_boundary_lists_to_GPU();
    void send_surface_lists_to_GPU();
    void send_obstacle_lists_to_GPU();

    void surface_dominant_restriction(size_t level);
    Obstacle **obstacle_dominant_restriction(size_t level);

    void control();
    void print();

    size_t **m_data_boundary_patches_joined;
    // size_t **m_data_surfaces_patches_joined;
    size_t **m_data_obstacles_patches_joined;
    BoundaryDataController *m_bdc_boundary;
    BoundaryDataController **m_bdc_obstacle;

    void remove_boundary_lists_from_GPU();

    static bool control_obstacle_overlap(
            Obstacle* o,
            size_t *i1, size_t *i2,
            size_t *j1, size_t *j2,
            size_t *k1, size_t *k2);
};

#endif /* ARTSS_BOUNDARY_MULTIGRID_H_*/

