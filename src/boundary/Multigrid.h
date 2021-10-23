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

    size_t get_size_obstacle_list() const { return m_jl_obstacle_list.get_size(); }
    size_t *get_obstacle_list() const { return m_jl_obstacle_list.get_data(); }

    // getter -- domain inner cell
    size_t* get_domain_inner_cells_level_joined() const { return m_jl_domain_inner_list.get_data(); }
    size_t get_slice_size_domain_inner_cells_level_joined(size_t level) const { return m_jl_domain_inner_list.get_slice_size(level); };
    size_t get_size_domain_inner_cells_level_joined() const { return m_jl_domain_inner_list.get_size(); }
    size_t get_start_index_domain_inner_cells_level_joined(size_t level) const { return m_jl_domain_inner_list.get_first_index(level); }
    size_t get_end_index_domain_inner_cells_level_joined(size_t level) const { return m_jl_domain_inner_list.get_last_index(level); }

    // getter -- domain boundary cell
    size_t* get_domain_boundary_cells_level_joined() const { return m_jl_domain_boundary_list.get_data(); }
    size_t get_size_domain_boundary_cells_level_joined() const { return m_jl_domain_boundary_list.get_size(); }
    size_t get_start_index_domain_boundary_cells_level_joined(size_t level) const { return m_jl_domain_boundary_list.get_first_index(level); }
    size_t get_end_index_domain_boundary_cells_level_joined(size_t level) const { return m_jl_domain_boundary_list.get_last_index(level); }

    void update_lists();

    void apply_boundary_condition(Field &field, bool sync = false);

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

    // surface indices divided by level
    size_t** m_MG_surface_index_list;

    // start index of each level in level joined list
    size_t* m_size_MG_surface_index_list_level;

    //---- all level joined / arrays for GPU -----
    SimpleJoinedList m_jl_domain_inner_list;
    SimpleJoinedList m_jl_domain_boundary_list;
    SimpleJoinedList **m_jl_domain_boundary_list_patch_divided;  // [Patch]

    SimpleJoinedList m_jl_obstacle_list;
    ObstacleJoinedList **m_jl_obstacle_boundary_list_patch_divided;  // [Patch]

    size_t* m_data_MG_boundary_list_level_joined;
    size_t* m_data_MG_surface_list_level_joined;

    size_t get_length_of_surface_index_list_joined() const;
    size_t get_first_index_of_surface_index_list(size_t level) const;
    size_t get_last_index_of_surface_index_list(size_t level) const;

    void init_domain();
    void init_obstacles(Obstacle **obstacle_list);
    void init_surfaces(Surface **surface_list);

    void add_MG_lists();
    size_t calc_obstacles(Obstacle** obstacle_object_list, PatchObject *patch_sum) const;
    void calc_surfaces(Surface** surface_list);
    void send_lists_to_GPU();
    void create_domain_lists_for_GPU();
    void send_surface_lists_to_GPU();
    void send_obstacle_lists_to_GPU(const size_t *size_MG_obstacle_index_list_level, size_t **MG_obstacle_index_list);

    void surface_dominant_restriction(size_t level);
    size_t obstacle_dominant_restriction(size_t level, PatchObject *sum_patches);

    void control();
    void print();

    // size_t **m_data_surfaces_patches_joined;
    BoundaryDataController *m_bdc_boundary;
    BoundaryDataController **m_bdc_obstacle;

    void remove_boundary_lists_from_GPU();
};

#endif /* ARTSS_BOUNDARY_MULTIGRID_H_*/

