/// \file       Multigrid.h
/// \brief      Creates all lists needed for multigrid
/// \date       Oct 01, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DOMAIN_MULTIGRID_H_
#define ARTSS_DOMAIN_MULTIGRID_H_

#include "Domain.h"
#include "BoundaryDataController.h"
#include "Obstacle.h"
#include "Surface.h"
#include "DomainData.h"
#include "../field/Field.h"
#include "../utility/Utility.h"
#include "../GPULists/SingleJoinedList.h"
#include "../GPULists/MultipleJoinedList.h"

class Multigrid {
 public:
    Multigrid(const std::vector<Surface> &surfaces, const std::vector<BoundaryDataController> &bdc_surfaces,
              const std::vector<Obstacle> &obstacles, const std::vector<BoundaryDataController> &bdc_obstacles,
              const BoundaryDataController &bdc_domain,
              size_t multigrid_level);
    ~Multigrid();

    // getter -- domain inner cells
    size_t* get_domain_inner_cells_level_joined() const { return m_jl_domain_inner_list.get_data(); }
    size_t get_slice_size_domain_inner_cells_level_joined(size_t level) const { return m_jl_domain_inner_list.get_slice_size(level); };
    size_t get_size_domain_inner_cells_level_joined() const { return m_jl_domain_inner_list.get_size(); }
    size_t get_start_index_domain_inner_cells_level_joined(size_t level) const { return m_jl_domain_inner_list.get_first_index(level); }
    size_t get_end_index_domain_inner_cells_level_joined(size_t level) const { return m_jl_domain_inner_list.get_last_index(level); }

    // getter -- domain cells (boundary + inner)
    size_t* get_domain_cells_level_joined() const { return m_jl_domain_list.get_data(); }
    size_t get_slice_size_domain_cells(size_t level) const { return m_jl_domain_list.get_slice_size(level); }
    size_t get_size_domain_cells() const { return m_jl_domain_list.get_size(); }
    size_t get_start_index_domain_cells_level_joined(size_t level) const { return m_jl_domain_list.get_first_index(level); }
    size_t get_end_index_domain_cells_level_joined(size_t level) const { return m_jl_domain_list.get_last_index(level); }

    // getter -- obstacle cells (boundary + inner)
    size_t* get_obstacle_cells_level_joined() const { return m_jl_obstacle_list.get_data(); }
    size_t get_slice_size_obstacle_cells(size_t level) const { return m_jl_obstacle_list.get_slice_size(level); }
    size_t get_size_obstacle_cells() const { return m_jl_obstacle_list.get_size(); }
    size_t get_start_index_obstacle_cells_level_joined(size_t level) const { return m_jl_obstacle_list.get_first_index(level); }
    size_t get_end_index_obstacle_cells_level_joined(size_t level) const { return m_jl_obstacle_list.get_last_index(level); }

    void replace_obstacles(const std::vector<Obstacle> &obstacles, const std::vector<BoundaryDataController> &bdc_obstacles);

    void update_lists();

    void apply_boundary_condition(Field &field, bool sync = false);

    bool is_obstacle_cell(size_t level, const Coordinate<size_t> &coords);

    std::vector<FieldType> get_used_fields();

    bool is_blocked_by_obstacle(const Coordinate<size_t> &start, const Coordinate<size_t> &end) const;

private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    size_t m_multigrid_levels;
    // all surfaces divided by level
    std::vector<std::vector<Surface>> m_MG_surface_object_list;  // m_MG_surface_object_list[level][surfaceID]
    size_t m_number_of_surface_objects = 0;
    // all obstacles divided by level
    std::vector<std::vector<Obstacle>> m_MG_obstacle_object_list;  // m_MG_obstacle_object_list[level][obstacleID]
    size_t m_number_of_obstacle_objects = 0;
    // boundary for each level
    std::vector<Domain> m_MG_domain_object_list;  // m_MG_boundary_object_list[level]

    //---- all level joined / arrays for GPU -----
    SingleJoinedList m_jl_domain_list;
    SingleJoinedList m_jl_domain_inner_list;
    SingleJoinedList **m_jl_domain_boundary_list_patch_divided;  // [Patch]

    SingleJoinedList m_jl_obstacle_list;  // all obstacle boundary cells
    MultipleJoinedList **m_jl_obstacle_boundary_list_patch_divided;  // [Patch]

    MultipleJoinedList **m_jl_surface_list_patch_divided;  // [Patch]

    void create_multigrid_obstacle_lists();
    void create_multigrid_surface_lists();
    void create_multigrid_domain_lists();

    void send_surface_lists_to_GPU();
    void send_obstacle_lists_to_GPU();
    void send_domain_lists_to_GPU();

    size_t surface_dominant_restriction(size_t level, PatchObject *sum_patches);
    size_t obstacle_dominant_restriction(size_t level, PatchObject *sum_patches, size_t **tmp_store_obstacle);

    void print();

    BoundaryDataController m_bdc_domain;
    std::vector<BoundaryDataController> m_bdc_obstacle;
    std::vector<BoundaryDataController> m_bdc_surface;

    void remove_domain_lists_from_GPU();

};

#endif /* ARTSS_DOMAIN_MULTIGRID_H_*/

