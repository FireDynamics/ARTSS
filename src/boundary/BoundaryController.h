/// \file       BoundaryController.h
/// \brief      Controll class for boundary
/// \date       Oct 01, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_BOUNDARYCONTROLLER_H_
#define ARTSS_BOUNDARY_BOUNDARYCONTROLLER_H_

#include <vector>
#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"
#include "BoundaryDataController.h"
#include "Coordinate.h"
#include "Multigrid.h"
#include "Obstacle.h"
#include "Surface.h"

class BoundaryController {
 public:
    static BoundaryController* getInstance() { return singleton; }
    static BoundaryController* getInstance(Settings::Settings const &settings);

    ~BoundaryController();

    void apply_boundary(Field &field, bool sync = true);

    void print_boundaries();
    void update_lists();

    /// \brief get array of joined domain inner list
    size_t* get_domain_inner_list_level_joined() const { return m_multigrid->get_domain_inner_cells_level_joined(); }
    /// \brief get size of domain inner list for specified level
    size_t get_size_domain_inner_list_level_joined(size_t level) const { return get_slice_size_domain_inner_list_level_joined(level); }
    size_t get_domain_inner_list_level_joined_start(size_t level) const { return m_multigrid->get_start_index_domain_inner_cells_level_joined(level); }
    size_t get_domain_inner_list_level_joined_end(size_t level) const { return m_multigrid->get_end_index_domain_inner_cells_level_joined(level); }

    /// \bierf get array of joined domain list
    size_t* get_domain_list_level_joined() const { return m_multigrid->get_domain_cells_level_joined(); }
    /// \brief get size of domain list for specified level
    size_t get_slice_size_domain_list_level_joined(size_t level) const { return m_multigrid->get_slice_size_domain_cells(level); }
    size_t get_domain_list_level_joined_start(size_t level) const { return m_multigrid->get_start_index_domain_cells_level_joined(level); }
    size_t get_domain_list_level_joined_end(size_t level) const { return m_multigrid->get_end_index_domain_cells_level_joined(level); }

    size_t* get_obstacle_list_level_joined() const { return m_multigrid->get_obstacle_cells_level_joined(); }
    size_t get_slice_size_obstacle_list_level_joined(size_t level) const { return m_multigrid->get_slice_size_obstacle_cells(level); }
    size_t get_obstacle_list_level_joined_start(size_t level) const { return m_multigrid->get_start_index_obstacle_cells_level_joined(level); }
    size_t get_obstacle_list_level_joined_end(size_t level) const { return m_multigrid->get_end_index_obstacle_cells_level_joined(level); }

    size_t get_size_surfaceList() const {return m_size_surface_list;}

    std::vector<FieldType> get_used_fields() const;

    bool inline is_obstacle_cell(const size_t level, const Coordinate<size_t> &coords) {
        return m_multigrid->is_obstacle_cell(level, coords);
    }

 private:
    explicit BoundaryController(Settings::Settings const &settings);
    size_t get_slice_size_domain_inner_list_level_joined(size_t level) const { return m_multigrid->get_slice_size_domain_inner_cells_level_joined(level); }  // get size of domain inner list

    Settings::Settings const &m_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    static BoundaryController* singleton;

    BoundaryDataController *m_bdc_boundary;
    BoundaryDataController **m_bdc_obstacles;
    Multigrid* m_multigrid;

    Surface** m_surface_list;
    size_t m_number_of_surfaces = 0;
    Obstacle** m_obstacle_list;
    size_t m_number_of_obstacles = 0;

    size_t m_size_surface_list = 0;

    bool m_has_obstacles;
    bool m_has_surfaces;

    void read_XML();
    void parse_boundary_parameter(const std::vector<Settings::BoundarySetting>& boundaries);
    void parse_obstacle_parameter(const std::vector<Settings::ObstacleSetting>& obstacles);
    void parse_surface_parameter(const std::vector<Settings::SurfaceSetting>& surfaces);

    void detect_neighbouring_obstacles();
};

#endif /* ARTSS_BOUNDARY_BOUNDARYCONTROLLER_H_ */

