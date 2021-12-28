/// \file       DomainController.h
/// \brief      Controller class for domain, grants access to indices lists
/// \date       Oct 01, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_DOMAIN_DOMAINCONTROLLER_H_
#define ARTSS_DOMAIN_DOMAINCONTROLLER_H_

#include <vector>
#include "../field/Field.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Utility.h"
#include "BoundaryDataController.h"
#include "Coordinate.h"
#include "Multigrid.h"
#include "Obstacle.h"
#include "Surface.h"
using return_surface = std::tuple<std::vector<Surface>, std::vector<BoundaryDataController>>;
using return_obstacle = std::tuple<std::vector<Obstacle>, std::vector<BoundaryDataController>>;
using return_xml_objects = std::tuple<BoundaryDataController, std::vector<Surface>, std::vector<BoundaryDataController>, std::vector<Obstacle>, std::vector<BoundaryDataController>>;

class DomainController {
 public:
    static DomainController* getInstance() { return singleton; }
    static DomainController* getInstance(Settings::Settings const &settings);

    ~DomainController();

    void apply_boundary(Field &field, bool sync = true);

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

    [[nodiscard]] std::vector<FieldType> get_used_fields() const;

    bool inline is_obstacle_cell(const size_t level, const Coordinate<size_t> &coords) {
        return m_multigrid->is_obstacle_cell(level, coords);
    }

 private:
    explicit DomainController(Settings::Settings const &settings);
    size_t get_slice_size_domain_inner_list_level_joined(size_t level) const { return m_multigrid->get_slice_size_domain_inner_cells_level_joined(level); }  // get size of domain inner list

    Settings::Settings const &m_settings;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    static DomainController* singleton;

    Multigrid* m_multigrid;

    bool m_has_obstacles;
    bool m_has_surfaces;

    return_xml_objects read_XML();
    return_obstacle parse_obstacle_parameter(const std::vector<Settings::ObstacleSetting>& obstacle_setting);
    return_surface parse_surface_parameter(const std::vector<Settings::SurfaceSetting>& surface_setting);

    void detect_neighbouring_obstacles(std::vector<Obstacle> &obstacle_list);
    void print_boundaries(const BoundaryDataController &bdc_domain, const std::vector<BoundaryDataController> &bdc_surfaces, const std::vector<BoundaryDataController> &bdc_obstacles);

};

#endif /* ARTSS_DOMAIN_DOMAINCONTROLLER_H_ */

