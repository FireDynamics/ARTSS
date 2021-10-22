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
#include "../utility/tinyxml2.h"
#include "../utility/Utility.h"
#include "BoundaryDataController.h"
#include "Multigrid.h"
#include "Obstacle.h"
#include "Surface.h"

class BoundaryController {
 public:
    static BoundaryController* getInstance();
    ~BoundaryController();

    void apply_boundary(Field &field, bool sync = true);

    void print_boundaries();
    void update_lists();

    size_t get_size_inner_list() const;
    size_t get_size_boundary_list() const;
    size_t* get_obstacle_list() const;
    size_t get_size_obstacle_list() const;

    size_t* get_inner_list_level_joined() const;
    size_t get_size_inner_list_level_joined() const;
    size_t get_inner_list_level_joined_start(size_t level) const;
    size_t get_inner_list_level_joined_end(size_t level) const;

    size_t* get_boundary_list_level_joined() const;
    size_t get_size_boundary_list_level_joined() const;
    size_t get_boundary_list_level_joined_start(size_t level) const;
    size_t get_boundary_list_level_joined_end(size_t level) const;

    size_t get_size_surfaceList() const {return m_size_surface_list;}

    size_t get_obstacle_stride_x(size_t id, size_t level) const;
    size_t get_obstacle_stride_y(size_t id, size_t level) const;
    size_t get_obstacle_stride_z(size_t id, size_t level) const;

    std::vector<FieldType> get_used_fields() const;

    bool inline is_obstacle_cell(const size_t level, const size_t idx) {
        return m_multigrid->is_obstacle_cell(level, idx);
    }

    bool inline is_obstacle_cell(const size_t level,
                          const size_t i, const size_t j, const size_t k) {
        return m_multigrid->is_obstacle_cell(level, i, j, k);
    }

private:
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

    BoundaryController();
    void read_XML();
    void parse_boundary_parameter(tinyxml2::XMLElement *xml_parameter);
    void parse_obstacle_parameter(tinyxml2::XMLElement *xml_parameter);
    void parse_surface_parameter(tinyxml2::XMLElement *xml_parameter);

    void detect_neighbouring_obstacles();
};

#endif /* ARTSS_BOUNDARY_BOUNDARYCONTROLLER_H_ */

