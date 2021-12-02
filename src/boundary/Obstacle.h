/// \file       Obstacle.h
/// \brief      Data class of obstacle object
/// \date       Oct 01, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_OBSTACLE_H_
#define ARTSS_BOUNDARY_OBSTACLE_H_

#include <cmath>
#include <string>

#include "Coordinate.h"
#include "BoundaryData.h"
#include "BoundaryDataController.h"
#include "PatchObject.h"
#include "../DomainData.h"
#include "../utility/Utility.h"
#include "../utility/settings/Settings.h"


class Obstacle {
 public:
    Obstacle(Settings::Settings const &settings,
             real x1, real x2,
             real y1, real y2,
             real z1, real z2,
             const std::string &name);
    Obstacle(Settings::Settings const &settings,
             Coordinate<size_t> &coords_start, Coordinate<size_t> &coords_end,
             size_t level,
             const std::string &name);
    ~Obstacle();

    size_t * get_obstacle_list() const { return m_obstacle_list; }
    size_t get_size_obstacle_list() const { return m_size_obstacle_list; }

    size_t ** get_boundary_list() const { return m_boundary; }
    PatchObject *get_size_boundary_list() { return &m_size_boundary; }

    bool is_obstacle_cell(const Coordinate<size_t> &coords) const;

    void print();

    Coordinate<size_t> *get_start_coordinates() { return &m_start; }
    Coordinate<size_t> *get_end_coordinates() { return &m_end; }
    Coordinate<size_t> *get_strides() { return &m_strides; }

    std::string get_name() { return m_name; }

    void replace_patch(size_t *indices, size_t size, Patch p);
    void control();

    bool is_corner_cell(const Coordinate<size_t> &coord) const;
    bool is_edge_cell(const Coordinate<size_t> &coord) const;
    bool has_overlap(size_t i1, size_t i2, size_t j1, size_t j2, size_t k1, size_t k2) const;

    bool static remove_circular_constraints(Obstacle *o1, Obstacle* o2);
private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    Settings::Settings const &m_settings;
    std::string m_name;
    Coordinate<size_t> m_start;
    Coordinate<size_t> m_end;
    Coordinate<size_t> m_strides;

    size_t m_level = 0;

    size_t* m_obstacle_list;
    size_t ** m_boundary;

    size_t m_size_obstacle_list;
    PatchObject m_size_boundary;

    void init();
    void create_obstacle(size_t Nx, size_t Ny);

    void print_details();

    size_t get_size() const { return m_strides[CoordinateAxis::X] * m_strides[CoordinateAxis::Y] * m_strides[CoordinateAxis::Z]; }

    void remove_cells_at_boundary();

    static int get_matching_index(real obstacle_coordinate, real spacing, real start_coordinate);
    static bool has_overlap(size_t o1_coord1, size_t o1_coord2, size_t o2_coord1, size_t o2_coord2);
    static void calculate_area_index(
            Obstacle *o1, Obstacle *o2,
            size_t *o1_coordinate, size_t *o2_coordinate,
            CoordinateAxis coordinate_axis,
            bool start);
    static bool circular_constraints(Obstacle *o1, Obstacle *o2, CoordinateAxis coordinate_axis);
};

#endif /* ARTSS_BOUNDARY_OBSTACLE_H_ */
