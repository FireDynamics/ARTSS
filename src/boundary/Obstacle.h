/// \file       Obstacle.h
/// \brief      Data class of obstacle object
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_OBSTACLE_H_
#define ARTSS_BOUNDARY_OBSTACLE_H_

#include <cmath>
#include <string>

#include "BoundaryData.h"
#include "BoundaryDataController.h"
#include "../Domain.h"
#include "../utility/tinyxml2.h"
#include "../utility/Utility.h"


class Obstacle {
 public:
    Obstacle(real x1, real x2, real y1, real y2, real z1, real z2, const std::string &name);
    Obstacle(
            size_t coords_i1, size_t coords_j1, size_t coords_k1,
            size_t coords_i2, size_t coords_j2, size_t coords_k2,
            size_t level,
            const std::string& name);
    ~Obstacle();

    size_t* get_obstacle_list() const { return m_obstacle_list; }
    size_t* get_obstacle_front() const { return m_obstacle_front; }
    size_t* get_obstacle_back() const { return m_obstacle_back; }
    size_t* get_obstacle_top() const { return m_obstacle_top; }
    size_t* get_obstacle_bottom() const { return m_obstacle_bottom; }
    size_t* get_obstacle_left() const { return m_obstacle_left; }
    size_t* get_obstacle_right() const { return m_obstacle_right; }

    size_t get_size_obstacle_list() const { return m_size_obstacle_list; }

    size_t get_size_obstacle_front() const { return m_size_obstacle_front; }
    size_t get_size_obstacle_back() const { return m_size_obstacle_back; }
    size_t get_size_obstacle_bottom() const { return m_size_obstacle_bottom; }
    size_t get_size_obstacle_top() const { return m_size_obstacle_top; }
    size_t get_size_obstacle_left() const { return m_size_obstacle_left; }
    size_t get_size_obstacle_right() const { return m_size_obstacle_right; }

    bool is_obstacle_cell(size_t i, size_t j, size_t k) const;

    void print();

    size_t get_coordinates_i1() const { return m_i1; }
    size_t get_coordinates_j1() const { return m_j1; }
    size_t get_coordinates_k1() const { return m_k1; }
    size_t get_coordinates_i2() const { return m_i2; }
    size_t get_coordinates_j2() const { return m_j2; }
    size_t get_coordinates_k2() const { return m_k2; }

    size_t get_stride_x() const { return m_i2 - m_i1 + 1; }
    size_t get_stride_y() const { return m_j2 - m_j1 + 1; }
    size_t get_stride_z() const { return m_k2 - m_k1 + 1; }

    std::string get_name() { return m_name; }

    bool static remove_circular_constraints(Obstacle *o1, Obstacle* o2);
    void replace_patch(size_t *indices, size_t size, Patch p);
    void control();

    void set_inner_cells(Field *f, real value);

    bool is_corner_cell(size_t i, size_t j, size_t k);
    bool is_edge_cell(size_t i, size_t j, size_t k);
    bool has_overlap(size_t i1, size_t i2, size_t j1, size_t j2, size_t k1, size_t k2);

 private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    std::string m_name;
    size_t m_i1, m_j1, m_k1;
    size_t m_i2, m_j2, m_k2;

    size_t m_level = 0;

    size_t* m_obstacle_list;
    size_t* m_obstacle_front;
    size_t* m_obstacle_back;
    size_t* m_obstacle_top;
    size_t* m_obstacle_bottom;
    size_t* m_obstacle_left;
    size_t* m_obstacle_right;

    size_t m_size_obstacle_list;

    // actual size
    size_t m_size_obstacle_front;
    size_t m_size_obstacle_back;
    size_t m_size_obstacle_top;
    size_t m_size_obstacle_bottom;
    size_t m_size_obstacle_left;
    size_t m_size_obstacle_right;

    void init(size_t level);
    void create_obstacle(size_t Nx, size_t Ny);

    static real match_grid(real obstacle_coordinate, real spacing, real start_coordinate);
    static int get_matching_index(real obstacle_coordinate, real spacing, real start_coordinate);

    void print_details();

    static bool has_overlap(size_t o1_coord1, size_t o1_coord2, size_t o2_coord1, size_t o2_coord2);

    size_t get_size() const { return get_stride_z() * get_stride_y() * get_stride_x(); }

    void remove_patch(Patch patch);
    void remove_cells_at_boundary(size_t level);

    static void calculate_area_index(
            Obstacle *o1, Obstacle *o2,
            size_t *o1_coordinate, size_t *o2_coordinate,
            CoordinateAxis direction,
            bool start);
    static bool circular_constraints_x_direction(Obstacle *o1, Obstacle *o2);
    static bool circular_constraints_y_direction(Obstacle *o1, Obstacle *o2);
    static bool circular_constraints_z_direction(Obstacle *o1, Obstacle *o2);
};

#endif /* ARTSS_BOUNDARY_OBSTACLE_H_ */
