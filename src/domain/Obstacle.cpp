/// \file       Obstacle.cpp
/// \brief      Data class of obstacle object
/// \date       Oct 01, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Obstacle.h"
#include <algorithm>
#include <utility>
#include <vector>


real constexpr eps = 10E-6;

Obstacle::Obstacle(const Coordinate<real> &coords_start,
                   const Coordinate<real> &coords_end,
                   const std::string &name) :
        m_name(name),
        m_level(0),
        m_size_boundary() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    auto domain_data = DomainData::getInstance();
    for (CoordinateAxis axis: {CoordinateAxis::X, CoordinateAxis::Y, CoordinateAxis::Z}) {
        m_start[axis] = get_matching_index(coords_start[axis], domain_data->get_spacing(axis),
                                           domain_data->get_start_coord_PD(axis)) + 1;  // plus 1 for ghost cell
        m_end[axis] = get_matching_index(coords_end[axis], domain_data->get_spacing(axis),
                                         domain_data->get_start_coord_PD(axis));
    }
    init();
}

Obstacle::Obstacle(Coordinate<size_t> &coords_start, Coordinate<size_t> &coords_end,
                   size_t level,
                   const std::string &name) :
        m_name(name),
        m_start(coords_start), m_end(coords_end),
        m_level(level), m_size_boundary() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    init();
}


//======================================== Init ====================================================
// *************************************************************************************************
/// \brief  Initialize member variables (arrays)
/// \param  level Multigrid level
// *************************************************************************************************
void Obstacle::init() {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx(m_level);
    size_t Ny = domain_data->get_Ny(m_level);

    for (size_t axis = 0; axis < number_of_axes; axis++) {
        m_strides[axis] = m_end[axis] - m_start[axis] + 1;
    }

    m_size_obstacle_list = m_strides[CoordinateAxis::X] * m_strides[CoordinateAxis::Y] * m_strides[CoordinateAxis::Z];
    m_obstacle_list.resize(m_size_obstacle_list);

    m_size_boundary[FRONT] = m_strides[CoordinateAxis::Y] * m_strides[CoordinateAxis::X];
    m_size_boundary[BACK] = m_strides[CoordinateAxis::Y] * m_strides[CoordinateAxis::X];
    m_size_boundary[BOTTOM] = m_strides[CoordinateAxis::Z] * m_strides[CoordinateAxis::X];
    m_size_boundary[TOP] = m_strides[CoordinateAxis::Z] * m_strides[CoordinateAxis::X];
    m_size_boundary[LEFT] = m_strides[CoordinateAxis::Z] * m_strides[CoordinateAxis::Y];
    m_size_boundary[RIGHT] = m_strides[CoordinateAxis::Z] * m_strides[CoordinateAxis::Y];
    remove_cells_at_boundary();

    m_boundary.resize(number_of_patches);
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_boundary[patch].resize(m_size_boundary[patch]);
    }

    create_obstacle(Nx, Ny);

    control();
    print_details();
}

//===================================== Create obstacle ============================================
// *************************************************************************************************
/// \brief  Creates lists of indices of obstacle cells
// *************************************************************************************************
void Obstacle::create_obstacle(size_t Nx, size_t Ny) {
    size_t counter = 0;
    // fill obstacleList with corresponding indices
    for (size_t k = m_start[CoordinateAxis::Z]; k <= m_end[CoordinateAxis::Z]; ++k) {
        for (size_t j = m_start[CoordinateAxis::Y]; j <= m_end[CoordinateAxis::Y]; ++j) {
            for (size_t i = m_start[CoordinateAxis::X]; i <= m_end[CoordinateAxis::X]; ++i) {
                size_t idx = IX(i, j, k, Nx, Ny);
                m_obstacle_list[counter] = idx;
                counter++;
            }
        }
    }

#ifndef BENCHMARKING
    m_logger->debug("added {} cells to obstacle, array size: {}", counter, m_size_obstacle_list);
#endif

    // DETAILED OBSTACLE LISTS
    // FRONT and BACK of OBSTACLE
    // fill oFront list with front indices of obstacle and oBack list with back indices of obstacle
    if (m_size_boundary[FRONT] > 0) {
        counter = 0;
        for (size_t j = 0; j < m_strides[CoordinateAxis::Y]; ++j) {
            for (size_t i = 0; i < m_strides[CoordinateAxis::X]; ++i) {
                size_t idx_front = IX(i, j, 0, m_strides[CoordinateAxis::X], m_strides[CoordinateAxis::Y]);
                m_boundary[FRONT][counter++] = m_obstacle_list[idx_front];
            }
        }
#ifndef BENCHMARKING
        m_logger->debug("added {} front patch cells to obstacle, array size: {}", counter,
                        m_size_boundary[Patch::FRONT]);
#endif
    }
    if (m_size_boundary[BACK] > 0) {
        counter = 0;
        for (size_t j = 0; j < m_strides[CoordinateAxis::Y]; ++j) {
            for (size_t i = 0; i < m_strides[CoordinateAxis::X]; ++i) {
                size_t idx_back = IX(i, j, m_strides[CoordinateAxis::Z] - 1, m_strides[CoordinateAxis::X],
                                     m_strides[CoordinateAxis::Y]);
                m_boundary[BACK][counter++] = m_obstacle_list[idx_back];
            }
        }
#ifndef BENCHMARKING
        m_logger->debug("added {} back patch cells to obstacle, array size: {}", counter, m_size_boundary[Patch::BACK]);
#endif
    }

    // TOP and BOTTOM of OBSTACLE
    // fill m_boundary_patch_divided[TOP] list with top indices of obstacle and oBottom list with bottom indices
    // of obstacle
    if (m_size_boundary[BOTTOM] > 0) {
        counter = 0;
        for (size_t k = 0; k < m_strides[CoordinateAxis::Z]; ++k) {
            for (size_t i = 0; i < m_strides[CoordinateAxis::X]; ++i) {
                size_t idx_bottom = IX(i, 0, k, m_strides[CoordinateAxis::X], m_strides[CoordinateAxis::Y]);
                m_boundary[BOTTOM][counter++] = m_obstacle_list[idx_bottom];
            }
        }
#ifndef BENCHMARKING
        m_logger->debug("added {} bottom patch cells to obstacle, array size: {}", counter,
                        m_size_boundary[Patch::BOTTOM]);
#endif
    }
    if (m_size_boundary[TOP] > 0) {
        counter = 0;
        for (size_t k = 0; k < m_strides[CoordinateAxis::Z]; ++k) {
            for (size_t i = 0; i < m_strides[CoordinateAxis::X]; ++i) {
                size_t idx_top = IX(i, m_strides[CoordinateAxis::Y] - 1, k, m_strides[CoordinateAxis::X],
                                    m_strides[CoordinateAxis::Y]);
                m_boundary[TOP][counter++] = m_obstacle_list[idx_top];
            }
        }
#ifndef BENCHMARKING
        m_logger->debug("added {} top patch cells to obstacle, array size: {}", counter, m_size_boundary[Patch::TOP]);
#endif
    }

    // LEFT and RIGHT of OBSTACLE
    // fill oLeft list with left indices of obstacle and oRight list with right indices of obstacle
    if (m_size_boundary[LEFT] > 0) {
        counter = 0;
        for (size_t k = 0; k < m_strides[CoordinateAxis::Z]; ++k) {
            for (size_t j = 0; j < m_strides[CoordinateAxis::Y]; ++j) {
                size_t idx_left = IX(0, j, k, m_strides[CoordinateAxis::X], m_strides[CoordinateAxis::Y]);
                m_boundary[LEFT][counter++] = m_obstacle_list[idx_left];
            }
        }
#ifndef BENCHMARKING
        m_logger->debug("added {} left patch cells to obstacle, array size: {}", counter, m_size_boundary[Patch::LEFT]);
#endif
    }
    if (m_size_boundary[RIGHT] > 0) {
        counter = 0;
        for (size_t k = 0; k < m_strides[CoordinateAxis::Z]; ++k) {
            for (size_t j = 0; j < m_strides[CoordinateAxis::Y]; ++j) {
                size_t idx_right = IX(m_strides[CoordinateAxis::X] - 1, j, k, m_strides[CoordinateAxis::X],
                                      m_strides[CoordinateAxis::Y]);
                m_boundary[RIGHT][counter++] = m_obstacle_list[idx_right];
            }
        }
#ifndef BENCHMARKING
        m_logger->debug("added {} right patch cells to obstacle, array size: {}", counter,
                        m_size_boundary[Patch::RIGHT]);
#endif
    }
}

//======================================== Print ===================================================
// *************************************************************************************************
/// \brief  Print obstacle infos
// *************************************************************************************************
void Obstacle::print() const {
#ifndef BENCHMARKING
    m_logger->info("-- Obstacle {}", m_name);
    m_logger->info("\t strides (x y z): {} {} {}", m_strides[CoordinateAxis::X], m_strides[CoordinateAxis::Y],
                   m_strides[CoordinateAxis::Z]);
    m_logger->info("\t size of slices  (Front|Back Bottom|Top Left|Right): {}|{} {}|{} {}|{}",
                   m_size_boundary[FRONT], m_size_boundary[BACK],
                   m_size_boundary[BOTTOM], m_size_boundary[TOP],
                   m_size_boundary[LEFT], m_size_boundary[RIGHT]);
    m_logger->info("\t size of Obstacle: {}", m_size_obstacle_list);
    m_logger->info("\t coords (x y z): ({}|{}) ({}|{}) ({}|{})", m_start[CoordinateAxis::X], m_end[CoordinateAxis::X],
                   m_start[CoordinateAxis::Y], m_end[CoordinateAxis::Y],
                   m_start[CoordinateAxis::Z], m_end[CoordinateAxis::Z]);
#endif
}

//======================================== Print ===================================================
// *************************************************************************************************
/// \brief  Print detailed obstacle infos
// *************************************************************************************************
void Obstacle::print_details() {
#ifndef BENCHMARKING
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx(m_level);
    size_t Ny = domain_data->get_Ny(m_level);
    size_t coords_i, coords_j, coords_k;

    m_logger->debug("############### OBSTACLE {} ###############", m_name);
    m_logger->debug("level: {}", m_level);
    m_logger->debug("strides (x y z): {} {} {}", m_strides[CoordinateAxis::X], m_strides[CoordinateAxis::Y],
                    m_strides[CoordinateAxis::Z]);
    m_logger->debug("size of slices  (Front|Back Bottom|Top Left|Right): {}|{} {}|{} {}|{}",
                    m_size_boundary[FRONT], m_size_boundary[BACK],
                    m_size_boundary[BOTTOM], m_size_boundary[TOP],
                    m_size_boundary[LEFT], m_size_boundary[RIGHT]);
    m_logger->debug("size of Obstacle: {}", m_size_obstacle_list);
    m_logger->debug("coords (x y z): ({}|{}) ({}|{}) ({}|{})", m_start[CoordinateAxis::X], m_end[CoordinateAxis::X],
                    m_start[CoordinateAxis::Y], m_end[CoordinateAxis::Y],
                    m_start[CoordinateAxis::Z], m_end[CoordinateAxis::Z]);

    std::vector<size_t> coords;
    size_t size_front = (*get_size_boundary_list())[FRONT];
    if (size_front > 0) {
        m_logger->debug("Front: {} | {}", m_boundary[FRONT][0],
                        m_boundary[FRONT][size_front - 1]);

        coords_k = getCoordinateK(m_boundary[FRONT][0], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[FRONT][0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[FRONT][0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Front start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_boundary[FRONT][size_front - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[FRONT][size_front - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[FRONT][size_front - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Front end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Front size = 0");
    }

    size_t size_back = (*get_size_boundary_list())[BACK];
    if (size_back > 0) {
        m_logger->debug("Back: {} | {}", m_boundary[BACK][0], m_boundary[BACK][size_back - 1]);

        coords_k = getCoordinateK(m_boundary[BACK][0], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[BACK][0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[BACK][0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Back start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_boundary[BACK][size_back - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[BACK][size_back - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[BACK][size_back - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Back end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Back size = 0");
    }

    size_t size_top = (*get_size_boundary_list())[TOP];
    if (size_top > 0) {
        m_logger->debug("Top: {} | {}", m_boundary[TOP][0], m_boundary[TOP][size_top - 1]);

        coords_k = getCoordinateK(m_boundary[TOP][0], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[TOP][0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[TOP][0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Top start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_boundary[TOP][size_top - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[TOP][size_top - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[TOP][size_top - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Top end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Top size = 0");
    }

    size_t size_bottom = (*get_size_boundary_list())[BOTTOM];
    if (size_bottom > 0) {
        m_logger->debug("Bottom: {} | {}", m_boundary[BOTTOM][0], m_boundary[BOTTOM][size_bottom - 1]);

        coords_k = getCoordinateK(m_boundary[BOTTOM][0], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[BOTTOM][0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[BOTTOM][0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Bottom start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_boundary[BOTTOM][size_bottom - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[BOTTOM][size_bottom - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[BOTTOM][size_bottom - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Bottom end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Bottom size = 0");
    }

    size_t size_left = (*get_size_boundary_list())[LEFT];
    if (size_left > 0) {
        m_logger->debug("Left: {} | {}", m_boundary[LEFT][0], m_boundary[LEFT][size_left - 1]);

        coords_k = getCoordinateK(m_boundary[LEFT][0], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[LEFT][0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[LEFT][0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Left start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_boundary[LEFT][size_left - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[LEFT][size_left - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[LEFT][size_left - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Left end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Left size = 0");
    }

    size_t size_right = (*get_size_boundary_list())[RIGHT];
    if (size_right > 0) {
        m_logger->debug("Right: {} | {}", m_boundary[RIGHT][0], m_boundary[RIGHT][size_right - 1]);

        coords_k = getCoordinateK(m_boundary[RIGHT][0], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[RIGHT][0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[RIGHT][0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Right start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_boundary[RIGHT][size_right - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_boundary[RIGHT][size_right - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_boundary[RIGHT][size_right - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Right end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Right size = 0");
    }
    m_logger->debug("############### OBSTACLE END ###############");
#endif
}

//======================================== Control =================================================
// *************************************************************************************************
/// \brief  Units test emergency solution
// *************************************************************************************************
void Obstacle::control() {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_Nx(m_level);
    size_t Ny = domain_data->get_Ny(m_level);

    std::string message;
    for (size_t i = 1; i < m_size_obstacle_list; i++) {
        long int diff = static_cast<long int>(m_obstacle_list[i]) -
                        static_cast<long int>(m_obstacle_list[i - 1]);
        if (diff < 0) {
            message += "sorting error at index "
                       + std::to_string(i - 1) + "|" + std::to_string(i) + " with values "
                       + std::to_string(m_obstacle_list[i - 1]) + "|"
                       + std::to_string(m_obstacle_list[i]) + "\n";
        }
    }

    size_t start_index = IX(m_start[CoordinateAxis::X], m_start[CoordinateAxis::Y], m_start[CoordinateAxis::Z], Nx, Ny);
    size_t end_index = IX(m_end[CoordinateAxis::X], m_end[CoordinateAxis::Y], m_end[CoordinateAxis::Z], Nx, Ny);

    if (m_size_boundary[FRONT] > 0) {
        size_t front_end = IX(m_end[CoordinateAxis::X], m_end[CoordinateAxis::Y], m_start[CoordinateAxis::Z], Nx, Ny);
        if (start_index != m_boundary[FRONT][0] ||
            front_end != m_boundary[FRONT][m_size_boundary[FRONT] - 1]) {
            message += "first or last index of obstacle front list not correct ("
                       + std::to_string(start_index) + "|" + std::to_string(m_boundary[FRONT][0])
                       + ")(" + std::to_string(front_end) + "|"
                       + std::to_string(m_boundary[FRONT][m_size_boundary[FRONT] - 1]) + ")\n";
        }
    }
    if (m_size_boundary[BACK] > 0) {
        size_t back_start = IX(m_start[CoordinateAxis::X], m_start[CoordinateAxis::Y], m_end[CoordinateAxis::Z], Nx,
                               Ny);
        if (back_start != m_boundary[BACK][0] ||
            end_index != m_boundary[BACK][m_size_boundary[BACK] - 1]) {
            message += "first or last index of obstacle back list not correct ("
                       + std::to_string(back_start) + "|" + std::to_string(m_boundary[BACK][0])
                       + ")(" + std::to_string(end_index) + "|"
                       + std::to_string(m_boundary[BACK][m_size_boundary[BACK] - 1]) + ")\n";
        }
    }
    if (m_size_boundary[BOTTOM] > 0) {
        size_t bottom_end = IX(m_end[CoordinateAxis::X], m_start[CoordinateAxis::Y], m_end[CoordinateAxis::Z], Nx, Ny);
        if (start_index != m_boundary[BOTTOM][0] ||
            bottom_end != m_boundary[BOTTOM][m_size_boundary[BOTTOM] - 1]) {
            message += "first or last index of obstacle bottom list not correct ("
                       + std::to_string(start_index) + "|" + std::to_string(m_boundary[BOTTOM][0])
                       + ")(" + std::to_string(bottom_end) + "|"
                       + std::to_string(m_boundary[BOTTOM][m_size_boundary[BOTTOM] - 1]) + ")\n";
        }
    }
    if (m_size_boundary[TOP] > 0) {
        size_t top_start = IX(m_start[CoordinateAxis::X], m_end[CoordinateAxis::Y], m_start[CoordinateAxis::Z], Nx, Ny);
        if (top_start != m_boundary[TOP][0] ||
            end_index != m_boundary[TOP][m_size_boundary[TOP] - 1]) {
            message += "first or last index of obstacle top list not correct ("
                       + std::to_string(top_start) + "|" + std::to_string(m_boundary[TOP][0])
                       + ")(" + std::to_string(end_index) + "|"
                       + std::to_string(m_boundary[TOP][m_size_boundary[TOP] - 1]) + ")\n";
        }
    }
    if (m_size_boundary[LEFT] > 0) {
        size_t left_end = IX(m_start[CoordinateAxis::X], m_end[CoordinateAxis::Y], m_end[CoordinateAxis::Z], Nx, Ny);
        if (start_index != m_boundary[LEFT][0] ||
            left_end != m_boundary[LEFT][m_size_boundary[LEFT] - 1]) {
            message += "first or last index of obstacle left list not correct ("
                       + std::to_string(start_index) + "|" + std::to_string(m_boundary[LEFT][0])
                       + ")(" + std::to_string(left_end) + "|"
                       + std::to_string(m_boundary[LEFT][m_size_boundary[LEFT] - 1]) + ")\n";
        }
    }
    if (m_size_boundary[RIGHT] > 0) {
        size_t right_start = IX(m_end[CoordinateAxis::X], m_start[CoordinateAxis::Y], m_start[CoordinateAxis::Z], Nx,
                                Ny);
        if (right_start != m_boundary[RIGHT][0] ||
            end_index != m_boundary[RIGHT][m_size_boundary[RIGHT] - 1]) {
            message += "first or last index of obstacle right list not correct ("
                       + std::to_string(right_start) + "|" + std::to_string(m_boundary[RIGHT][0])
                       + ")(" + std::to_string(end_index) + "|"
                       + std::to_string(m_boundary[RIGHT][m_size_boundary[RIGHT] - 1]) + ")\n";
        }
    }
    if (!message.empty()) {
        message = "################ OBSTACLE CONTROL ################\n-- name "
                  + m_name + "\n-- level " + std::to_string(m_level) + "\n" + message
                  + "---------------- OBSTACLE CONTROL END ----------------";
#ifndef BENCHMARKING
        m_logger->warn(message);
#endif
    }
}

//======================================== Match grid ==============================================
// *************************************************************************************************
/// \brief  Snaps value to grid discretisation
/// \param  obstacle_coordinate Coordinate<size_t> of obstacle
/// \param  spacing dx/dy/dz
/// \param  start_coordinate X1/Y1/Z1
/// \return real Calculated index on grid
// *************************************************************************************************
int Obstacle::get_matching_index(real obstacle_coordinate, real spacing, real start_coordinate) {
    return static_cast<int>(round((-start_coordinate + obstacle_coordinate) / spacing));
}

//======================================== Remove cells at boundary ================================
// *************************************************************************************************
/// \brief  Remove obstacle patch facing the boundary
/// \param  level Multigrid level
// *************************************************************************************************
void Obstacle::remove_cells_at_boundary() {
    auto domain_data = DomainData::getInstance();
    if (m_start[CoordinateAxis::Z] <= domain_data->get_index_z1(m_level)) {
        m_size_boundary[FRONT] = 0;
    }
    if (m_end[CoordinateAxis::Z] >= domain_data->get_index_z2(m_level)) {
        m_size_boundary[BACK] = 0;
    }
    if (m_start[CoordinateAxis::Y] <= domain_data->get_index_y1(m_level)) {
        m_size_boundary[BOTTOM] = 0;
    }
    if (m_end[CoordinateAxis::Y] >= domain_data->get_index_y2(m_level)) {
        m_size_boundary[TOP] = 0;
    }
    if (m_start[CoordinateAxis::X] <= domain_data->get_index_x1(m_level)) {
        m_size_boundary[LEFT] = 0;
    }
    if (m_end[CoordinateAxis::X] >= domain_data->get_index_x2(m_level)) {
        m_size_boundary[RIGHT] = 0;
    }
}

//======================================== has overlap =============================================
// *************************************************************************************************
/// checks, if two intervals defined by two coordinates are overlapping
/// \param o1_coord1 starting coordinate of obstacle 1
/// \param o1_coord2 ending coordinate of obstacle 1
/// \param o2_coord1 starting coordinate of obstacle 2
/// \param o2_coord2 starting coordinate of obstacle 2
/// \return true if its overlapping, false otherwise
// *************************************************************************************************
bool Obstacle::has_overlap(size_t o1_coord1, size_t o1_coord2, size_t o2_coord1, size_t o2_coord2) {
    return o1_coord1 <= o2_coord2 && o1_coord2 >= o2_coord1;
}

void Obstacle::replace_patch(size_t *indices, size_t size, Patch p) {
    m_boundary[p].clear();
    m_boundary[p].resize(size);
    for (size_t i = 0; i < size; i++) {
        m_boundary[p][i] = indices[i];
    }
    m_size_boundary[p] = size;
}

//======================================== has overlap =============================================
// *************************************************************************************************
/// \brief calculate indices of area to be excluded. o1_coordinate == o2_coordinate only if the length
/// of the patch of both obstacles are the same
/// \param o1 Obstacle 1
/// \param o2 Obstacle 2
/// \param o1_coordinate calculated coordinate of obstacle 1
/// \param o2_coordinate calculated coordinate of obstacle 2
/// \param coordinate_axis X/Y/Z axis
/// \param start true = (i/j/k)1 or false = (i/j/k)2
// *************************************************************************************************
void Obstacle::calculate_area_index(
        Obstacle &o1, Obstacle &o2,
        size_t *o1_coordinate, size_t *o2_coordinate,
        CoordinateAxis coordinate_axis,
        bool start) {
    Coordinate<size_t> coords_start_o1 = o1.get_start_coordinates();
    Coordinate<size_t> coords_end_o1 = o1.get_end_coordinates();
    Coordinate<size_t> coords_start_o2 = o2.get_start_coordinates();
    Coordinate<size_t> coords_end_o2 = o2.get_end_coordinates();

    if (start) {
        *o1_coordinate = coords_start_o1[coordinate_axis];
        *o2_coordinate = coords_start_o2[coordinate_axis];
        if (coords_start_o1[coordinate_axis] > coords_start_o2[coordinate_axis]) {
            // do not remove inner edge, can be accessed by SL Advection Solver
            *o2_coordinate = coords_start_o1[coordinate_axis] + 1;
        } else if (coords_start_o1[coordinate_axis] < coords_start_o2[coordinate_axis]) {
            // do not remove inner edge, can be accessed by SL Advection Solver
            *o1_coordinate = coords_start_o2[coordinate_axis] + 1;
        }
    } else {
        *o1_coordinate = coords_end_o1[coordinate_axis];
        *o2_coordinate = coords_end_o2[coordinate_axis];
        if (coords_end_o1[coordinate_axis] < coords_end_o2[coordinate_axis]) {
            // do not remove inner edge, can be accessed by SL Advection Solver
            *o2_coordinate = coords_end_o1[coordinate_axis] - 1;
        } else if (coords_end_o1[coordinate_axis] > coords_end_o2[coordinate_axis]) {
            // do not remove inner edge, can be accessed by SL Advection Solver
            *o1_coordinate = coords_end_o2[coordinate_axis] - 1;
        }
    }
}

//======================================== remove circular constraints =============================
// *************************************************************************************************
/// \brief removes boundary cells of obstacle o1 and obstacle o2 which are facing each other and
/// lead to circular constraints. If the two obstacles are not next to each other nothing happens.
/// \details In worst case these constraints are overwriting important boundary cells. To prevent
/// this, cells which does not change their own value because of the circular constraints will be
/// removed. There are two different cases:
/// 1. both patches are of the same size and facing each other: all cells will be removed
/// 2. both patches are not facing each other completely: if obstacle o1 is longer than obstacle
///    o2 the last shared cell will be only removed in o2. In o2 this cell is a edge (or corner)
///    cell therefore the cell value can change through the other neighbouring patches. For this
///    reason the cell of o1 may also change its value.
/// \param o1 Obstacle 1
/// \param o2 Obstacle 2
/// \return true, if cells were removed otherwise false
// *************************************************************************************************
bool Obstacle::remove_circular_constraints(Obstacle &o1, Obstacle &o2) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger(typeid(Obstacle).name());
#endif
    bool any_overlap = false;
    for (size_t axis = 0; axis < number_of_axes; axis++) {
        auto coordinate_axis = CoordinateAxis(axis);
#ifndef BENCHMARKING
        logger->debug("neighbouring obstacles ! comparing {} {} with {} {} in {} direction",
                      o1.get_name(), o1.get_end_coordinate(coordinate_axis),
                      o2.get_name(), o2.get_start_coordinate(coordinate_axis),
                      Mapping::get_axis_name(coordinate_axis));
#endif
        bool overlap;
        if (o1.get_end_coordinate(coordinate_axis) + 1 == o2.get_start_coordinate(coordinate_axis)) {
            overlap = circular_constraints(o2, o1, coordinate_axis);
        } else {
            overlap = circular_constraints(o1, o2, coordinate_axis);
        }
        any_overlap = any_overlap || overlap;
    }
#ifndef BENCHMARKING
    if (any_overlap) {
        logger->debug("neighbouring obstacles ! {} is next to {}", o1.get_name(), o2.get_name());
    }
#endif
    return any_overlap;
}

//============================ circular constraints =================================
// *************************************************************************************************
/// \brief removes circular constraints in specified direction if two obstacles are next to each other
/// \param o1 obstacle 1
/// \param o2 obstacle 2
// *************************************************************************************************
bool Obstacle::circular_constraints(Obstacle &o1, Obstacle &o2, CoordinateAxis coordinate_axis) {
    const Coordinate<size_t> &coords_start_o1 = (o1.get_start_coordinates());
    const Coordinate<size_t> &coords_end_o1 = (o1.get_end_coordinates());
    const Coordinate<size_t> &coords_start_o2 = (o2.get_start_coordinates());
    const Coordinate<size_t> &coords_end_o2 = (o2.get_end_coordinates());

    bool overlap = false;

    auto domain_data = DomainData::getInstance();
    auto Nx = domain_data->get_Nx();
    auto Ny = domain_data->get_Ny();

#ifndef BENCHMARKING
    auto logger = Utility::create_logger(typeid(Obstacle).name());
    logger->debug("neighbouring obstacles ! comparing {} {} with {} {} in {} direction",
                  o1.get_name(), coords_start_o1[coordinate_axis],
                  o2.get_name(), coords_end_o2[coordinate_axis],
                  Mapping::get_axis_name(coordinate_axis));
#endif
    if (coords_start_o1[coordinate_axis] - 1 == coords_end_o2[coordinate_axis]) {
        auto other_axes = new CoordinateAxis[2];
        if (coordinate_axis == CoordinateAxis::X) {
            other_axes[0] = CoordinateAxis::Y;
            other_axes[1] = CoordinateAxis::Z;
        } else if (coordinate_axis == CoordinateAxis::Y) {
            other_axes[0] = CoordinateAxis::X;
            other_axes[1] = CoordinateAxis::Z;
        } else if (coordinate_axis == CoordinateAxis::Z) {
            other_axes[0] = CoordinateAxis::X;
            other_axes[1] = CoordinateAxis::Y;
        }
        // for constraints in x-direction, check for overlapping in y-direction
        bool overlap1 = has_overlap(
                coords_start_o1[other_axes[0]], coords_end_o1[other_axes[0]],
                coords_start_o2[other_axes[0]], coords_end_o2[other_axes[0]]);
        // for constraints in x-direction, check for overlapping in z-direction
        bool overlap2 = has_overlap(
                coords_start_o1[other_axes[1]], coords_end_o1[other_axes[1]],
                coords_start_o2[other_axes[1]], coords_end_o2[other_axes[1]]);
        if (overlap1 && overlap2) {
            Coordinate<size_t> tmp;
            auto o1_patch = Mapping::to_patch(coordinate_axis, true);
            auto o2_patch = Mapping::to_patch(coordinate_axis, false);
#ifndef BENCHMARKING
            // for constraints in x-direction: front patch of o1 and back patch of o2
            logger->debug(
                    "neighbouring obstacles ! obstacles are next to each other. Working on '{}' {} side and on '{}' {} side",
                    o1.get_name(), Mapping::get_patch_name(o1_patch),
                    o2.get_name(), Mapping::get_patch_name(o2_patch));
#endif
            overlap = true;
            // calculate coordinates of area which should be removed
            // the area is for both obstacle the same only if there are equally long

            Coordinate<size_t> o1_remove_start;
            Coordinate<size_t> o1_remove_end;
            Coordinate<size_t> o2_remove_start;
            Coordinate<size_t> o2_remove_end;

            o1_remove_start[coordinate_axis] = coords_start_o1[coordinate_axis];
            o1_remove_end[coordinate_axis] = coords_start_o1[coordinate_axis];
            o2_remove_start[coordinate_axis] = coords_end_o2[coordinate_axis];
            o2_remove_end[coordinate_axis] = coords_end_o2[coordinate_axis];

            size_t o1_y1;
            size_t o2_y1;
            Obstacle::calculate_area_index(o1, o2, &o1_y1, &o2_y1, other_axes[0], true);
            o1_remove_start[other_axes[0]] = o1_y1;
            o2_remove_start[other_axes[0]] = o2_y1;

            size_t o1_y2;
            size_t o2_y2;
            Obstacle::calculate_area_index(o1, o2, &o1_y2, &o2_y2, other_axes[0], false);
            o1_remove_end[other_axes[0]] = o1_y2;
            o2_remove_end[other_axes[0]] = o2_y2;

            size_t o1_z1;
            size_t o2_z1;
            Obstacle::calculate_area_index(o1, o2, &o1_z1, &o2_z1, other_axes[1], true);
            o1_remove_start[other_axes[1]] = o1_z1;
            o2_remove_start[other_axes[1]] = o2_z1;

            size_t o1_z2;
            size_t o2_z2;
            Obstacle::calculate_area_index(o1, o2, &o1_z2, &o2_z2, other_axes[1], false);
            o1_remove_end[other_axes[1]] = o1_z2;
            o2_remove_end[other_axes[1]] = o2_z2;

#ifndef BENCHMARKING
            logger->debug("neighbouring obstacles ! removing indices in area ({}|{}) ({}|{}) ({}|{}) for {}",
                          o1_remove_start[CoordinateAxis::X], o1_remove_end[CoordinateAxis::X],
                          o1_remove_start[CoordinateAxis::Y], o1_remove_end[CoordinateAxis::Y],
                          o1_remove_start[CoordinateAxis::Z], o1_remove_end[CoordinateAxis::Z],
                          o1.get_name());
            logger->debug("neighbouring obstacles ! removing indices in area ({}|{}) ({}|{}) ({}|{}) for {}",
                          o2_remove_start[CoordinateAxis::X], o2_remove_end[CoordinateAxis::X],
                          o2_remove_start[CoordinateAxis::Y], o2_remove_end[CoordinateAxis::Y],
                          o2_remove_start[CoordinateAxis::Z], o2_remove_end[CoordinateAxis::Z],
                          o2.get_name());
#endif

            std::vector<size_t> o1_new;
            o1_new.reserve((*o1.get_size_boundary_list())[o1_patch]);
            std::vector<size_t> o2_new;
            o2_new.reserve((*o2.get_size_boundary_list())[o2_patch]);

            // add all cells which are smaller than the smallest removing index of o1
            size_t o1_counter_old = 0;
            size_t o1_smallest_removing_index = o1_remove_start.get_index(Nx, Ny);
            size_t o1_current_index = o1.get_boundary_list()[o1_patch][o1_counter_old];
            while (o1_current_index < o1_smallest_removing_index) {
                o1_new.push_back(o1_current_index);
                o1_counter_old++;
                o1_current_index = o1.get_boundary_list()[o1_patch][o1_counter_old];
            }
            size_t o1_new_size = o1_counter_old;

            // add all cells which are smaller than the smallest removing index of o2
            size_t o2_counter_old = 0;
            size_t o2_smallest_removing_index = o2_remove_start.get_index(Nx, Ny);
            size_t o2_current_index = o2.get_boundary_list()[o2_patch][o2_counter_old];
            while (o2_current_index < o2_smallest_removing_index) {
                o2_new.push_back(o2_current_index);
                o2_counter_old++;
                o2_current_index = o2.get_boundary_list()[o2_patch][o2_counter_old];
            }
            size_t o2_new_size = o2_counter_old;

            size_t o1_current_axis_0 = o1_remove_start[other_axes[0]];
            size_t o1_current_axis_1 = o1_remove_start[other_axes[1]];
            size_t o1_removing_index = o1_smallest_removing_index;
            bool o1_end = false;

            size_t o2_current_axis_0 = o2_remove_start[other_axes[0]];
            size_t o2_current_axis_1 = o2_remove_start[other_axes[1]];
            size_t o2_removing_index = o2_smallest_removing_index;
            bool o2_end = false;
            for (; o1_counter_old < (*o1.get_size_boundary_list())[o1_patch]
                   && o2_counter_old < (*o2.get_size_boundary_list())[o2_patch]
                   && !o1_end && !o2_end;
                   o1_counter_old++, o2_counter_old++) {
                o1_current_index = o1.get_boundary_list()[o1_patch][o1_counter_old];
                o2_current_index = o2.get_boundary_list()[o2_patch][o2_counter_old];
                if (o1_current_index != o1_removing_index) {
                    o1_new.push_back(o1_current_index);
                    o1_new_size++;
                } else {
                    o1_current_axis_0++;
                    if (o1_current_axis_0 > o1_remove_end[other_axes[0]]) {
                        o1_current_axis_0 = o1_remove_start[other_axes[0]];
                        o1_current_axis_1++;
                        if (o1_current_axis_1 > o1_remove_end[other_axes[1]]) {
                            o1_end = true;
                        }
                    }
                    tmp[other_axes[0]] = o1_current_axis_0;
                    tmp[other_axes[1]] = o1_current_axis_1;
                    tmp[coordinate_axis] = o1_remove_start[coordinate_axis];
                    o1_removing_index = tmp.get_index(Nx, Ny);
                }
                if (o2_current_index != o2_removing_index) {
                    o2_new.push_back(o2_current_index);
                    o2_new_size++;
                } else {
                    o2_current_axis_0++;
                    if (o2_current_axis_0 > o2_remove_end[other_axes[0]]) {
                        o2_current_axis_0 = o2_remove_start[other_axes[0]];
                        o2_current_axis_1++;
                        if (o2_current_axis_1 > o2_remove_end[other_axes[1]]) {
                            o2_end = true;
                        }
                    }
                    tmp[other_axes[0]] = o2_current_axis_0;
                    tmp[other_axes[1]] = o2_current_axis_1;
                    tmp[coordinate_axis] = o2_remove_start[coordinate_axis];
                    o2_removing_index = tmp.get_index(Nx, Ny);
                }
            }

            if (!o1_end) {
                for (; o1_counter_old < (*o1.get_size_boundary_list())[o1_patch] &&
                       o1_current_axis_1 <= o1_remove_end[other_axes[1]]; o1_counter_old++) {
                    o1_current_index = o1.get_boundary_list()[o1_patch][o1_counter_old];
                    if (o1_current_index != o1_removing_index) {
                        o1_new.push_back(o1_current_index);
                        o1_new_size++;
                    } else {
                        o1_current_axis_0++;
                        if (o1_current_axis_0 > o1_remove_end[other_axes[0]]) {
                            o1_current_axis_0 = o1_remove_start[other_axes[0]];
                            o1_current_axis_1++;
                        }
                        tmp[other_axes[0]] = o1_current_axis_0;
                        tmp[other_axes[1]] = o1_current_axis_1;
                        tmp[coordinate_axis] = o1_remove_start[coordinate_axis];
                        o1_removing_index = tmp.get_index(Nx, Ny);
                    }
                }
            }

            if (!o2_end) {
                for (; o2_counter_old < (*o2.get_size_boundary_list())[o2_patch] &&
                       o2_current_axis_1 <= o2_remove_end[other_axes[1]]; o2_counter_old++) {
                    o2_current_index = o2.get_boundary_list()[o2_patch][o2_counter_old];
                    if (o2_current_index != o2_removing_index) {
                        o2_new.push_back(o2_current_index);
                        o2_new_size++;
                    } else {
                        o2_current_axis_0++;
                        if (o2_current_axis_0 > o2_remove_end[other_axes[0]]) {
                            o2_current_axis_0 = o2_remove_start[other_axes[0]];
                            o2_current_axis_1++;
                        }
                        tmp[other_axes[0]] = o2_current_axis_0;
                        tmp[other_axes[1]] = o2_current_axis_1;
                        tmp[coordinate_axis] = o2_remove_start[coordinate_axis];
                        o2_removing_index = tmp.get_index(Nx, Ny);
                    }
                }
            }

            for (; o1_counter_old < (*o1.get_size_boundary_list())[o1_patch]; o1_counter_old++) {
                o1_new.push_back(o1.get_boundary_list()[o1_patch][o1_counter_old]);
                o1_new_size++;
            }
            o1_new.resize(o1_new_size);

            size_t o1_diff_target = (o1_remove_end[other_axes[0]] - o1_remove_start[other_axes[0]] + 1) *
                                    (o1_remove_end[other_axes[1]] - o1_remove_start[other_axes[1]] + 1);
            size_t o1_diff_actual = (*o1.get_size_boundary_list())[o1_patch] - o1_new_size;
#ifndef BENCHMARKING
            logger->debug("neighbouring obstacles ! new size of obstacle '{}' {} patch: {} -> {} ({}|{})",
                          o1.get_name(), Mapping::get_patch_name(o1_patch),
                          (*o1.get_size_boundary_list())[o1_patch], o1_new_size,
                          o1_diff_target, o1_diff_actual);
#endif
            auto o1_new_data = new size_t[o1_new_size];
            std::copy(o1_new.begin(), o1_new.end(), o1_new_data);
            o1.replace_patch(o1_new_data, o1_new_size, o1_patch);

            for (; o2_counter_old < (*o2.get_size_boundary_list())[o2_patch]; o2_counter_old++) {
                o2_new.push_back(o2.get_boundary_list()[o2_patch][o2_counter_old]);
                o2_new_size++;
            }
            o2_new.resize(o2_new_size);

            size_t o2_diff_target = (o2_remove_end[other_axes[0]] - o2_remove_start[other_axes[0]] + 1) *
                                    (o2_remove_end[other_axes[1]] - o2_remove_start[other_axes[1]] + 1);
            size_t o2_diff_actual = (*o2.get_size_boundary_list())[o2_patch] - o2_new_size;
#ifndef BENCHMARKING
            logger->debug("neighbouring obstacles ! new size of obstacle '{}' {} patch: {} -> {} ({}|{})",
                          o2.get_name(), Mapping::get_patch_name(o2_patch),
                          (*o2.get_size_boundary_list())[o2_patch], o2_new_size,
                          o2_diff_target, o2_diff_actual);
#endif
            auto o2_new_data = new size_t[o2_new_size];
            std::copy(o2_new.begin(), o2_new.end(), o2_new_data);
            o2.replace_patch(o2_new_data, o2_new_size, o2_patch);
        }
        delete[] other_axes;
    }
    return overlap;
}

//======================================== is corner cell ==========================================
// *************************************************************************************************
/// \brief return whether cell is a corner cell
/// \param coord Coordinate<size_t> triple
/// \return true if cell is a corner cell, otherwise false
// *************************************************************************************************
bool Obstacle::is_corner_cell(const Coordinate<size_t> &coord) const {
    bool on_x = (coord[CoordinateAxis::X] == m_start[CoordinateAxis::X] ||
                 coord[CoordinateAxis::X] == m_end[CoordinateAxis::X]);
    bool on_y = (coord[CoordinateAxis::Y] == m_start[CoordinateAxis::Y] ||
                 coord[CoordinateAxis::Y] == m_end[CoordinateAxis::Y]);
    bool on_z = (coord[CoordinateAxis::Z] == m_start[CoordinateAxis::Z] ||
                 coord[CoordinateAxis::Z] == m_end[CoordinateAxis::Z]);
    return on_x && on_y && on_z;
}

//======================================== is edge cell ============================================
// *************************************************************************************************
/// \brief return whether cell is a edge cell
/// \param coord Coordinate<size_t> triple
/// \return true if cell is a edge cell, otherwise false
// *************************************************************************************************
bool Obstacle::is_edge_cell(const Coordinate<size_t> &coord) const {
    bool on_x = (coord[CoordinateAxis::X] == m_start[CoordinateAxis::X] ||
                 coord[CoordinateAxis::X] == m_end[CoordinateAxis::X]);
    bool on_y = (coord[CoordinateAxis::Y] == m_start[CoordinateAxis::Y] ||
                 coord[CoordinateAxis::Y] == m_end[CoordinateAxis::Y]);
    bool on_z = (coord[CoordinateAxis::Z] == m_start[CoordinateAxis::Z] ||
                 coord[CoordinateAxis::Z] == m_end[CoordinateAxis::Z]);
    return (on_x && on_y) || (on_y || on_z) || (on_x && on_z);
}

bool Obstacle::has_overlap(size_t i1, size_t i2, size_t j1, size_t j2, size_t k1, size_t k2) const {
    return has_overlap(m_start[CoordinateAxis::X], m_end[CoordinateAxis::X], i1, i2) &&
           has_overlap(m_start[CoordinateAxis::Y], m_end[CoordinateAxis::Y], j1, j2) &&
           has_overlap(m_start[CoordinateAxis::Z], m_end[CoordinateAxis::Z], k1, k2);
}

// ================================ intersection ===================================================
// *************************************************************************************************
/// \brief Calculates whether the obstacle is between the two points or not.
/// \details Done by calculation the intersection point between area spanned by the individual
/// patches and the line spanned by the two given points. Algorithm is optimised due to the patches
/// being axis aligned.
/// \param start start point
/// \param end end point
/// \return bool, whether the obstacle intersects with the line from start to end
// *************************************************************************************************
bool Obstacle::intersection(const Coordinate<size_t> &start, const Coordinate<size_t> &end) const {
    if (is_obstacle_cell(end)) {
        return true;
    }
    Coordinate<real> tmp_start_line = {
            static_cast<real>(start[CoordinateAxis::X]) + 0.5,
            static_cast<real>(start[CoordinateAxis::Y]) + 0.5,
            static_cast<real>(start[CoordinateAxis::Z]) + 0.5};
    Coordinate<real> tmp_end_line = {
            static_cast<real>(end[CoordinateAxis::X] ) + 0.5,
            static_cast<real>(end[CoordinateAxis::Y] ) + 0.5,
            static_cast<real>(end[CoordinateAxis::Z] ) + 0.5};
    Coordinate<real> direction_vector;
    for (size_t axis = 0; axis < number_of_axes; axis++) {
        direction_vector[axis] = tmp_end_line[axis] - tmp_start_line[axis];
    }
    // line equation : (tmp_start_line) + a * (direction_vector)

    for (CoordinateAxis coordinate_axis: {CoordinateAxis::X, CoordinateAxis::Y, CoordinateAxis::Z}) {
        CoordinateAxis other_axes[2];
        if (coordinate_axis == CoordinateAxis::X) {
            other_axes[0] = CoordinateAxis::Y;
            other_axes[1] = CoordinateAxis::Z;
        } else if (coordinate_axis == CoordinateAxis::Y) {
            other_axes[0] = CoordinateAxis::X;
            other_axes[1] = CoordinateAxis::Z;
        } else if (coordinate_axis == CoordinateAxis::Z) {
            other_axes[0] = CoordinateAxis::X;
            other_axes[1] = CoordinateAxis::Y;
        }
        Coordinate<real> tmp_start_area = {
                static_cast<real>(m_start[CoordinateAxis::X]) + 0.5,
                static_cast<real>(m_start[CoordinateAxis::Y]) + 0.5,
                static_cast<real>(m_start[CoordinateAxis::Z]) + 0.5};
        Coordinate<real> tmp_end_area = {
                static_cast<real>(m_end[CoordinateAxis::X]) + 0.5,
                static_cast<real>(m_end[CoordinateAxis::Y]) + 0.5,
                static_cast<real>(m_end[CoordinateAxis::Z]) + 0.5};
        for (auto patch: Mapping::get_patches(CoordinateAxis(coordinate_axis))) {
            if (patch % 2 == 1) {  // equals patch == BACK || patch == RIGHT || patch == TOP
                // set point of area.
                // for FRONT,LEFT,TOP it is m_start
                // for bACK, RIGHT, TOP it is m_end
                std::swap(tmp_start_area, tmp_end_area);
            }

            // calculate parameter form of area by using the start point and the
            // respective edges of the area, resulting in orthogonal direction vectors
            // e.g. Patch FRONT has the starting point in the lower left corner,
            // one direction vectors goes to the x-direction and the other
            // direction vector goes to the y-direction.

            // normal vector of area equation :
            Coordinate<real> normal_vector = {0, 0, 0};
            normal_vector[coordinate_axis] =
                    (tmp_end_area[other_axes[0]] - tmp_start_area[other_axes[0]]) *
                    (tmp_end_area[other_axes[1]] - tmp_start_area[other_axes[1]]);
            real res_dot = dot(normal_vector, direction_vector);
            if (fabs(res_dot) > eps) {  // else area and line are parallel
                // intersection
                real lambda = (tmp_start_area[coordinate_axis] - tmp_start_line[coordinate_axis]) /
                              direction_vector[coordinate_axis];
                if (lambda <= 1 && lambda >= 0) {
                    Coordinate<real> intersection = direction_vector;
                    intersection *= lambda;
                    intersection += tmp_start_line;
                    for (size_t axis = 0; axis < number_of_axes; axis++) {
                        intersection[axis] -= 0.5;
                    }
                    if (is_obstacle_cell(intersection)) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}
