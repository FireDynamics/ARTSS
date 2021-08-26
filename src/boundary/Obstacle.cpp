/// \file       Obstacle.cpp
/// \brief      Data class of obstacle object
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Obstacle.h"
#include <algorithm>
#include <utility>
#include <vector>


Obstacle::Obstacle(
        real x1, real x2, real y1,
        real y2, real z1, real z2,
        const std::string &name) :
    m_domain(*(Domain::getInstance())), m_name(name) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    real dx = m_domain.get_dx();
    real dy = m_domain.get_dy();
    real dz = m_domain.get_dz();

    real X1 = m_domain.get_X1();
    real Y1 = m_domain.get_Y1();
    real Z1 = m_domain.get_Z1();

    m_i1 = get_matching_index(x1, dx, X1) + 1;  // plus 1 for ghost cell
    m_j1 = get_matching_index(y1, dy, Y1) + 1;
    m_k1 = get_matching_index(z1, dz, Z1) + 1;

    m_i2 = get_matching_index(x2, dx, X1);
    m_j2 = get_matching_index(y2, dy, Y1);
    m_k2 = get_matching_index(z2, dz, Z1);

    init(0);
}

Obstacle::Obstacle(
        size_t coords_i1, size_t coords_j1, size_t coords_k1,
        size_t coords_i2, size_t coords_j2, size_t coords_k2,
        size_t level,
        const std::string& name,
        const Domain &domain) :
    m_i1(coords_i1), m_j1(coords_j1), m_k1(coords_k1),
    m_i2(coords_i2), m_j2(coords_j2), m_k2(coords_k2),
    m_level(level), m_name(name), m_domain(domain) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    init(m_level);
}

#ifndef BENCHMARKING
Obstacle::Obstacle(
        size_t coords_i1, size_t coords_j1, size_t coords_k1,
        size_t coords_i2, size_t coords_j2, size_t coords_k2,
        size_t level,
        std::shared_ptr<spdlog::logger> logger,
        const Domain &domain) :
    m_i1(coords_i1), m_j1(coords_j1), m_k1(coords_k1),
    m_i2(coords_i2), m_j2(coords_j2), m_k2(coords_k2),
    m_level(level),
    m_domain(domain), m_logger(logger) {
    init(m_level);
}
#endif

//======================================== Init ====================================================
// *************************************************************************************************
/// \brief  Initialize member variables (arrays)
/// \param  level Multigrid level
// *************************************************************************************************
void Obstacle::init(size_t level) {
    size_t Nx = m_domain.get_Nx(level);
    size_t Ny = m_domain.get_Ny(level);

    size_t strideX = get_stride_x();
    size_t strideY = get_stride_y();
    size_t strideZ = get_stride_z();

    m_size_obstacle_list = strideX * strideY * strideZ;
    m_obstacle_list = new size_t[m_size_obstacle_list];

    m_size_obstacle_front = strideY * strideX;
    m_size_obstacle_back = strideY * strideX;
    m_size_obstacle_bottom = strideZ * strideX;
    m_size_obstacle_top = strideZ * strideX;
    m_size_obstacle_left = strideZ * strideY;
    m_size_obstacle_right = strideZ * strideY;
    remove_cells_at_boundary(level);

    m_obstacle_front = new size_t[m_size_obstacle_front];
    m_obstacle_back = new size_t[m_size_obstacle_back];

    m_obstacle_top = new size_t[m_size_obstacle_top];
    m_obstacle_bottom = new size_t[m_size_obstacle_bottom];

    m_obstacle_left = new size_t[m_size_obstacle_left];
    m_obstacle_right = new size_t[m_size_obstacle_right];

    create_obstacle(Nx, Ny);

    control();
    print_details();
}

Obstacle::~Obstacle() {
    delete (m_obstacle_list);
    delete (m_obstacle_front);
    delete (m_obstacle_back);
    delete (m_obstacle_top);
    delete (m_obstacle_bottom);
    delete (m_obstacle_left);
    delete (m_obstacle_right);
}

//===================================== Create obstacle ============================================
// *************************************************************************************************
/// \brief  Creates lists of indices of obstacle cells
// *************************************************************************************************
void Obstacle::create_obstacle(size_t Nx, size_t Ny) {
    size_t strideX = get_stride_x();
    size_t strideY = get_stride_y();
    size_t strideZ = get_stride_z();

    size_t counter = 0;
    // fill obstacleList with corresponding indices
    for (size_t k = m_k1; k <= m_k2; ++k) {
        for (size_t j = m_j1; j <= m_j2; ++j) {
            for (size_t i = m_i1; i <= m_i2; ++i) {
                size_t idx = IX(i, j, k, Nx, Ny);
                *(m_obstacle_list + counter) = idx;
                counter++;
            }
        }
    }

    // DETAILED OBSTACLE LISTS
    // FRONT and BACK of OBSTACLE
    // fill oFront list with front indices of obstacle and oBack list with back indices of obstacle
    if (m_size_obstacle_front > 0) {
        for (size_t j = 0; j < strideY; ++j) {
            for (size_t i = 0; i < strideX; ++i) {
                size_t index = i + strideX * j;
                size_t idx_front = IX(i, j, 0, strideX, strideY);
                *(m_obstacle_front + index) = m_obstacle_list[idx_front];
            }
        }
    }
    if (m_size_obstacle_back > 0) {
        for (size_t j = 0; j < strideY; ++j) {
            for (size_t i = 0; i < strideX; ++i) {
                size_t index = i + strideX * j;
                size_t idx_back = IX(i, j, strideZ - 1, strideX, strideY);
                *(m_obstacle_back + index) = m_obstacle_list[idx_back];
            }
        }
    }

    // TOP and BOTTOM of OBSTACLE
    // fill m_obstacle_top list with top indices of obstacle and oBottom list with bottom indices
    // of obstacle
    if (m_size_obstacle_bottom > 0) {
        for (size_t k = 0; k < strideZ; ++k) {
            for (size_t i = 0; i < strideX; ++i) {
                size_t index = i + strideX * k;
                size_t idx_bottom = IX(i, 0, k, strideX, strideY);
                *(m_obstacle_bottom + index) = m_obstacle_list[idx_bottom];
            }
        }
    }
    if (m_size_obstacle_top > 0) {
        for (size_t k = 0; k < strideZ; ++k) {
            for (size_t i = 0; i < strideX; ++i) {
                size_t index = i + strideX * k;
                size_t idx_top = IX(i, strideY - 1, k, strideX, strideY);
                *(m_obstacle_top + index) = m_obstacle_list[idx_top];
            }
        }
    }

    // LEFT and RIGHT of OBSTACLE
    // fill oLeft list with left indices of obstacle and oRight list with right indices of obstacle
    if (m_size_obstacle_left > 0) {
        for (size_t k = 0; k < strideZ; ++k) {
            for (size_t j = 0; j < strideY; ++j) {
                size_t index = j + strideY * k;
                size_t idx_left = IX(0, j, k, strideX, strideY);
                *(m_obstacle_left + index) = m_obstacle_list[idx_left];
            }
        }
    }
    if (m_size_obstacle_right > 0) {
        for (size_t k = 0; k < strideZ; ++k) {
            for (size_t j = 0; j < strideY; ++j) {
                size_t index = j + strideY * k;
                size_t idx_right = IX(strideX - 1, j, k, strideX, strideY);
                *(m_obstacle_right + index) = m_obstacle_list[idx_right];
            }
        }
    }
}

//======================================== Print ===================================================
// *************************************************************************************************
/// \brief  Print obstacle infos
// *************************************************************************************************
void Obstacle::print() {
#ifndef BENCHMARKING
    size_t strideX = get_stride_x();
    size_t strideY = get_stride_y();
    size_t strideZ = get_stride_z();

    m_logger->info("-- Obstacle {}", m_name);
    m_logger->info("\t strides (x y z): {} {} {}", strideX, strideY, strideZ);
    m_logger->info("\t size of slices  (Front|Back Bottom|Top Left|Right): {}|{} {}|{} {}|{}",
                   m_size_obstacle_front, m_size_obstacle_back,
                   m_size_obstacle_bottom, m_size_obstacle_top,
                   m_size_obstacle_left, m_size_obstacle_right);
    m_logger->info("\t size of Obstacle: {}", m_size_obstacle_list);
    m_logger->info("\t coords (x y z): ({}|{}) ({}|{}) ({}|{})", m_i1, m_i2, m_j1, m_j2,
                   m_k1, m_k2);
#endif
}

//======================================== Print ===================================================
// *************************************************************************************************
/// \brief  Print detailed obstacle infos
// *************************************************************************************************
void Obstacle::print_details() {
#ifndef BENCHMARKING
    size_t Nx = m_domain.get_Nx(m_level);
    size_t Ny = m_domain.get_Ny(m_level);
    size_t strideX = get_stride_x();
    size_t strideY = get_stride_y();
    size_t strideZ = get_stride_z();
    size_t coords_i, coords_j, coords_k;

    m_logger->debug("############### OBSTACLE {} ###############", m_name);
    m_logger->debug("level: {}", m_level);
    m_logger->debug("strides (x y z): {} {} {}", strideX, strideY, strideZ);
    m_logger->debug("size of slices  (Front|Back Bottom|Top Left|Right): {}|{} {}|{} {}|{}",
                    m_size_obstacle_front, m_size_obstacle_back,
                    m_size_obstacle_bottom, m_size_obstacle_top,
                    m_size_obstacle_left, m_size_obstacle_right);
    m_logger->debug("size of Obstacle: {}", m_size_obstacle_list);
    m_logger->debug("coords (x y z): ({}|{}) ({}|{}) ({}|{})", m_i1, m_i2, m_j1, m_j2,
                    m_k1, m_k2);

    std::vector<size_t> coords;
    size_t size_front = get_size_obstacle_front();
    if (size_front > 0) {
        m_logger->debug("Front: {} | {}", m_obstacle_front[0],
                        m_obstacle_front[size_front - 1]);

        coords_k = getCoordinateK(m_obstacle_front[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_front[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_front[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Front start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacle_front[size_front - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_front[size_front - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_front[size_front - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Front end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Front size = 0");
    }

    size_t size_back = get_size_obstacle_back();
    if (size_back > 0) {
        m_logger->debug("Back: {} | {}", m_obstacle_back[0], m_obstacle_back[size_back - 1]);

        coords_k = getCoordinateK(m_obstacle_back[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_back[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_back[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Back start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacle_back[size_back - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_back[size_back - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_back[size_back - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Back end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Back size = 0");
    }

    size_t size_top = get_size_obstacle_top();
    if (size_top > 0) {
        m_logger->debug("Top: {} | {}", m_obstacle_top[0], m_obstacle_top[size_top - 1]);

        coords_k = getCoordinateK(m_obstacle_top[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_top[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_top[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Top start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacle_top[size_top - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_top[size_top - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_top[size_top - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Top end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Top size = 0");
    }

    size_t size_bottom = get_size_obstacle_bottom();
    if (size_bottom > 0) {
        m_logger->debug("Bottom: {} | {}", m_obstacle_bottom[0], m_obstacle_bottom[size_bottom - 1]);

        coords_k = getCoordinateK(m_obstacle_bottom[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_bottom[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_bottom[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Bottom start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacle_bottom[size_bottom - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_bottom[size_bottom - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_bottom[size_bottom - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Bottom end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Bottom size = 0");
    }

    size_t size_left = get_size_obstacle_left();
    if (size_left > 0) {
        m_logger->debug("Left: {} | {}", m_obstacle_left[0], m_obstacle_left[size_left - 1]);

        coords_k = getCoordinateK(m_obstacle_left[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_left[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_left[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Left start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacle_left[size_left - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_left[size_left - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_left[size_left - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Left end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Left size = 0");
    }

    size_t size_right = get_size_obstacle_right();
    if (size_right > 0) {
        m_logger->debug("Right: {} | {}", m_obstacle_right[0], m_obstacle_right[size_right - 1]);

        coords_k = getCoordinateK(m_obstacle_right[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_right[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_right[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Right start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacle_right[size_right - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacle_right[size_right - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacle_right[size_right - 1], Nx, Ny, coords_j, coords_k);
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
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(m_level);
    size_t Ny = domain->get_Ny(m_level);

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

    size_t start_index = IX(m_i1, m_j1, m_k1, Nx, Ny);
    size_t end_index = IX(m_i2, m_j2, m_k2, Nx, Ny);

    if (m_size_obstacle_front > 0) {
        size_t front_end = IX(m_i2, m_j2, m_k1, Nx, Ny);
        if (start_index != m_obstacle_front[0] ||
            front_end != m_obstacle_front[m_size_obstacle_front - 1]) {
            message += "first or last index of obstacle front list not correct ("
                       + std::to_string(start_index) + "|" + std::to_string(m_obstacle_front[0])
                       + ")(" + std::to_string(front_end) + "|"
                       + std::to_string(m_obstacle_front[m_size_obstacle_front - 1]) + ")\n";
        }
    }
    if (m_size_obstacle_back > 0) {
        size_t back_start = IX(m_i1, m_j1, m_k2, Nx, Ny);
        if (back_start != m_obstacle_back[0] ||
            end_index != m_obstacle_back[m_size_obstacle_back - 1]) {
            message += "first or last index of obstacle back list not correct ("
                       + std::to_string(back_start) + "|" + std::to_string(m_obstacle_back[0])
                       + ")(" + std::to_string(end_index) + "|"
                       + std::to_string(m_obstacle_back[m_size_obstacle_back - 1]) + ")\n";
        }
    }
    if (m_size_obstacle_bottom > 0) {
        size_t bottom_end = IX(m_i2, m_j1, m_k2, Nx, Ny);
        if (start_index != m_obstacle_bottom[0] ||
            bottom_end != m_obstacle_bottom[m_size_obstacle_bottom - 1]) {
            message += "first or last index of obstacle bottom list not correct ("
                       + std::to_string(start_index) + "|" + std::to_string(m_obstacle_bottom[0])
                       + ")(" + std::to_string(bottom_end) + "|"
                       + std::to_string(m_obstacle_bottom[m_size_obstacle_bottom - 1]) + ")\n";
        }
    }
    if (m_size_obstacle_top > 0) {
        size_t top_start = IX(m_i1, m_j2, m_k1, Nx, Ny);
        if (top_start != m_obstacle_top[0] ||
            end_index != m_obstacle_top[m_size_obstacle_top - 1]) {
            message += "first or last index of obstacle top list not correct ("
                       + std::to_string(top_start) + "|" + std::to_string(m_obstacle_top[0])
                       + ")(" + std::to_string(end_index) + "|"
                       + std::to_string(m_obstacle_top[m_size_obstacle_top - 1]) + ")\n";
        }
    }
    if (m_size_obstacle_left > 0) {
        size_t left_end = IX(m_i1, m_j2, m_k2, Nx, Ny);
        if (start_index != m_obstacle_left[0] ||
            left_end != m_obstacle_left[m_size_obstacle_left - 1]) {
            message += "first or last index of obstacle left list not correct ("
                       + std::to_string(start_index) + "|" + std::to_string(m_obstacle_left[0])
                       + ")(" + std::to_string(left_end) + "|"
                       + std::to_string(m_obstacle_left[m_size_obstacle_left - 1]) + ")\n";
        }
    }
    if (m_size_obstacle_right > 0) {
        size_t right_start = IX(m_i2, m_j1, m_k1, Nx, Ny);
        if (right_start != m_obstacle_right[0] ||
            end_index != m_obstacle_right[m_size_obstacle_right - 1]) {
            message += "first or last index of obstacle right list not correct ("
                       + std::to_string(right_start) + "|" + std::to_string(m_obstacle_right[0])
                       + ")(" + std::to_string(end_index) + "|"
                       + std::to_string(m_obstacle_right[m_size_obstacle_right - 1]) + ")\n";
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

//======================================== Is obstacle cell ========================================
// *************************************************************************************************
/// \brief  Check if cell is an obstacle cell
/// \param  i x-coordinate
/// \param  j y-coordinate
/// \param  k z-coordinate
/// \return  bool true if yes false if no
// ***************************************************************************************
bool Obstacle::is_obstacle_cell(size_t i, size_t j, size_t k) const {
    return m_i1 <= i && i <= m_i2 &&
        m_j1 <= j && j <= m_j2 &&
        m_k1 <= k && k <= m_k2;
}

real constexpr det3(real a1, real a2, real a3,
        real b1, real b2, real b3,
        real c1, real c2, real c3) {
    return a1 * (b2 * c3 - b3 * c2)
        - b1 * (a2 * c3 - a3 * c2)
        + c1 * (a2 * b3 - a3 * b2);
}

real constexpr eps = 10E-6;

constexpr int cuboid_points[8][3] = {
    {0, 2, 4},
    {1, 2, 4},
    {1, 3, 4},
    {0, 3, 4},
    {0, 2, 5},
    {1, 2, 5},
    {1, 3, 5},
    {0, 3, 5}
};

constexpr int surface_points[6][3] = {
    {0, 1, 4},
    {1, 2, 5},
    {3, 2, 6},
    {1, 3, 7},
    {0, 1, 3},
    {4, 5, 7},
};

// ======================== Tests if line crosses obstacle ====================
// ****************************************************************************
/// \brief  Check if any cell  between two pins is inside the obstacle
/// \param  i x-coordinate (start_point)
/// \param  j y-coordinate (start_point)
/// \param  k z-coordinate (start_point)
/// \param  i x-coordinate (end_point)
/// \param  j y-coordinate (end_point)
/// \param  k z-coordinate (end_point)
/// \return  bool true if yes false if no
// ****************************************************************************
bool Obstacle::line_crosses(int i0, int j0, int k0, int i, int j, int k) const {
    // m_logger->error("i:{} j:{} k:{}", i0, j0, k0);
    // m_logger->warn("i:{} j:{} k:{}", i, j, k);
    bool blocked;
    const auto di = i - i0;
    const auto dj = j - j0;
    const auto dk = k - k0;
    // m_logger->warn("i:{} j:{} k:{}", di, dj, dk);
    const auto point_i1 = get_coordinates_i1();
    const auto point_i2 = get_coordinates_i2();
    const auto point_j1 = get_coordinates_j1();
    const auto point_j2 = get_coordinates_j2();
    const auto point_k1 = get_coordinates_k1();
    const auto point_k2 = get_coordinates_k2();
    real indeces[6] = {static_cast<real>(point_i1), static_cast<real>(point_i2),
                       static_cast<real>(point_j1), static_cast<real>(point_j2),
                       static_cast<real>(point_k1), static_cast<real>(point_k2)};

    if (is_obstacle_cell(i, j, k))
        return true;

    for (int surface_id=0; surface_id < 6; ++surface_id) {
        // m_logger->debug("surface_id: {}", surface_id);
        auto s = surface_points[surface_id];
        auto p1 = cuboid_points[s[0]];
        auto p2 = cuboid_points[s[1]];
        auto p3 = cuboid_points[s[2]];
        // m_logger->debug("surface_point: ({},{},{})",
        //        indeces[p1[0]], indeces[p1[1]], indeces[p1[2]]);
        // m_logger->debug("surface_point: ({},{},{})",
        //        indeces[p2[0]], indeces[p2[1]], indeces[p2[2]]);
        // m_logger->debug("surface_point: ({},{},{})",
        //        indeces[p3[0]], indeces[p3[1]], indeces[p3[2]]);

        // surface vector 1
        auto svi1 = indeces[p2[0]] - indeces[p1[0]];  // saving dx
        auto svj1 = indeces[p2[1]] - indeces[p1[1]];
        auto svk1 = indeces[p2[2]] - indeces[p1[2]];
        // m_logger->debug("surface_vector1: ({},{},{})", svi1, svj1, svk1);

        // surface vector 2
        auto svi2 = indeces[p3[0]] - indeces[p1[0]];
        auto svj2 = indeces[p3[1]] - indeces[p1[1]];
        auto svk2 = indeces[p3[2]] - indeces[p1[2]];
        // m_logger->debug("surface_vector1: ({},{},{})", svi2, svj2, svk2);

        // p1 + l*sv1 + m*sv2 = n*d + c0 <=>
        // l*sv1 + m*sv2 - n*d = c0 - p1
        auto det_A = det3(di, -svi1, -svi2,
            dj, -svj1, -svj2,
            dk, -svk1, -svk2);
        // auto det_A = det3(svi1, svj1, svk1,
        //     svi2, svj2, svk2,
        //     -di, -dj, -dk);

        // okay det is small, its never gonna meet
        // m_logger->info("| {}, {}, {} |", di, -svi1, -svi2);
        // m_logger->info("| {}, {}, {} |", dj, -svj1, -svj2);
        // m_logger->info("| {}, {}, {} |", dk, -svk1, -svk2);
        // m_logger->error("|det|: {}", det_A);
        if (fabs(det_A) < eps) {
            continue;
        }

        // rhs of les (c - p1)
        auto ddi = indeces[p1[0]] - static_cast<real>(i0);
        auto ddj = indeces[p1[1]] - static_cast<real>(j0);
        auto ddk = indeces[p1[2]] - static_cast<real>(k0);
        // m_logger->debug("rhs: ({},{},{})", ddi, ddj, ddk);

        auto det_Ax = det3(ddi, -svi1, -svi2,
            ddj, -svj1, -svj2,
            ddk, -svk1, -svk2);
        auto det_Ay = det3(di, ddi, -svi2,
            dj, ddj, -svj2,
            dk, ddk, -svk2);
        auto det_Az = det3(di, -svi1, ddi,
            dj, -svj1, ddj,
            dk, -svk1, ddk);
        // auto det_Ax = det3(ddi, svi1, -di,
        //     ddj, svj2, -dj,
        //     ddk, svk2, -dk);
        // auto det_Ay = det3(svi1, ddi, -di,
        //     svj1, ddj, -dj,
        //     svk1, ddk, -dk);
        // auto det_Az = det3(svi1, svi2, ddi,
        //     svj1, svj2, ddj,
        //     svk1, svk2, ddk);

        auto sx = det_Ax / det_A;  // sx : l
        auto sy = det_Ay / det_A;  // sy : m
        auto sz = det_Az / det_A;  // sz : n

        // m_logger->debug("{:6f} < eps < {6f}", -eps, 1.0 + eps);
        // m_logger->debug("x,y,z: ({:6f},{:6f},{:6f})", sx, sy, sz);
        // m_logger->debug("xp: ({},{},{})",
        //        indeces[p1[0]] + svi1*sy + svi2*sz,
        //        indeces[p1[1]] + svj1*sy + svj2*sz,
        //        indeces[p1[2]] + svk1*sy + svk2*sz);
        // m_logger->debug("xp: ({},{},{})",
        //        static_cast<real>(i0) + di*sx,
        //        static_cast<real>(j0) + dj*sx,
        //        static_cast<real>(k0) + dk*sx);

        // m_logger->debug("{}", sx >= -eps && sx <= 1.0 + eps);
        // m_logger->debug("{}", sy >= -eps && sy <= 1.0 + eps);
        // m_logger->debug("{}", sz >= -eps && sz <= 1.0 + eps);
        blocked = sx >= -eps && sx <= 1.0 + eps
                && sy >= -eps && sy <= 1.0 + eps
                && sz >= -eps && sz <= 1.0 + eps;

        // std::cout << di << "," << i << std::endl;
        // std::cout << dj << "," << j << std::endl;
        // std::cout << dk << "," << k << std::endl;

        // std::cout << surface_id << std::endl;

        // std::cout << sx << " : " << sy << " : " << sz << std::endl;
        // std::cout << sy << std::endl;
        // std::cout << sz << std::endl;

        if (blocked) {
            return true;
        }
    }

    return false;
}


//======================================== Match grid ====================================
// ***************************************************************************************
/// \brief  Snaps value to grid discretisation
/// \param  obstacle_coordinate Coordinate of obstacle
/// \param  spacing dx/dy/dz
/// \param  start_coordinate X1/Y1/Z1
/// \return real Calculated real grid coordinate
// *************************************************************************************************
real Obstacle::match_grid(real obstacle_coordinate, real spacing, real start_coordinate) {
    return get_matching_index(obstacle_coordinate, spacing, start_coordinate)
           * spacing + start_coordinate;
}

int Obstacle::get_matching_index(real obstacle_coordinate, real spacing, real start_coordinate) {
    return static_cast<int>(round((-start_coordinate + obstacle_coordinate) / spacing));
}

//======================================== Remove cells at boundary ================================
// *************************************************************************************************
/// \brief  Remove obstacle patch facing the boundary
/// \param  level Multigrid level
// *************************************************************************************************
void Obstacle::remove_cells_at_boundary(size_t level) {
    Domain *domain = Domain::getInstance();
    if (m_k1 <= domain->get_index_z1(level)) {
        m_size_obstacle_front = 0;
    }
    if (m_k2 >= domain->get_index_z2(level)) {
        m_size_obstacle_back = 0;
    }
    if (m_j1 <= domain->get_index_y1(level)) {
        m_size_obstacle_bottom = 0;
    }
    if (m_j2 >= domain->get_index_y2(level)) {
        m_size_obstacle_top = 0;
    }
    if (m_i1 <= domain->get_index_x1(level)) {
        m_size_obstacle_left = 0;
    }
    if (m_i2 >= domain->get_index_x2(level)) {
        m_size_obstacle_right = 0;
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
    switch (p) {
        case FRONT:
            delete[] m_obstacle_front;
            m_obstacle_front = indices;
            m_size_obstacle_front = size;
            break;
        case BACK:
            delete[] m_obstacle_back;
            m_obstacle_back = indices;
            m_size_obstacle_back = size;
            break;
        case BOTTOM:
            delete[] m_obstacle_bottom;
            m_obstacle_bottom = indices;
            m_size_obstacle_bottom = size;
            break;
        case TOP:
            delete[] m_obstacle_top;
            m_obstacle_top = indices;
            m_size_obstacle_top = size;
            break;
        case LEFT:
            delete[] m_obstacle_left;
            m_obstacle_left = indices;
            m_size_obstacle_left = size;
            break;
        case RIGHT:
            delete[] m_obstacle_right;
            m_obstacle_right = indices;
            m_size_obstacle_right = size;
            break;
        default:
#ifndef BENCHMARKING
            m_logger->warn("wrong patch: {}", p);
#endif
            break;
    }
}

//======================================== has overlap =============================================
// *************************************************************************************************
/// calculate indices of area to be excluded. o1_coordinate == o2_coordinate only if the length
/// of the patch of both obstacles are the same
/// \param o1 Obstacle 1
/// \param o2 Obstacle 2
/// \param o1_coordinate calculated coordinate of obstacle 1
/// \param o2_coordinate calculated coordinate of obstacle 2
/// \param direction X/Y/Z axis
/// \param start true = (i/j/k)1 or false = (i/j/k)2
// *************************************************************************************************
void Obstacle::calculate_area_index(
        Obstacle *o1, Obstacle *o2,
        size_t *o1_coordinate, size_t *o2_coordinate,
        CoordinateAxis direction,
        bool start) {
    if (direction == CoordinateAxis::X) {
        if (start) {
            *o1_coordinate = o1->get_coordinates_i1();
            *o2_coordinate = o2->get_coordinates_i1();
            if (o1->get_coordinates_i1() > o2->get_coordinates_i1()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o2_coordinate = o1->get_coordinates_i1() + 1;
            } else if (o1->get_coordinates_i1() < o2->get_coordinates_i1()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o1_coordinate = o2->get_coordinates_i1() + 1;
            }
        } else {
            *o1_coordinate = o1->get_coordinates_i2();
            *o2_coordinate = o2->get_coordinates_i2();
            if (o1->get_coordinates_i2() < o2->get_coordinates_i2()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o2_coordinate = o1->get_coordinates_i2() - 1;
            } else if (o1->get_coordinates_i2() > o2->get_coordinates_i2()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o1_coordinate = o2->get_coordinates_i2() - 1;
            }
        }
    }
    if (direction == CoordinateAxis::Y) {
        if (start) {
            *o1_coordinate = o1->get_coordinates_j1();
            *o2_coordinate = o2->get_coordinates_j1();
            if (o1->get_coordinates_j1() > o2->get_coordinates_j1()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o2_coordinate = o1->get_coordinates_j1() + 1;
            } else if (o1->get_coordinates_j1() < o2->get_coordinates_j1()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o1_coordinate = o2->get_coordinates_j1() + 1;
            }
        } else {
            *o1_coordinate = o1->get_coordinates_j2();
            *o2_coordinate = o2->get_coordinates_j2();
            if (o1->get_coordinates_j2() < o2->get_coordinates_j2()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o2_coordinate = o1->get_coordinates_j2() - 1;
            } else if (o1->get_coordinates_j2() > o2->get_coordinates_j2()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o1_coordinate = o2->get_coordinates_j2() - 1;
            }
        }
    }
    if (direction == CoordinateAxis::Z) {
        if (start) {
            *o1_coordinate = o1->get_coordinates_k1();
            *o2_coordinate = o2->get_coordinates_k1();
            if (o1->get_coordinates_k1() > o2->get_coordinates_k1()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o2_coordinate = o1->get_coordinates_k1() + 1;
            } else if (o1->get_coordinates_k1() < o2->get_coordinates_k1()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o1_coordinate = o2->get_coordinates_k1() + 1;
            }
        } else {
            *o1_coordinate = o1->get_coordinates_k2();
            *o2_coordinate = o2->get_coordinates_k2();
            if (o1->get_coordinates_k2() < o2->get_coordinates_k2()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o2_coordinate = o1->get_coordinates_k2() - 1;
            } else if (o1->get_coordinates_k2() > o2->get_coordinates_k2()) {
                // do not remove inner edge, can be accessed by SL Advection Solver
                *o1_coordinate = o2->get_coordinates_k2() - 1;
            }
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
bool Obstacle::remove_circular_constraints(Obstacle *o1, Obstacle *o2) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger("Obstacle");
#endif

    bool overlap = circular_constraints_x_direction(o1, o2);
    overlap = overlap || circular_constraints_y_direction(o1, o2);
    overlap = overlap || circular_constraints_z_direction(o1, o2);
#ifndef BENCHMARKING
    if (overlap) {
        logger->debug("{} is next to {}", o1->get_name(), o2->get_name());
    }
#endif
    return overlap;
}

//============================ circular constraints in x direction =================================
// *************************************************************************************************
/// \brief removes circular constraints in x direction if two obstacles are next to each other
/// \param o1 obstacle 1
/// \param o2 obstacle 2
// *************************************************************************************************
bool Obstacle::circular_constraints_x_direction(Obstacle *o1, Obstacle *o2) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger("Obstacle");
#endif
    bool overlap = false;

    auto domain = Domain::getInstance();
    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();

    if (o1->get_coordinates_i2() + 1 == o2->get_coordinates_i1()) {
        std::swap(o1, o2);
    }

    if (o1->get_coordinates_i1() - 1 == o2->get_coordinates_i2()) {
        bool j_overlap = has_overlap(o1->get_coordinates_j1(), o1->get_coordinates_j2(),
                                     o2->get_coordinates_j1(), o2->get_coordinates_j2());
        bool k_overlap = has_overlap(o1->get_coordinates_k1(), o1->get_coordinates_k2(),
                                     o2->get_coordinates_k1(), o2->get_coordinates_k2());
        if (j_overlap && k_overlap) {
#ifndef BENCHMARKING
            logger->debug("obstacles are next to each other. Working on {} left side and on {} right side.",
                          o1->get_name(), o2->get_name());
#endif
            // another obstacle (o2) at the left side of o1
            overlap = true;
            // calculate coordinates of area which should be removed
            // the area is for both obstacle the same only if there are equally long
            size_t o1_x1 = o1->get_coordinates_i1();
            size_t o2_x2 = o2->get_coordinates_i2();

            size_t o1_y1;
            size_t o2_y1;
            Obstacle::calculate_area_index(o1, o2, &o1_y1, &o2_y1, CoordinateAxis::Y, true);

            size_t o1_y2;
            size_t o2_y2;
            Obstacle::calculate_area_index(o1, o2, &o1_y2, &o2_y2, CoordinateAxis::Y, false);

            size_t o1_z1;
            size_t o2_z1;
            Obstacle::calculate_area_index(o1, o2, &o1_z1, &o2_z1, CoordinateAxis::Z, true);

            size_t o1_z2;
            size_t o2_z2;
            Obstacle::calculate_area_index(o1, o2, &o1_z2, &o2_z2, CoordinateAxis::Z, false);

#ifndef BENCHMARKING
            logger->debug("removing indices in area ({}) ({}|{}) ({}|{}) for {}",
                          o1_x1, o1_y1, o1_y2, o1_z1, o1_z2, o1->get_name());
            logger->debug("removing indices in area ({}) ({}|{}) ({}|{}) for {}",
                          o2_x2, o2_y1, o2_y2, o2_z1, o2_z2, o2->get_name());
#endif

            std::vector<size_t> o1_new;
            o1_new.reserve(o1->get_size_obstacle_left());
            std::vector<size_t> o2_new;
            o2_new.reserve(o1->get_size_obstacle_right());

            size_t o1_counter_old = 0;
            size_t o1_smallest_removing_index = IX(o1_x1, o1_y1, o1_z1, Nx, Ny);
            size_t o1_current_index = o1->get_obstacle_left()[o1_counter_old];
            while (o1_current_index < o1_smallest_removing_index) {
                o1_new.push_back(o1_current_index);
                o1_counter_old++;
                o1_current_index = o1->get_obstacle_left()[o1_counter_old];
            }
            size_t o1_new_size_left = o1_counter_old;

            size_t o2_counter_old = 0;
            size_t o2_smallest_removing_index = IX(o2_x2, o2_y1, o2_z1, Nx, Ny);
            size_t o2_current_index = o2->get_obstacle_right()[o2_counter_old];
            while (o2_current_index < o2_smallest_removing_index) {
                o2_new.push_back(o2_current_index);
                o2_counter_old++;
                o2_current_index = o2->get_obstacle_right()[o2_counter_old];
            }
            size_t o2_new_size_right = o2_counter_old;

            size_t o1_current_y = o1_y1;
            size_t o1_current_z = o1_z1;
            size_t o1_removing_index = IX(o1_x1, o1_current_y, o1_current_z, Nx, Ny);
            bool o1_end = false;

            size_t o2_current_y = o2_y1;
            size_t o2_current_z = o2_z1;
            size_t o2_removing_index = IX(o2_x2, o2_current_y, o2_current_z, Nx, Ny);
            bool o2_end = false;
            for (; o1_counter_old < o1->get_size_obstacle_left()
                   && o2_counter_old < o2->get_size_obstacle_right() && !o1_end && !o2_end;
                   o1_counter_old++, o2_counter_old++) {
                o1_current_index = o1->get_obstacle_left()[o1_counter_old];
                o2_current_index = o2->get_obstacle_right()[o2_counter_old];
                if (o1_current_index != o1_removing_index) {
                    o1_new.push_back(o1_current_index);
                    o1_new_size_left++;
                } else {
                    o1_current_y++;
                    if (o1_current_y > o1_y2) {
                        o1_current_y = o1_y1;
                        o1_current_z++;
                        if (o1_current_z > o1_z2) {
                            o1_end = true;
                        }
                    }
                    o1_removing_index = IX(o1_x1, o1_current_y, o1_current_z, Nx, Ny);
                }
                if (o2_current_index != o2_removing_index) {
                    o2_new.push_back(o2_current_index);
                    o2_new_size_right++;
                } else {
                    o2_current_y++;
                    if (o2_current_y > o2_y2) {
                        o2_current_y = o2_y1;
                        o2_current_z++;
                        if (o2_current_z > o2_z2) {
                            o2_end = true;
                        }
                    }
                    o2_removing_index = IX(o2_x2, o2_current_y, o2_current_z, Nx, Ny);
                }
            }

            if (!o1_end) {
                for (; o1_counter_old < o1->get_size_obstacle_left() && o1_current_z <= o1_z2; o1_counter_old++) {
                    o1_current_index = o1->get_obstacle_left()[o1_counter_old];
                    if (o1_current_index != o1_removing_index) {
                        o1_new.push_back(o1_current_index);
                        o1_new_size_left++;
                    } else {
                        o1_current_y++;
                        if (o1_current_y > o1_y2) {
                            o1_current_y = o1_y1;
                            o1_current_z++;
                        }
                        o1_removing_index = IX(o1_x1, o1_current_y, o1_current_z, Nx, Ny);
                    }
                }
            }

            if (!o2_end) {
                for (; o2_counter_old < o2->get_size_obstacle_right() && o2_current_z <= o2_z2; o2_counter_old++) {
                    o2_current_index = o2->get_obstacle_right()[o2_counter_old];
                    if (o2_current_index != o2_removing_index) {
                        o2_new.push_back(o2_current_index);
                        o2_new_size_right++;
                    } else {
                        o2_current_y++;
                        if (o2_current_y > o2_y2) {
                            o2_current_y = o2_y1;
                            o2_current_z++;
                        }
                        o2_removing_index = IX(o2_x2, o2_current_y, o2_current_z, Nx, Ny);
                    }
                }
            }

            for (; o1_counter_old < o1->get_size_obstacle_left(); o1_counter_old++) {
                o1_new.push_back(o1->get_obstacle_left()[o1_counter_old]);
                o1_new_size_left++;
            }
            o1_new.resize(o1_new_size_left);

            size_t o1_diff_target = (o1_z2 - o1_z1 + 1) * (o1_y2 - o1_y1 + 1);
            size_t o1_diff_actual = o1->get_size_obstacle_left() - o1_new_size_left;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} left patch: {} -> {} ({}|{})",
                          o1->get_name(), o1->get_size_obstacle_left(), o1_new_size_left,
                          o1_diff_target, o1_diff_actual);
#endif
            size_t *o1_new_data = new size_t[o1_new_size_left];
            std::copy(o1_new.begin(), o1_new.end(), o1_new_data);
            o1->replace_patch(o1_new_data, o1_new_size_left, Patch::LEFT);

            for (; o2_counter_old < o2->get_size_obstacle_right(); o2_counter_old++) {
                o2_new.push_back(o2->get_obstacle_right()[o2_counter_old]);
                o2_new_size_right++;
            }
            o2_new.resize(o2_new_size_right);

            size_t o2_diff_target = (o2_z2 - o2_z1 + 1) * (o2_y2 - o2_y1 + 1);
            size_t o2_diff_actual = o2->get_size_obstacle_right() - o2_new_size_right;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} right patch: {} -> {} ({}|{})",
                          o2->get_name(), o2->get_size_obstacle_right(), o2_new_size_right,
                          o2_diff_target, o2_diff_actual);
#endif
            size_t *o2_new_data = new size_t[o2_new_size_right];
            std::copy(o2_new.begin(), o2_new.end(), o2_new_data);
            o2->replace_patch(o2_new_data, o2_new_size_right, Patch::RIGHT);
        }
    }
    return overlap;
}

//============================ circular constraints in y direction =================================
// *************************************************************************************************
/// \brief removes circular constraints in y direction if two obstacles are next to each other
/// \param o1 obstacle 1
/// \param o2 obstacle 2
// *************************************************************************************************
bool Obstacle::circular_constraints_y_direction(Obstacle *o1, Obstacle *o2) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger("Obstacle");
#endif
    bool overlap = false;

    auto domain = Domain::getInstance();
    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();

    if (o1->get_coordinates_j2() + 1 == o2->get_coordinates_j1()) {
        std::swap(o1, o2);
    }

    if (o1->get_coordinates_j1() - 1 == o2->get_coordinates_j2()) {
        bool i_overlap = has_overlap(
                o1->get_coordinates_i1(), o1->get_coordinates_i2(),
                o2->get_coordinates_i1(), o2->get_coordinates_i2());
        bool k_overlap = has_overlap(
                o1->get_coordinates_k1(), o1->get_coordinates_k2(),
                o2->get_coordinates_k1(), o2->get_coordinates_k2());
        if (i_overlap && k_overlap) {
#ifndef BENCHMARKING
            logger->debug("obstacles are next to each other. Working on {} bottom side and on {} top side",
                          o1->get_name(), o2->get_name());
#endif
            overlap = true;
            // calculate coordinates of area which should be removed
            // the area is for both obstacle the same only if there are equally long

            size_t o1_x1;
            size_t o2_x1;
            Obstacle::calculate_area_index(o1, o2, &o1_x1, &o2_x1, CoordinateAxis::X, true);

            size_t o1_x2;
            size_t o2_x2;
            Obstacle::calculate_area_index(o1, o2, &o1_x2, &o2_x2, CoordinateAxis::X, false);

            size_t o1_y1 = o1->get_coordinates_j1();
            size_t o2_y2 = o2->get_coordinates_j2();

            size_t o1_z1;
            size_t o2_z1;
            Obstacle::calculate_area_index(o1, o2, &o1_z1, &o2_z1, CoordinateAxis::Z, true);

            size_t o1_z2;
            size_t o2_z2;
            Obstacle::calculate_area_index(o1, o2, &o1_z2, &o2_z2, CoordinateAxis::Z, false);

#ifndef BENCHMARKING
            logger->debug("removing indices in area ({}|{}) ({}) ({}|{}) for {}",
                          o1_x1, o1_x2, o1_y1, o1_z1, o1_z2, o1->get_name());
            logger->debug("removing indices in area ({}|{}) ({}) ({}|{}) for {}",
                          o2_x1, o2_x2, o2_y2, o2_z1, o2_z2, o2->get_name());
#endif

            std::vector<size_t> o1_new;
            o1_new.reserve(o1->get_size_obstacle_bottom());
            std::vector<size_t> o2_new;
            o2_new.reserve(o2->get_size_obstacle_top());

            size_t o1_counter_old = 0;
            size_t o1_smallest_removing_index = IX(o1_x1, o1_y1, o1_z1, Nx, Ny);
            size_t o1_current_index = o1->get_obstacle_bottom()[o1_counter_old];
            while (o1_current_index < o1_smallest_removing_index) {
                o1_new.push_back(o1_current_index);
                o1_counter_old++;
                o1_current_index = o1->get_obstacle_bottom()[o1_counter_old];
            }
            size_t o1_new_size_bottom = o1_counter_old;

            size_t o2_counter_old = 0;
            size_t o2_smallest_removing_index = IX(o2_x1, o2_y2, o2_z1, Nx, Ny);
            size_t o2_current_index = o2->get_obstacle_top()[o2_counter_old];
            while (o2_current_index < o2_smallest_removing_index) {
                o2_new.push_back(o2_current_index);
                o2_counter_old++;
                o2_current_index = o2->get_obstacle_top()[o2_counter_old];
            }
            size_t o2_new_size_top = o2_counter_old;

            size_t o1_current_x = o1_x1;
            size_t o1_current_z = o1_z1;
            size_t o1_removing_index = IX(o1_current_x, o1_y1, o1_current_z, Nx, Ny);  // equals smallest removing index
            bool o1_end = false;

            size_t o2_current_x = o2_x1;
            size_t o2_current_z = o2_z1;
            size_t o2_removing_index = IX(o2_current_x, o2_y2, o2_current_z, Nx, Ny);
            bool o2_end = false;
            for (; o1_counter_old < o1->get_size_obstacle_bottom()
                   && o2_counter_old < o2->get_size_obstacle_top() && !o1_end && !o2_end;
                   o1_counter_old++, o2_counter_old++) {
                o1_current_index = o1->get_obstacle_bottom()[o1_counter_old];
                o2_current_index = o2->get_obstacle_top()[o2_counter_old];
                if (o1_current_index != o1_removing_index) {
                    o1_new.push_back(o1_current_index);
                    o1_new_size_bottom++;
                } else {
                    o1_current_x++;
                    if (o1_current_x > o1_x2) {
                        o1_current_x = o1_x1;
                        o1_current_z++;
                        if (o1_current_z > o1_z2) {
                            o1_end = true;
                        }
                    }
                    o1_removing_index = IX(o1_current_x, o1_y1, o1_current_z, Nx, Ny);
                }
                if (o2_current_index != o2_removing_index) {
                    o2_new.push_back(o2_current_index);
                    o2_new_size_top++;
                } else {
                    o2_current_x++;
                    if (o2_current_x > o2_x2) {
                        o2_current_x = o2_x1;
                        o2_current_z++;
                        if (o2_current_z > o2_z2) {
                            o2_end = true;
                        }
                    }
                    o2_removing_index = IX(o2_current_x, o2_y2, o2_current_z, Nx, Ny);
                }
            }

            if (!o1_end) {
                for (; o1_counter_old < o1->get_size_obstacle_bottom() && o1_current_z <= o1_z2; o1_counter_old++) {
                    o1_current_index = o1->get_obstacle_bottom()[o1_counter_old];
                    if (o1_current_index != o1_removing_index) {
                        o1_new.push_back(o1_current_index);
                        o1_new_size_bottom++;
                    } else {
                        o1_current_x++;
                        if (o1_current_x > o1_x2) {
                            o1_current_x = o1_x1;
                            o1_current_z++;
                        }
                        o1_removing_index = IX(o1_current_x, o1_y1, o1_current_z, Nx, Ny);
                    }
                }
            }

            if (!o2_end) {
                for (; o2_counter_old < o2->get_size_obstacle_top() && o2_current_z <= o2_z2; o2_counter_old++) {
                    o2_current_index = o2->get_obstacle_top()[o2_counter_old];
                    if (o2_current_index != o2_removing_index) {
                        o2_new.push_back(o2_current_index);
                        o2_new_size_top++;
                    } else {
                        o2_current_x++;
                        if (o2_current_x > o2_x2) {
                            o2_current_x = o2_x1;
                            o2_current_z++;
                        }
                        o2_removing_index = IX(o2_current_x, o2_y2, o2_current_z, Nx, Ny);
                    }
                }
            }

            for (; o1_counter_old < o1->get_size_obstacle_bottom(); o1_counter_old++) {
                o1_new.push_back(o1->get_obstacle_bottom()[o1_counter_old]);
                o1_new_size_bottom++;
            }
            o1_new.resize(o1_new_size_bottom);

            size_t o1_diff_target = (o1_x2 - o1_x1 + 1) * (o1_z2 - o1_z1 + 1);
            size_t o1_diff_actual = o1->get_size_obstacle_bottom() - o1_new_size_bottom;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} bottom patch: {} -> {} ({}|{})",
                          o1->get_name(), o1->get_size_obstacle_bottom(), o1_new_size_bottom,
                          o1_diff_target, o1_diff_actual);
#endif
            size_t *o1_new_data = new size_t[o1_new_size_bottom];
            std::copy(o1_new.begin(), o1_new.end(), o1_new_data);
            o1->replace_patch(o1_new_data, o1_new_size_bottom, Patch::BOTTOM);

            for (; o2_counter_old < o2->get_size_obstacle_top(); o2_counter_old++) {
                o2_new.push_back(o2->get_obstacle_top()[o2_counter_old]);
                o2_new_size_top++;
            }
            o2_new.resize(o2_new_size_top);

            size_t o2_diff_target = (o2_x2 - o2_x1 + 1) * (o2_z2 - o2_z1 + 1);
            size_t o2_diff_actual = o2->get_size_obstacle_top() - o2_new_size_top;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} top patch: {} -> {} ({}|{})",
                          o2->get_name(), o2->get_size_obstacle_top(), o2_new_size_top,
                          o2_diff_target, o2_diff_actual);
#endif
            size_t *o2_new_data = new size_t[o2_new_size_top];
            std::copy(o2_new.begin(), o2_new.end(), o2_new_data);
            o2->replace_patch(o2_new_data, o2_new_size_top, Patch::TOP);
        }
    }
    return overlap;
}

//============================ circular constraints in z direction =================================
// *************************************************************************************************
/// \brief removes circular constraints in z direction if two obstacles are next to each other
/// \param o1 obstacle 1
/// \param o2 obstacle 2
// *************************************************************************************************
bool Obstacle::circular_constraints_z_direction(Obstacle *o1, Obstacle *o2) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger("Obstacle");
#endif
    bool overlap = false;

    auto domain = Domain::getInstance();
    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();

    if (o1->get_coordinates_k2() + 1 == o2->get_coordinates_k1()) {
        std::swap(o1, o2);
    }

    if (o1->get_coordinates_k1() - 1 == o2->get_coordinates_k2()) {
        bool i_overlap = has_overlap(o1->get_coordinates_i1(), o1->get_coordinates_i2(),
                                     o2->get_coordinates_i1(), o2->get_coordinates_i2());
        bool j_overlap = has_overlap(o1->get_coordinates_j1(), o1->get_coordinates_j2(),
                                     o2->get_coordinates_j1(), o2->get_coordinates_j2());
        if (i_overlap && j_overlap) {
#ifndef BENCHMARKING
            logger->debug("obstacles are next to each other."
                          "Working on {} front side and on {} back side.",
                          o1->get_name(), o2->get_name());
#endif
            // another obstacle (o2) at the front side of o1
            overlap = true;
            // calculate coordinates of area which should be removed
            // the area is for both obstacle the same only if there are equally long

            size_t o1_x1;
            size_t o2_x1;
            Obstacle::calculate_area_index(o1, o2, &o1_x1, &o2_x1, CoordinateAxis::X, true);

            size_t o1_x2;
            size_t o2_x2;
            Obstacle::calculate_area_index(o1, o2, &o1_x2, &o2_x2, CoordinateAxis::X, false);

            size_t o1_y1;
            size_t o2_y1;
            Obstacle::calculate_area_index(o1, o2, &o1_y1, &o2_y1, CoordinateAxis::Y, true);

            size_t o1_y2;
            size_t o2_y2;
            Obstacle::calculate_area_index(o1, o2, &o1_y2, &o2_y2, CoordinateAxis::Y, false);

            size_t o1_z1 = o1->get_coordinates_k1();
            size_t o2_z2 = o2->get_coordinates_k2();

#ifndef BENCHMARKING
            logger->debug("removing indices in area ({}|{}) ({}|{}) ({}) for {}",
                          o1_x1, o1_x2, o1_y1, o1_y2, o1_z1, o1->get_name());
            logger->debug("removing indices in area ({}|{}) ({}|{}) ({}) for {}",
                          o2_x1, o2_x2, o2_y1, o2_y2, o2_z2, o2->get_name());
#endif

            std::vector<size_t> o1_new;
            o1_new.reserve(o1->get_size_obstacle_front());
            std::vector<size_t> o2_new;
            o2_new.reserve(o1->get_size_obstacle_back());

            size_t o1_counter_old = 0;
            size_t o1_smallest_removing_index = IX(o1_x1, o1_y1, o1_z1, Nx, Ny);
            size_t o1_current_index = o1->get_obstacle_front()[o1_counter_old];
            while (o1_current_index < o1_smallest_removing_index) {
                o1_new.push_back(o1_current_index);
                o1_counter_old++;
                o1_current_index = o1->get_obstacle_front()[o1_counter_old];
            }
            size_t o1_new_size_front = o1_counter_old;

            size_t o2_counter_old = 0;
            size_t o2_smallest_removing_index = IX(o2_x1, o2_y1, o2_z2, Nx, Ny);
            size_t o2_current_index = o2->get_obstacle_back()[o2_counter_old];
            while (o2_current_index < o2_smallest_removing_index) {
                o2_new.push_back(o2_current_index);
                o2_counter_old++;
                o2_current_index = o2->get_obstacle_back()[o2_counter_old];
            }
            size_t o2_new_size_back = o2_counter_old;

            size_t o1_current_x = o1_x1;
            size_t o1_current_y = o1_y1;
            size_t o1_removing_index = IX(o1_current_x, o1_current_y, o1_z1, Nx, Ny);
            bool o1_end = false;

            size_t o2_current_x = o2_x1;
            size_t o2_current_y = o2_y1;
            size_t o2_removing_index = IX(o2_current_x, o2_current_y, o2_z2, Nx, Ny);
            bool o2_end = false;
            for (; o1_counter_old < o1->get_size_obstacle_front()
                   && o2_counter_old < o2->get_size_obstacle_back() && !o1_end && !o2_end;
                   o1_counter_old++, o2_counter_old++) {
                o1_current_index = o1->get_obstacle_front()[o1_counter_old];
                o2_current_index = o2->get_obstacle_back()[o2_counter_old];
                if (o1_current_index != o1_removing_index) {
                    o1_new.push_back(o1_current_index);
                    o1_new_size_front++;
                } else {
                    o1_current_x++;
                    if (o1_current_x > o1_x2) {
                        o1_current_x = o1_x1;
                        o1_current_y++;
                        if (o1_current_y > o1_y2) {
                            o1_end = true;
                        }
                    }
                    o1_removing_index = IX(o1_current_x, o1_current_y, o1_z1, Nx, Ny);
                }
                if (o2_current_index != o2_removing_index) {
                    o2_new.push_back(o2_current_index);
                    o2_new_size_back++;
                } else {
                    o2_current_x++;
                    if (o2_current_x > o2_x2) {
                        o2_current_x = o2_x1;
                        o2_current_y++;
                        if (o2_current_y > o2_y2) {
                            o2_end = true;
                        }
                    }
                    o2_removing_index = IX(o2_current_x, o2_current_y, o2_z2, Nx, Ny);
                }
            }

            if (!o1_end) {
                for (; o1_counter_old < o1->get_size_obstacle_front() && o1_current_y <= o1_y2;
                       o1_counter_old++) {
                    o1_current_index = o1->get_obstacle_front()[o1_counter_old];
                    if (o1_current_index != o1_removing_index) {
                        o1_new.push_back(o1_current_index);
                        o1_new_size_front++;
                    } else {
                        o1_current_x++;
                        if (o1_current_x > o1_x2) {
                            o1_current_x = o1_x1;
                            o1_current_y++;
                        }
                        o1_removing_index = IX(o1_current_x, o1_current_y, o1_z1, Nx, Ny);
                    }
                }
            }

            if (!o2_end) {
                for (; o2_counter_old < o2->get_size_obstacle_back() && o2_current_y <= o2_y2;
                       o2_counter_old++) {
                    o2_current_index = o2->get_obstacle_back()[o2_counter_old];
                    if (o2_current_index != o2_removing_index) {
                        o2_new.push_back(o2_current_index);
                        o2_new_size_back++;
                    } else {
                        o2_current_x++;
                        if (o2_current_x > o2_x2) {
                            o2_current_x = o2_x1;
                            o2_current_y++;
                        }
                        o2_removing_index = IX(o2_current_x, o2_current_y, o2_z2, Nx, Ny);
                    }
                }
            }

            for (; o1_counter_old < o1->get_size_obstacle_front(); o1_counter_old++) {
                o1_new.push_back(o1->get_obstacle_front()[o1_counter_old]);
                o1_new_size_front++;
            }
            o1_new.resize(o1_new_size_front);

            size_t o1_diff_target = (o1_x2 - o1_x1 + 1) * (o1_y2 - o1_y1 + 1);
            size_t o1_diff_actual = o1->get_size_obstacle_front() - o1_new_size_front;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} front patch: {} -> {} ({}|{})",
                          o1->get_name(), o1->get_size_obstacle_front(), o1_new_size_front,
                          o1_diff_target, o1_diff_actual);
#endif
            size_t *o1_new_data = new size_t[o1_new_size_front];
            std::copy(o1_new.begin(), o1_new.end(), o1_new_data);
            o1->replace_patch(o1_new_data, o1_new_size_front, Patch::FRONT);

            for (; o2_counter_old < o2->get_size_obstacle_back(); o2_counter_old++) {
                o2_new.push_back(o2->get_obstacle_back()[o2_counter_old]);
                o2_new_size_back++;
            }
            o2_new.resize(o2_new_size_back);

            size_t o2_diff_target = (o2_x2 - o2_x1 + 1) * (o2_y2 - o2_y1 + 1);
            size_t o2_diff_actual = o2->get_size_obstacle_back() - o2_new_size_back;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} back patch: {} -> {} ({}|{})",
                          o2->get_name(), o2->get_size_obstacle_back(), o2_new_size_back,
                          o2_diff_target, o2_diff_actual);
#endif
            size_t *o2_new_data = new size_t[o2_new_size_back];
            std::copy(o2_new.begin(), o2_new.end(), o2_new_data);
            o2->replace_patch(o2_new_data, o2_new_size_back, Patch::BACK);
        }
    }
    return overlap;
}

//======================================== set inner cells =========================================
// *************************************************************************************************
/// \brief set inner cells of obstacle in the specified field to the specified value
/// \param f Field where the cells should be changed
/// \param value value to which the cells should be set
// *************************************************************************************************
void Obstacle::set_inner_cells(Field *f, real value) {
    auto Nx = Domain::getInstance()->get_Nx();
    auto Ny = Domain::getInstance()->get_Ny();

    auto data = f->data;
    for (size_t i = m_i1 + 1; i < m_i2; i++) {
        for (size_t j = m_j1 + 1; j < m_j2; j++) {
            for (size_t k = m_k1 + 1; k < m_k2; k++) {
                size_t index = IX(i, j, k, Nx, Ny);
                data[index] = value;
            }
        }
    }
}

//======================================== remove patch ============================================
// *************************************************************************************************
/// \brief remove the specified patch via deleting the array and setting the size to 0
/// \param patch Patch to be removed
// *************************************************************************************************
void Obstacle::remove_patch(Patch patch) {
    switch (patch) {
        case FRONT:
            delete[] m_obstacle_front;
            m_size_obstacle_front = 0;
            break;
        case BACK:
            delete[] m_obstacle_back;
            m_size_obstacle_back = 0;
            break;
        case BOTTOM:
            delete[] m_obstacle_bottom;
            m_size_obstacle_bottom = 0;
            break;
        case TOP:
            delete[] m_obstacle_top;
            m_size_obstacle_top = 0;
            break;
        case LEFT:
            delete[] m_obstacle_left;
            m_size_obstacle_left = 0;
            break;
        case RIGHT:
            delete[] m_obstacle_right;
            m_size_obstacle_right = 0;
            break;
        default:
#ifndef BENCHMARKING
            m_logger->warn("wrong patch: {}", patch);
#endif
            break;
    }
}

//======================================== is corner cell ==========================================
// *************************************************************************************************
/// \brief return whether cell is a corner cell
/// \param i coordinate in x direction
/// \param j coordinate in y direction
/// \param k coordinate in z direction
/// \return true if cell is a corner cell, otherwise false
// *************************************************************************************************
bool Obstacle::is_corner_cell(size_t i, size_t j, size_t k) const {
    return (i == m_i1 || i == m_i2) && (j == m_j1 || j == m_j2) && (k == m_k1 || k == m_k2);
}

//======================================== is edge cell ============================================
// *************************************************************************************************
/// \brief return whether cell is a edge cell
/// \param i coordinate in x direction
/// \param j coordinate in y direction
/// \param k coordinate in z direction
/// \return true if cell is a edge cell, otherwise false
// *************************************************************************************************
bool Obstacle::is_edge_cell(size_t i, size_t j, size_t k) const {
    bool on_x = (i == m_i1 || i == m_i2);
    bool on_y = (j == m_j1 || j == m_j2);
    bool on_z = (k == m_k1 || k == m_k2);

    return (on_x && on_y) || (on_y || on_z) || (on_x && on_z);
}

bool Obstacle::has_overlap(size_t i1, size_t i2, size_t j1, size_t j2, size_t k1, size_t k2) const {
    return has_overlap(m_i1, m_i2, i1, i2) && has_overlap(m_j1, m_j2, j1, j2) && has_overlap(m_k1, m_k2, k1, k2);
}
