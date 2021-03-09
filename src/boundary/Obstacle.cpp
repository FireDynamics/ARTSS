/// \file       Obstacle.cpp
/// \brief      Data class of obstacle object
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Obstacle.h"
#include <algorithm>
#include <utility>
#include <vector>


Obstacle::Obstacle(real x1, real x2, real y1, real y2, real z1, real z2, const std::string &name) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_name = name;
    Domain *domain = Domain::getInstance();

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    m_i1 = get_matching_index(x1, dx, X1) + 1;  // plus 1 for ghost cell
    m_j1 = get_matching_index(y1, dy, Y1) + 1;
    m_k1 = get_matching_index(z1, dz, Z1) + 1;

    m_i2 = get_matching_index(x2, dx, X1);
    m_j2 = get_matching_index(y2, dy, Y1);
    m_k2 = get_matching_index(z2, dz, Z1);

    init(0);
}


Obstacle::Obstacle(size_t coords_i1, size_t coords_j1, size_t coords_k1, size_t coords_i2, size_t coords_j2, size_t coords_k2, size_t level, const std::string& name) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_name = name;
    m_level = level;

    m_i1 = coords_i1;
    m_j1 = coords_j1;
    m_k1 = coords_k1;

    m_i2 = coords_i2;
    m_j2 = coords_j2;
    m_k2 = coords_k2;

    init(level);
}


//======================================== Init ====================================
// ***************************************************************************************
/// \brief  Initialize member variables (arrays)
/// \param  level Multigrid level
// ***************************************************************************************
void Obstacle::init(size_t level) {
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);

    size_t strideX = getStrideX();
    size_t strideY = getStrideY();
    size_t strideZ = getStrideZ();

    m_size_obstacleList = strideX * strideY * strideZ;
    m_obstacleList = new size_t[m_size_obstacleList];

    m_size_obstacleFront = strideY * strideX;
    m_size_obstacleBack = strideY * strideX;
    m_size_obstacleBottom = strideZ * strideX;
    m_size_obstacleTop = strideZ * strideX;
    m_size_obstacleLeft = strideZ * strideY;
    m_size_obstacleRight = strideZ * strideY;
    removeCellsAtBoundary(level);

    m_obstacleFront = new size_t[m_size_obstacleFront];
    m_obstacleBack = new size_t[m_size_obstacleBack];

    m_obstacleTop = new size_t[m_size_obstacleTop];
    m_obstacleBottom = new size_t[m_size_obstacleBottom];

    m_obstacleLeft = new size_t[m_size_obstacleLeft];
    m_obstacleRight = new size_t[m_size_obstacleRight];

    // m_size_obstacleInner = (m_strideX - 2) * (m_strideY - 2) * (m_strideZ - 2);
    // m_obstacleInner = new size_t[m_size_obstacleInner];
    createObstacle(Nx, Ny);

    control();
    printDetails();
}

Obstacle::~Obstacle() {
    delete (m_obstacleList);
    delete (m_obstacleFront);
    delete (m_obstacleBack);
    delete (m_obstacleTop);
    delete (m_obstacleBottom);
    delete (m_obstacleLeft);
    delete (m_obstacleRight);
    // delete (m_obstacleInner);
}

//===================================== Create obstacle ==================================
// ***************************************************************************************
/// \brief  Creates lists of indices of obstacle cells
// ***************************************************************************************
void Obstacle::createObstacle(size_t Nx, size_t Ny) {
    size_t strideX = getStrideX();
    size_t strideY = getStrideY();
    size_t strideZ = getStrideZ();

    size_t counter = 0;
    // fill obstacleList with corresponding indices
    for (size_t k = m_k1; k <= m_k2; ++k) {
        for (size_t j = m_j1; j <= m_j2; ++j) {
            for (size_t i = m_i1; i <= m_i2; ++i) {
                size_t idx = IX(i, j, k, Nx, Ny);
                *(m_obstacleList + counter) = idx;
                counter++;
            }
        }
    }

    // DETAILED OBSTACLE LISTS
    // FRONT and BACK of OBSTACLE
    // fill oFront list with front indices of obstacle and oBack list with back indices of obstacle
    if (m_size_obstacleFront > 0) {
        for (size_t j = 0; j < strideY; ++j) {
            for (size_t i = 0; i < strideX; ++i) {
                size_t index = i + strideX * j;
                size_t idx_front = IX(i, j, 0, strideX, strideY);
                *(m_obstacleFront + index) = m_obstacleList[idx_front];
            }
        }
    }
    if (m_size_obstacleBack > 0) {
        for (size_t j = 0; j < strideY; ++j) {
            for (size_t i = 0; i < strideX; ++i) {
                size_t index = i + strideX * j;
                size_t idx_back = IX(i, j, strideZ - 1, strideX, strideY);
                *(m_obstacleBack + index) = m_obstacleList[idx_back];
            }
        }
    }

    // TOP and BOTTOM of OBSTACLE
    // fill m_obstacleTop list with top indices of obstacle and oBottom list with bottom indices of obstacle
    if (m_size_obstacleBottom > 0) {
        for (size_t k = 0; k < strideZ; ++k) {
            for (size_t i = 0; i < strideX; ++i) {
                size_t index = i + strideX * k;
                size_t idx_bottom = IX(i, 0, k, strideX, strideY);
                *(m_obstacleBottom + index) = m_obstacleList[idx_bottom];
            }
        }
    }
    if (m_size_obstacleTop > 0) {
        for (size_t k = 0; k < strideZ; ++k) {
            for (size_t i = 0; i < strideX; ++i) {
                size_t index = i + strideX * k;
                size_t idx_top = IX(i, strideY - 1, k, strideX, strideY);
                *(m_obstacleTop + index) = m_obstacleList[idx_top];
            }
        }
    }

    // LEFT and RIGHT of OBSTACLE
    // fill oLeft list with left indices of obstacle and oRight list with right indices of obstacle
    if (m_size_obstacleLeft > 0) {
        for (size_t k = 0; k < strideZ; ++k) {
            for (size_t j = 0; j < strideY; ++j) {
                size_t index = j + strideY * k;
                size_t idx_left = IX(0, j, k, strideX, strideY);
                *(m_obstacleLeft + index) = m_obstacleList[idx_left];
            }
        }
    }
    if (m_size_obstacleRight > 0) {
        for (size_t k = 0; k < strideZ; ++k) {
            for (size_t j = 0; j < strideY; ++j) {
                size_t index = j + strideY * k;
                size_t idx_right = IX(strideX - 1, j, k, strideX, strideY);
                *(m_obstacleRight + index) = m_obstacleList[idx_right];
            }
        }
    }

    // // INNER of OBSTACLE
    // // fill oInner list with inner indices of obstacles
    // for (size_t k = 1; k < strideZ - 1; ++k) {
    //     for (size_t j = 1; j < strideY - 1; ++j) {
    //         for (size_t i = 1; i < strideX - 1; ++i) {
    //             size_t index = (i - 1) + (strideX - 2) * (j - 1) + (strideX - 2) * (strideY - 2) * (k - 1);

    //             size_t idx = IX(i, j, k, strideX, strideY);
    //             *(m_obstacleInner + index) = m_obstacle_list[idx];
    //         }
    //     }
    // }
}

//======================================== Print ====================================
// ***************************************************************************************
/// \brief  Print obstacle infos
// ***************************************************************************************
void Obstacle::print() {
#ifndef BENCHMARKING
    size_t strideX = getStrideX();
    size_t strideY = getStrideY();
    size_t strideZ = getStrideZ();

    m_logger->info("-- Obstacle {}", m_name);
    m_logger->info("\t strides (x y z): {} {} {}", strideX, strideY, strideZ);
    m_logger->info("\t size of slices  (Front|Back Bottom|Top Left|Right): {}|{} {}|{} {}|{}", m_size_obstacleFront, m_size_obstacleBack, m_size_obstacleBottom, m_size_obstacleTop, m_size_obstacleLeft, m_size_obstacleRight);
    m_logger->info("\t size of Obstacle: {}", m_size_obstacleList);
    m_logger->info("\t coords (x y z): ({}|{}) ({}|{}) ({}|{})", m_i1, m_i2, m_j1, m_j2, m_k1, m_k2);
#endif
}

//======================================== Print ====================================
// ***************************************************************************************
/// \brief  Print detailed obstacle infos
// ***************************************************************************************
void Obstacle::printDetails() {
#ifndef BENCHMARKING
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(m_level);
    size_t Ny = domain->get_Ny(m_level);
    size_t strideX = getStrideX();
    size_t strideY = getStrideY();
    size_t strideZ = getStrideZ();
    size_t coords_i, coords_j, coords_k;

    m_logger->debug("############### OBSTACLE {} ###############", m_name);
    m_logger->debug("level: {}", m_level);
    m_logger->debug("strides (x y z): {} {} {}", strideX, strideY, strideZ);
    m_logger->debug("size of slices  (Front|Back Bottom|Top Left|Right): {}|{} {}|{} {}|{}", m_size_obstacleFront, m_size_obstacleBack, m_size_obstacleBottom, m_size_obstacleTop, m_size_obstacleLeft, m_size_obstacleRight);
    m_logger->debug("size of Obstacle: {}", m_size_obstacleList);
    m_logger->debug("coords (x y z): ({}|{}) ({}|{}) ({}|{})", m_i1, m_i2, m_j1, m_j2, m_k1, m_k2);

    std::vector<size_t> coords;
    size_t size_front = getSize_obstacleFront();
    if (size_front > 0) {
        m_logger->debug("Front: {} | {}", m_obstacleFront[0], m_obstacleFront[size_front - 1]);

        coords_k = getCoordinateK(m_obstacleFront[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleFront[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleFront[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Front start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacleFront[size_front - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleFront[size_front - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleFront[size_front - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Front end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Front size = 0");
    }

    size_t size_back = getSize_obstacleBack();
    if (size_back > 0) {
        m_logger->debug("Back: {} | {}", m_obstacleBack[0], m_obstacleBack[size_back - 1]);

        coords_k = getCoordinateK(m_obstacleBack[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleBack[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleBack[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Back start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacleBack[size_back - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleBack[size_back - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleBack[size_back - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Back end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Back size = 0");
    }

    size_t size_top = getSize_obstacleTop();
    if (size_top > 0) {
        m_logger->debug("Top: {} | {}", m_obstacleTop[0], m_obstacleTop[size_top - 1]);

        coords_k = getCoordinateK(m_obstacleTop[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleTop[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleTop[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Top start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacleTop[size_top - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleTop[size_top - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleTop[size_top - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Top end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Top size = 0");
    }

    size_t size_bottom = getSize_obstacleBottom();
    if (size_bottom > 0) {
        m_logger->debug("Bottom: {} | {}", m_obstacleBottom[0], m_obstacleBottom[size_bottom - 1]);

        coords_k = getCoordinateK(m_obstacleBottom[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleBottom[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleBottom[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Bottom start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacleBottom[size_bottom - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleBottom[size_bottom - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleBottom[size_bottom - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Bottom end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Bottom size = 0");
    }

    size_t size_left = getSize_obstacleLeft();
    if (size_left > 0) {
        m_logger->debug("Left: {} | {}", m_obstacleLeft[0], m_obstacleLeft[size_left - 1]);

        coords_k = getCoordinateK(m_obstacleLeft[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleLeft[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleLeft[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Left start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacleLeft[size_left - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleLeft[size_left - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleLeft[size_left - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Left end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Left size = 0");
    }

    size_t size_right = getSize_obstacleRight();
    if (size_right > 0) {
        m_logger->debug("Right: {} | {}", m_obstacleRight[0], m_obstacleRight[size_right - 1]);

        coords_k = getCoordinateK(m_obstacleRight[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleRight[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleRight[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Right start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacleRight[size_right - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleRight[size_right - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleRight[size_right - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Right end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Right size = 0");
    }
    m_logger->debug("############### OBSTACLE END ###############");
#endif
}

//======================================== Control ====================================
// ***************************************************************************************
/// \brief  Units test emergency solution
// ***************************************************************************************
void Obstacle::control() {
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(m_level);
    size_t Ny = domain->get_Ny(m_level);

    std::string message;
    for (size_t i = 1; i < m_size_obstacleList; i++) {
        long int diff = static_cast<long int>(m_obstacleList[i]) - static_cast<long int>(m_obstacleList[i - 1]);
        if (diff < 0) {
            message += "sorting error at index " + std::to_string(i - 1) + "|" + std::to_string(i) + " with values " + std::to_string(m_obstacleList[i - 1]) + "|" + std::to_string(m_obstacleList[i]) +
                       "\n";
        }
    }

    size_t start_index = IX(m_i1, m_j1, m_k1, Nx, Ny);
    size_t end_index = IX(m_i2, m_j2, m_k2, Nx, Ny);

    if (m_size_obstacleFront > 0) {
        size_t front_end = IX(m_i2, m_j2, m_k1, Nx, Ny);
        if (start_index != m_obstacleFront[0] || front_end != m_obstacleFront[m_size_obstacleFront - 1]) {
            message += "first or last index of obstacle front list not correct (" + std::to_string(start_index) + "|" + std::to_string(m_obstacleFront[0]) + ")(" + std::to_string(front_end) + "|" +
                       std::to_string(m_obstacleFront[m_size_obstacleFront - 1]) + ")\n";
        }
    }
    if (m_size_obstacleBack > 0) {
        size_t back_start = IX(m_i1, m_j1, m_k2, Nx, Ny);
        if (back_start != m_obstacleBack[0] || end_index != m_obstacleBack[m_size_obstacleBack - 1]) {
            message += "first or last index of obstacle back list not correct (" + std::to_string(back_start) + "|" + std::to_string(m_obstacleBack[0]) + ")(" + std::to_string(end_index) + "|" +
                       std::to_string(m_obstacleBack[m_size_obstacleBack - 1]) + ")\n";
        }
    }
    if (m_size_obstacleBottom > 0) {
        size_t bottom_end = IX(m_i2, m_j1, m_k2, Nx, Ny);
        if (start_index != m_obstacleBottom[0] || bottom_end != m_obstacleBottom[m_size_obstacleBottom - 1]) {
            message += "first or last index of obstacle bottom list not correct (" + std::to_string(start_index) + "|" + std::to_string(m_obstacleBottom[0]) + ")(" + std::to_string(bottom_end) + "|" +
                       std::to_string(m_obstacleBottom[m_size_obstacleBottom - 1]) + ")\n";
        }
    }
    if (m_size_obstacleTop > 0) {
        size_t top_start = IX(m_i1, m_j2, m_k1, Nx, Ny);
        if (top_start != m_obstacleTop[0] || end_index != m_obstacleTop[m_size_obstacleTop - 1]) {
            message += "first or last index of obstacle top list not correct (" + std::to_string(top_start) + "|" + std::to_string(m_obstacleTop[0]) + ")(" + std::to_string(end_index) + "|" +
                       std::to_string(m_obstacleTop[m_size_obstacleTop - 1]) + ")\n";
        }
    }
    if (m_size_obstacleLeft > 0) {
        size_t left_end = IX(m_i1, m_j2, m_k2, Nx, Ny);
        if (start_index != m_obstacleLeft[0] || left_end != m_obstacleLeft[m_size_obstacleLeft - 1]) {
            message += "first or last index of obstacle left list not correct (" + std::to_string(start_index) + "|" + std::to_string(m_obstacleLeft[0]) + ")(" + std::to_string(left_end) + "|" +
                       std::to_string(m_obstacleLeft[m_size_obstacleLeft - 1]) + ")\n";
        }
    }
    if (m_size_obstacleRight > 0) {
        size_t right_start = IX(m_i2, m_j1, m_k1, Nx, Ny);
        if (right_start != m_obstacleRight[0] || end_index != m_obstacleRight[m_size_obstacleRight - 1]) {
            message += "first or last index of obstacle right list not correct (" + std::to_string(right_start) + "|" + std::to_string(m_obstacleRight[0]) + ")(" + std::to_string(end_index) + "|" +
                       std::to_string(m_obstacleRight[m_size_obstacleRight - 1]) + ")\n";
        }
    }
    if (!message.empty()) {
        message = "################ OBSTACLE CONTROL ################\n-- name " + m_name + "\n-- level " + std::to_string(m_level) + "\n" + message + "---------------- OBSTACLE CONTROL END ----------------";
#ifndef BENCHMARKING
        m_logger->warn(message);
#endif
    }
}

//======================================== Is obstacle cell ====================================
// ***************************************************************************************
/// \brief  Check if cell is an obstacle cell
/// \param  i x-coordinate
/// \param  j y-coordinate
/// \param  k z-coordinate
/// \return  bool true if yes false if no
// ***************************************************************************************
bool Obstacle::isObstacleCell(size_t i, size_t j, size_t k) const {
    return m_i1 <= i && i <= m_i2 && m_j1 <= j && j <= m_j2 && m_k1 <= k && k <= m_k2;
}

//======================================== Match grid ====================================
// ***************************************************************************************
/// \brief  Snaps value to grid discretisation
/// \param  obstacle_coordinate Coordinate of obstacle
/// \param  spacing dx/dy/dz
/// \param  start_coordinate X1/Y1/Z1
/// \return real Calculated real grid coordinate
// ***************************************************************************************
real Obstacle::match_grid(real obstacle_coordinate, real spacing, real start_coordinate) {
    return get_matching_index(obstacle_coordinate, spacing, start_coordinate) * spacing + start_coordinate;
}

int Obstacle::get_matching_index(real obstacle_coordinate, real spacing, real start_coordinate) {
    return static_cast<int>(round((-start_coordinate + obstacle_coordinate) / spacing));
}

//======================================== Remove cells at boundary ====================================
// ***************************************************************************************
/// \brief  Remove obstacle patch facing the boundary
/// \param  level Multigrid level
// ***************************************************************************************
void Obstacle::removeCellsAtBoundary(size_t level) {
    Domain *domain = Domain::getInstance();
    if (m_k1 <= domain->get_index_z1(level)) {
        m_size_obstacleFront = 0;
    }
    if (m_k2 >= domain->get_index_z2(level)) {
        m_size_obstacleBack = 0;
    }
    if (m_j1 <= domain->get_index_y1(level)) {
        m_size_obstacleBottom = 0;
    }
    if (m_j2 >= domain->get_index_y2(level)) {
        m_size_obstacleTop = 0;
    }
    if (m_i1 <= domain->get_index_x1(level)) {
        m_size_obstacleLeft = 0;
    }
    if (m_i2 >= domain->get_index_x2(level)) {
        m_size_obstacleRight = 0;
    }
}

/**
 * checks, if coordinates are overlapping
 * @param o1_coord1 starting coordinate of obstacle 1
 * @param o1_coord2 ending coordinate of obstacle 1
 * @param o2_coord1 starting coordinate of obstacle 2
 * @param o2_coord2 starting coordinate of obstacle 2
 * @return
 */
bool Obstacle::hasOverlap(size_t o1_coord1, size_t o1_coord2, size_t o2_coord1, size_t o2_coord2) {
    return o1_coord1 <= o2_coord2 && o1_coord2 >= o2_coord1;
}

void Obstacle::replace_patch(size_t *indices, size_t size, Patch p) {
    switch (p) {
        case FRONT:
            delete[] m_obstacleFront;
            m_obstacleFront = indices;
            m_size_obstacleFront = size;
            break;
        case BACK:
            delete[] m_obstacleBack;
            m_obstacleBack = indices;
            m_size_obstacleBack = size;
            break;
        case BOTTOM:
            delete[] m_obstacleBottom;
            m_obstacleBottom = indices;
            m_size_obstacleBottom = size;
            break;
        case TOP:
            delete[] m_obstacleTop;
            m_obstacleTop = indices;
            m_size_obstacleTop = size;
            break;
        case LEFT:
            delete[] m_obstacleLeft;
            m_obstacleLeft = indices;
            m_size_obstacleLeft = size;
            break;
        case RIGHT:
            delete[] m_obstacleRight;
            m_obstacleRight = indices;
            m_size_obstacleRight = size;
            break;
        case UNKNOWN_PATCH:
            m_logger->warn("wrong patch: {}", p);
    }
}

/**
 * calculate indices of area to be excluded. o1_coordinate == o2_coordinate only if the length of the patch of both obstacles are the same
 * @param o1 Obstacle 1
 * @param o2 Obstacle 2
 * @param o1_coordinate calculated coordinate of obstacle 1
 * @param o2_coordinate calculated coordinate of obstacle 2
 * @param direction X/Y/Z axis
 * @param start true = (i/j/k)1 or false = (i/j/k)2
 */
void Obstacle::calculate_area_index(Obstacle *o1, Obstacle *o2, size_t *o1_coordinate, size_t *o2_coordinate, CoordinateAxis direction, bool start) {
    if (direction == CoordinateAxis::X) {
        if (start) {
            *o1_coordinate = o1->getCoordinates_i1();
            *o2_coordinate = o2->getCoordinates_i1();
            if (o1->getCoordinates_i1() > o2->getCoordinates_i1()) {
                *o2_coordinate = o1->getCoordinates_i1() + 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            } else if (o1->getCoordinates_i1() < o2->getCoordinates_i1()) {
                *o1_coordinate = o2->getCoordinates_i1() + 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            }
        } else {
            *o1_coordinate = o1->getCoordinates_i2();
            *o2_coordinate = o2->getCoordinates_i2();
            if (o1->getCoordinates_i2() < o2->getCoordinates_i2()) {
                *o2_coordinate = o1->getCoordinates_i2() - 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            } else if (o1->getCoordinates_i2() > o2->getCoordinates_i2()) {
                *o1_coordinate = o2->getCoordinates_i2() - 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            }
        }
    }
    if (direction == CoordinateAxis::Y) {
        if (start) {
            *o1_coordinate = o1->getCoordinates_j1();
            *o2_coordinate = o2->getCoordinates_j1();
            if (o1->getCoordinates_j1() > o2->getCoordinates_j1()) {
                *o2_coordinate = o1->getCoordinates_j1() + 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            } else if (o1->getCoordinates_j1() < o2->getCoordinates_j1()) {
                *o1_coordinate = o2->getCoordinates_j1() + 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            }
        } else {
            *o1_coordinate = o1->getCoordinates_j2();
            *o2_coordinate = o2->getCoordinates_j2();
            if (o1->getCoordinates_j2() < o2->getCoordinates_j2()) {
                *o2_coordinate = o1->getCoordinates_j2() - 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            } else if (o1->getCoordinates_j2() > o2->getCoordinates_j2()) {
                *o1_coordinate = o2->getCoordinates_j2() - 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            }
        }
    }
    if (direction == CoordinateAxis::Z) {
        if (start) {
            *o1_coordinate = o1->getCoordinates_k1();
            *o2_coordinate = o2->getCoordinates_k1();
            if (o1->getCoordinates_k1() > o2->getCoordinates_k1()) {
                *o2_coordinate = o1->getCoordinates_k1() + 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            } else if (o1->getCoordinates_k1() < o2->getCoordinates_k1()) {
                *o1_coordinate = o2->getCoordinates_k1() + 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            }
        } else {
            *o1_coordinate = o1->getCoordinates_k2();
            *o2_coordinate = o2->getCoordinates_k2();
            if (o1->getCoordinates_k2() < o2->getCoordinates_k2()) {
                *o2_coordinate = o1->getCoordinates_k2() - 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            } else if (o1->getCoordinates_k2() > o2->getCoordinates_k2()) {
                *o1_coordinate = o2->getCoordinates_k2() - 1;  // do not remove inner edge, can be accessed by SL Advection Solver
            }
        }
    }
}

/**
 * removes boundary cells which are facing the other obstacle
 * @param o1
 * @param o2
 * @return True/False, whether obstacles are facing each other or not
 */
bool Obstacle::remove_circular_constraints(Obstacle *o1, Obstacle *o2) {
    auto domain = Domain::getInstance();
    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();

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

bool Obstacle::circular_constraints_x_direction(Obstacle *o1, Obstacle *o2) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger("Obstacle");
#endif
    bool overlap = false;

    auto domain = Domain::getInstance();
    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();

    if (o1->getCoordinates_i2() + 1 == o2->getCoordinates_i1()) {
        std::swap(o1, o2);
    }

    if (o1->getCoordinates_i1() - 1 == o2->getCoordinates_i2()) {
        bool j_overlap = hasOverlap(o1->getCoordinates_j1(), o1->getCoordinates_j2(), o2->getCoordinates_j1(), o2->getCoordinates_j2());
        bool k_overlap = hasOverlap(o1->getCoordinates_k1(), o1->getCoordinates_k2(), o2->getCoordinates_k1(), o2->getCoordinates_k2());
        if (j_overlap && k_overlap) {
            logger->debug("obstacles are next to each other. Working on {} left side and on {} right side.", o1->get_name(), o2->get_name());
            // another obstacle (o2) at the left side of o1
            overlap = true;
            // calculate coordinates of area which should be removed
            // the area is for both obstacle the same only if there are equally long
            size_t o1_x1 = o1->getCoordinates_i1();
            size_t o2_x2 = o2->getCoordinates_i2();

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
            logger->debug("removing indices in area ({}) ({}|{}) ({}|{}) for {}", o1_x1, o1_y1, o1_y2, o1_z1, o1_z2, o1->get_name());
            logger->debug("removing indices in area ({}) ({}|{}) ({}|{}) for {}", o2_x2, o2_y1, o2_y2, o2_z1, o2_z2, o2->get_name());
#endif

            size_t o1_size_removing_indices = (o1_y2 + 1 - o1_y1) * (o1_z2 + 1 - o1_z1);
            size_t o2_size_removing_indices = (o2_y2 + 1 - o2_y1) * (o2_z2 + 1 - o2_z1);

            std::vector<size_t> o1_new;
            o1_new.reserve(o1->getSize_obstacleLeft());
            std::vector<size_t> o2_new;
            o2_new.reserve(o1->getSize_obstacleRight());

            size_t o1_counter_old = 0;
            size_t o1_smallest_removing_index = IX(o1_x1, o1_y1, o1_z1, Nx, Ny);
            size_t o1_current_index = o1->getObstacleLeft()[o1_counter_old];
            while (o1_current_index < o1_smallest_removing_index) {
                o1_new.push_back(o1_current_index);
                o1_counter_old++;
                o1_current_index = o1->getObstacleLeft()[o1_counter_old];
            }
            size_t o1_new_size_left = o1_counter_old;

            size_t o2_counter_old = 0;
            size_t o2_smallest_removing_index = IX(o2_x2, o2_y1, o2_z1, Nx, Ny);
            size_t o2_current_index = o2->getObstacleRight()[o2_counter_old];
            while (o2_current_index < o2_smallest_removing_index) {
                o2_new.push_back(o2_current_index);
                o2_counter_old++;
                o2_current_index = o2->getObstacleRight()[o2_counter_old];
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
            for (; o1_counter_old < o1->getSize_obstacleLeft() && o2_counter_old < o2->getSize_obstacleRight() && !o1_end && !o2_end; o1_counter_old++, o2_counter_old++) {
                o1_current_index = o1->getObstacleLeft()[o1_counter_old];
                o2_current_index = o2->getObstacleRight()[o2_counter_old];
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
                for (; o1_counter_old < o1->getSize_obstacleLeft() && o1_current_z <= o1_z2; o1_counter_old++) {
                    o1_current_index = o1->getObstacleLeft()[o1_counter_old];
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
                for (; o2_counter_old < o2->getSize_obstacleRight() && o2_current_z <= o2_z2; o2_counter_old++) {
                    o2_current_index = o2->getObstacleRight()[o2_counter_old];
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

            for (; o1_counter_old < o1->getSize_obstacleLeft(); o1_counter_old++) {
                o1_new.push_back(o1->getObstacleLeft()[o1_counter_old]);
                o1_new_size_left++;
            }
            o1_new.resize(o1_new_size_left);

            size_t o1_diff_target = (o1_z2 - o1_z1 + 1) * (o1_y2 - o1_y1 + 1);
            size_t o1_diff_actual = o1->getSize_obstacleLeft() - o1_new_size_left;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} left patch: {} -> {} ({}|{})", o1->get_name(), o1->getSize_obstacleLeft(), o1_new_size_left, o1_diff_target, o1_diff_actual);
#endif
            size_t *o1_new_data = new size_t[o1_new_size_left];
            std::copy(o1_new.begin(), o1_new.end(), o1_new_data);
            o1->replace_patch(o1_new_data, o1_new_size_left, Patch::LEFT);

            for (; o2_counter_old < o2->getSize_obstacleRight(); o2_counter_old++) {
                o2_new.push_back(o2->getObstacleRight()[o2_counter_old]);
                o2_new_size_right++;
            }
            o2_new.resize(o2_new_size_right);

            size_t o2_diff_target = (o2_z2 - o2_z1 + 1) * (o2_y2 - o2_y1 + 1);
            size_t o2_diff_actual = o2->getSize_obstacleRight() - o2_new_size_right;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} right patch: {} -> {} ({}|{})", o2->get_name(), o2->getSize_obstacleRight(), o2_new_size_right, o2_diff_target, o2_diff_actual);
#endif
            size_t *o2_new_data = new size_t[o2_new_size_right];
            std::copy(o2_new.begin(), o2_new.end(), o2_new_data);
            o2->replace_patch(o2_new_data, o2_new_size_right, Patch::RIGHT);
        }
    }
    return overlap;
}

bool Obstacle::circular_constraints_y_direction(Obstacle *o1, Obstacle *o2) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger("Obstacle");
#endif
    bool overlap = false;

    auto domain = Domain::getInstance();
    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();

    if (o1->getCoordinates_j2() + 1 == o2->getCoordinates_j1()) {
        std::swap(o1, o2);
    }

    if (o1->getCoordinates_j1() - 1 == o2->getCoordinates_j2()) {
        bool i_overlap = hasOverlap(o1->getCoordinates_i1(), o1->getCoordinates_i2(), o2->getCoordinates_i1(), o2->getCoordinates_i2());
        bool k_overlap = hasOverlap(o1->getCoordinates_k1(), o1->getCoordinates_k2(), o2->getCoordinates_k1(), o2->getCoordinates_k2());
        if (i_overlap && k_overlap) {
#ifndef BENCHMARKING
            logger->debug("obstacles are next to each other. Working on {} bottom side and on {} top side", o1->get_name(), o2->get_name());
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

            size_t o1_y1 = o1->getCoordinates_j1();
            size_t o2_y2 = o2->getCoordinates_j2();

            size_t o1_z1;
            size_t o2_z1;
            Obstacle::calculate_area_index(o1, o2, &o1_z1, &o2_z1, CoordinateAxis::Z, true);

            size_t o1_z2;
            size_t o2_z2;
            Obstacle::calculate_area_index(o1, o2, &o1_z2, &o2_z2, CoordinateAxis::Z, false);

#ifndef BENCHMARKING
            logger->debug("removing indices in area ({}|{}) ({}) ({}|{}) for {}", o1_x1, o1_x2, o1_y1, o1_z1, o1_z2, o1->get_name());
            logger->debug("removing indices in area ({}|{}) ({}) ({}|{}) for {}", o2_x1, o2_x2, o2_y2, o2_z1, o2_z2, o2->get_name());
#endif

            std::vector<size_t> o1_new;
            o1_new.reserve(o1->getSize_obstacleBottom());
            std::vector<size_t> o2_new;
            o2_new.reserve(o2->getSize_obstacleTop());

            size_t o1_counter_old = 0;
            size_t o1_smallest_removing_index = IX(o1_x1, o1_y1, o1_z1, Nx, Ny);
            size_t o1_current_index = o1->getObstacleBottom()[o1_counter_old];
            while (o1_current_index < o1_smallest_removing_index) {
                o1_new.push_back(o1_current_index);
                o1_counter_old++;
                o1_current_index = o1->getObstacleBottom()[o1_counter_old];
            }
            size_t o1_new_size_bottom = o1_counter_old;

            size_t o2_counter_old = 0;
            size_t o2_smallest_removing_index = IX(o2_x1, o2_y2, o2_z1, Nx, Ny);
            size_t o2_current_index = o2->getObstacleTop()[o2_counter_old];
            while (o2_current_index < o2_smallest_removing_index) {
                o2_new.push_back(o2_current_index);
                o2_counter_old++;
                o2_current_index = o2->getObstacleTop()[o2_counter_old];
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
            for (; o1_counter_old < o1->getSize_obstacleBottom() && o2_counter_old < o2->getSize_obstacleTop() && !o1_end && !o2_end; o1_counter_old++, o2_counter_old++) {
                o1_current_index = o1->getObstacleBottom()[o1_counter_old];
                o2_current_index = o2->getObstacleTop()[o2_counter_old];
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
                for (; o1_counter_old < o1->getSize_obstacleBottom() && o1_current_z <= o1_z2; o1_counter_old++) {
                    o1_current_index = o1->getObstacleBottom()[o1_counter_old];
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
                for (; o2_counter_old < o2->getSize_obstacleTop() && o2_current_z <= o2_z2; o2_counter_old++) {
                    o2_current_index = o2->getObstacleTop()[o2_counter_old];
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

            for (; o1_counter_old < o1->getSize_obstacleBottom(); o1_counter_old++) {
                o1_new.push_back(o1->getObstacleBottom()[o1_counter_old]);
                o1_new_size_bottom++;
            }
            o1_new.resize(o1_new_size_bottom);

            size_t o1_diff_target = (o1_x2 - o1_x1 + 1) * (o1_z2 - o1_z1 + 1);
            size_t o1_diff_actual = o1->getSize_obstacleBottom() - o1_new_size_bottom;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} bottom patch: {} -> {} ({}|{})", o1->get_name(), o1->getSize_obstacleBottom(), o1_new_size_bottom, o1_diff_target, o1_diff_actual);
#endif
            size_t *o1_new_data = new size_t[o1_new_size_bottom];
            std::copy(o1_new.begin(), o1_new.end(), o1_new_data);
            o1->replace_patch(o1_new_data, o1_new_size_bottom, Patch::BOTTOM);

            for (; o2_counter_old < o2->getSize_obstacleTop(); o2_counter_old++) {
                o2_new.push_back(o2->getObstacleTop()[o2_counter_old]);
                o2_new_size_top++;
            }
            o2_new.resize(o2_new_size_top);

            size_t o2_diff_target = (o2_x2 - o2_x1 + 1) * (o2_z2 - o2_z1 + 1);
            size_t o2_diff_actual = o2->getSize_obstacleTop() - o2_new_size_top;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} top patch: {} -> {} ({}|{})", o2->get_name(), o2->getSize_obstacleTop(), o2_new_size_top, o2_diff_target, o2_diff_actual);
#endif
            size_t *o2_new_data = new size_t[o2_new_size_top];
            std::copy(o2_new.begin(), o2_new.end(), o2_new_data);
            o2->replace_patch(o2_new_data, o2_new_size_top, Patch::TOP);
        }
    }
    return overlap;
}

bool Obstacle::circular_constraints_z_direction(Obstacle *o1, Obstacle *o2) {
#ifndef BENCHMARKING
    auto logger = Utility::create_logger("Obstacle");
#endif
    bool overlap = false;

    auto domain = Domain::getInstance();
    auto Nx = domain->get_Nx();
    auto Ny = domain->get_Ny();

    if (o1->getCoordinates_k2() + 1 == o2->getCoordinates_k1()) {
        std::swap(o1, o2);
    }

    if (o1->getCoordinates_k1() - 1 == o2->getCoordinates_k2()) {
        bool i_overlap = hasOverlap(o1->getCoordinates_i1(), o1->getCoordinates_i2(), o2->getCoordinates_i1(), o2->getCoordinates_i2());
        bool j_overlap = hasOverlap(o1->getCoordinates_j1(), o1->getCoordinates_j2(), o2->getCoordinates_j1(), o2->getCoordinates_j2());
        if (i_overlap && j_overlap) {
#ifndef BENCHMARKING
            logger->debug("obstacles are next to each other. Working on {} front side and on {} back side.", o1->get_name(), o2->get_name());
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

            size_t o1_z1 = o1->getCoordinates_k1();
            size_t o2_z2 = o2->getCoordinates_k2();

#ifndef BENCHMARKING
            logger->debug("removing indices in area ({}|{}) ({}|{}) ({}) for {}", o1_x1, o1_x2, o1_y1, o1_y2, o1_z1, o1->get_name());
            logger->debug("removing indices in area ({}|{}) ({}|{}) ({}) for {}", o2_x1, o2_x2, o2_y1, o2_y2, o2_z2, o2->get_name());
#endif

            std::vector<size_t> o1_new;
            o1_new.reserve(o1->getSize_obstacleFront());
            std::vector<size_t> o2_new;
            o2_new.reserve(o1->getSize_obstacleBack());

            size_t o1_counter_old = 0;
            size_t o1_smallest_removing_index = IX(o1_x1, o1_y1, o1_z1, Nx, Ny);
            size_t o1_current_index = o1->getObstacleFront()[o1_counter_old];
            while (o1_current_index < o1_smallest_removing_index) {
                o1_new.push_back(o1_current_index);
                o1_counter_old++;
                o1_current_index = o1->getObstacleFront()[o1_counter_old];
            }
            size_t o1_new_size_front = o1_counter_old;

            size_t o2_counter_old = 0;
            size_t o2_smallest_removing_index = IX(o2_x1, o2_y1, o2_z2, Nx, Ny);
            size_t o2_current_index = o2->getObstacleBack()[o2_counter_old];
            while (o2_current_index < o2_smallest_removing_index) {
                o2_new.push_back(o2_current_index);
                o2_counter_old++;
                o2_current_index = o2->getObstacleBack()[o2_counter_old];
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
            for (; o1_counter_old < o1->getSize_obstacleFront() && o2_counter_old < o2->getSize_obstacleBack() && !o1_end && !o2_end; o1_counter_old++, o2_counter_old++) {
                o1_current_index = o1->getObstacleFront()[o1_counter_old];
                o2_current_index = o2->getObstacleBack()[o2_counter_old];
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
                for (; o1_counter_old < o1->getSize_obstacleFront() && o1_current_y <= o1_y2; o1_counter_old++) {
                    o1_current_index = o1->getObstacleFront()[o1_counter_old];
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
                for (; o2_counter_old < o2->getSize_obstacleBack() && o2_current_y <= o2_y2; o2_counter_old++) {
                    o2_current_index = o2->getObstacleBack()[o2_counter_old];
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

            for (; o1_counter_old < o1->getSize_obstacleFront(); o1_counter_old++) {
                o1_new.push_back(o1->getObstacleFront()[o1_counter_old]);
                o1_new_size_front++;
            }
            o1_new.resize(o1_new_size_front);

            size_t o1_diff_target = (o1_x2 - o1_x1 + 1) * (o1_y2 - o1_y1 + 1);
            size_t o1_diff_actual = o1->getSize_obstacleFront() - o1_new_size_front;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} front patch: {} -> {} ({}|{})", o1->get_name(), o1->getSize_obstacleFront(), o1_new_size_front, o1_diff_target, o1_diff_actual);
#endif
            size_t *o1_new_data = new size_t[o1_new_size_front];
            std::copy(o1_new.begin(), o1_new.end(), o1_new_data);
            o1->replace_patch(o1_new_data, o1_new_size_front, Patch::FRONT);

            for (; o2_counter_old < o2->getSize_obstacleBack(); o2_counter_old++) {
                o2_new.push_back(o2->getObstacleBack()[o2_counter_old]);
                o2_new_size_back++;
            }
            o2_new.resize(o2_new_size_back);

            size_t o2_diff_target = (o2_x2 - o2_x1 + 1) * (o2_y2 - o2_y1 + 1);
            size_t o2_diff_actual = o2->getSize_obstacleBack() - o2_new_size_back;
#ifndef BENCHMARKING
            logger->debug("new size of obstacle {} back patch: {} -> {} ({}|{})", o2->get_name(), o2->getSize_obstacleBack(), o2_new_size_back, o2_diff_target, o2_diff_actual);
#endif
            size_t *o2_new_data = new size_t[o2_new_size_back];
            std::copy(o2_new.begin(), o2_new.end(), o2_new_data);
            o2->replace_patch(o2_new_data, o2_new_size_back, Patch::BACK);
        }
    }
    return overlap;
}

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

void Obstacle::remove_patch(Patch patch) {
    switch (patch) {
        case FRONT:
            delete[] m_obstacleFront;
            m_size_obstacleFront = 0;
            break;
        case BACK:
            delete[] m_obstacleBack;
            m_size_obstacleBack = 0;
            break;
        case BOTTOM:
            delete[] m_obstacleBottom;
            m_size_obstacleBottom = 0;
            break;
        case TOP:
            delete[] m_obstacleTop;
            m_size_obstacleTop = 0;
            break;
        case LEFT:
            delete[] m_obstacleLeft;
            m_size_obstacleLeft = 0;
            break;
        case RIGHT:
            delete[] m_obstacleRight;
            m_size_obstacleRight = 0;
            break;
        case UNKNOWN_PATCH:
            m_logger->warn("wrong patch: {}", patch);
    }
}

bool Obstacle::is_corner_cell(size_t i, size_t j, size_t k) {
    return (i == m_i1 || i == m_i2) && (j == m_j1 || j == m_j2) && (k == m_k1 || k == m_k2);
}

bool Obstacle::is_edge_cell(size_t i, size_t j, size_t k) {
    bool on_x = (i == m_i1 || i == m_i2);
    bool on_y = (j == m_j1 || j == m_j2);
    bool on_z = (k == m_k1 || k == m_k2);

    return (on_x && on_y) || (on_y || on_z) || (on_x && on_z);
}
