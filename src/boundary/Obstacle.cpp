/// \file       Obstacle.cpp
/// \brief      Data class of obstacle object
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include <vector>
#include "Obstacle.h"
#include "../Domain.h"
#include "../utility/Utility.h"

Obstacle::Obstacle(real x1, real x2, real y1, real y2, real z1, real z2) {
    //std::cout << "################ OBSTACLE ################" << std::endl;
    Domain *domain = Domain::getInstance();

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    real rdx = 1. / dx;
    real rdy = 1. / dy;
    real rdz = 1. / dz;

    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

//    std::cout << "obstacle ID ";
//    m_obstacleID = element->IntAttribute("ID");
//    std::cout << m_obstacleID << std::endl;

    //real ox1 = match_grid(x1, dx, X1);
    //real ox2 = match_grid(x2, dx, X1);
    //real oy1 = match_grid(y1, dy, Y1);
    //real oy2 = match_grid(y2, dy, Y1);
    //real oz1 = match_grid(z1, dz, Z1);
    //real oz2 = match_grid(z2, dz, Z1);

    /* ox1 and ox2 as outer boundary
       #########################
       #       #       #       #
       #########################
       ^       ^               ^
      X1      ox1             ox2
                   ^       ^
                  m_i1    m_i2
    */

    m_i1 = get_matching_index(x1, dx, X1) + 1;//((ox1 - X1) * rdx) + 1; //plus 1 for ghost cell
    m_j1 = get_matching_index(y1, dy, Y1) + 1;//((oy1 - Y1) * rdy) + 1;
    m_k1 = get_matching_index(z1, dz, Z1) + 1;//((oz1 - Z1) * rdz) + 1;

    m_i2 = get_matching_index(x2, dx, X1);//((ox2 - X1) * rdx);
    m_j2 = get_matching_index(y2, dy, Y1);//((oy2 - Y1) * rdy);
    m_k2 = get_matching_index(z2, dz, Z1);//((oz2 - Z1) * rdz);

    init(0);
    //std::cout << "---------------- END OBSTACLE ----------------" << std::endl;
}


Obstacle::Obstacle(size_t coords_i1, size_t coords_j1, size_t coords_k1, size_t coords_i2, size_t coords_j2, size_t coords_k2, size_t level) {
    //std::cout << "################ OBSTACLE coarse ################" << std::endl;
    //std::cout << "LEVEL: " << level << std::endl;
    m_level = level;

    m_i1 = coords_i1;
    m_j1 = coords_j1;
    m_k1 = coords_k1;

    m_i2 = coords_i2;
    m_j2 = coords_j2;
    m_k2 = coords_k2;

    init(level);
    //std::cout << "---------------- END OBSTACLE coarse ----------------" << std::endl;
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

    //m_size_obstacleInner = (m_strideX - 2) * (m_strideY - 2) * (m_strideZ - 2);
    //m_obstacleInner = new size_t[m_size_obstacleInner];
    createObstacle(Nx, Ny);

    control();
    //printDetails();
}

Obstacle::~Obstacle() {
    delete (m_obstacleList);
    delete (m_obstacleFront);
    delete (m_obstacleBack);
    delete (m_obstacleTop);
    delete (m_obstacleBottom);
    delete (m_obstacleLeft);
    delete (m_obstacleRight);
    //delete (m_obstacleInner);
}

//======================================== Create obstacle ====================================
// ***************************************************************************************
/// \brief  Creates lists of indices of obstacle cells
// ***************************************************************************************
void Obstacle::createObstacle(size_t Nx, size_t Ny) {
    size_t strideX = getStrideX();
    size_t strideY = getStrideY();
    size_t strideZ = getStrideZ();

    size_t counter = 0;
    //fill obstacleList with corresponding indices
    for (size_t k = m_k1; k <= m_k2; ++k) {
        for (size_t j = m_j1; j <= m_j2; ++j) {
            for (size_t i = m_i1; i <= m_i2; ++i) {
                size_t idx = IX(i, j, k, Nx, Ny);
                *(m_obstacleList + counter) = idx;
                counter++;
            }
        }
    }

//DETAILED OBSTACLE LISTS
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

    // INNER of OBSTACLE
    // fill oInner list with inner indices of obstacles
    //for (size_t k = 1; k < strideZ - 1; ++k) {
    //    for (size_t j = 1; j < strideY - 1; ++j) {
    //        for (size_t i = 1; i < strideX - 1; ++i) {
    //            size_t index = (i - 1) + (strideX - 2) * (j - 1) + (strideX - 2) * (strideY - 2) * (k - 1);

    //            size_t idx = IX(i, j, k, strideX, strideY);
    //            *(m_obstacleInner + index) = m_obstacleList[idx];
    //        }
    //    }
    //}
}

//======================================== Print ====================================
// ***************************************************************************************
/// \brief  Print obstacle infos
// ***************************************************************************************
void Obstacle::print() {
    size_t strideX = getStrideX();
    size_t strideY = getStrideY();
    size_t strideZ = getStrideZ();

    std::cout << "-- Obstacle" << std::endl;
    std::cout << "\t strides (x y z): " << strideX << " " << strideY << " " << strideZ << std::endl;
    std::cout << "\t size of slices  (Front|Back Bottom|Top Left|Right): " << m_size_obstacleFront << "|" << m_size_obstacleBack << " " << m_size_obstacleBottom << "|" <<m_size_obstacleTop << " " << m_size_obstacleLeft << "|" << m_size_obstacleRight << std::endl;
    std::cout << "\t size of Obstacle: " << m_size_obstacleList << std::endl;
    //std::cout << "\t size of inner cells: " << m_size_obstacleInner << std::endl;
    std::cout << "\t coords (x y z): (" << m_i1 << "|" << m_i2 << ")(" << m_j1 << "|" << m_j2 << ")(" << m_k1 << "|" << m_k2 << ")" << std::endl;
}

//======================================== Print ====================================
// ***************************************************************************************
/// \brief  Print detailed obstacle infos
// ***************************************************************************************
void Obstacle::printDetails() {
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(m_level);
    size_t Ny = domain->get_Ny(m_level);

    std::vector<size_t> coords;


    size_t size_front = getSize_obstacleFront();
    if (size_front > 0) {
        std::cout << "Front: " << m_obstacleFront[0] << "|" << m_obstacleFront[size_front - 1] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleFront[0], Nx, Ny);
        std::cout << "Front start: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleFront[size_front - 1], Nx, Ny);
        std::cout << "Front end: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
    } else {
        std::cout << "Front size = 0" << std::endl;
    }

    size_t size_back = getSize_obstacleBack();
    if (size_back > 0) {
        std::cout << "Back: " << m_obstacleBack[0] << "|" << m_obstacleBack[size_back - 1] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleBack[0], Nx, Ny);
        std::cout << "Back start: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleBack[size_back - 1], Nx, Ny);
        std::cout << "Back end: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
    } else {
        std::cout << "Back size = 0" << std::endl;
    }

    size_t size_top = getSize_obstacleTop();
    if (size_top > 0) {
        std::cout << "Top: " << m_obstacleTop[0] << "|" << m_obstacleTop[size_top - 1] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleTop[0], Nx, Ny);
        std::cout << "Top start: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleTop[size_top - 1], Nx, Ny);
        std::cout << "Top end: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
    } else {
        std::cout << "Top size = 0" << std::endl;
    }

    size_t size_bottom = getSize_obstacleBottom();
    if (size_bottom > 0) {
        std::cout << "Bottom: " << m_obstacleBottom[0] << "|" << m_obstacleBottom[size_bottom - 1] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleBottom[0], Nx, Ny);
        std::cout << "Bottom start: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleBottom[size_bottom - 1], Nx, Ny);
        std::cout << "Bottom end: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
    } else {
        std::cout << "Bottom size = 0" << std::endl;
    }

    size_t size_left = getSize_obstacleLeft();
    if (size_left > 0) {
        std::cout << "Left: " << m_obstacleLeft[0] << "|" << m_obstacleLeft[size_left - 1] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleLeft[0], Nx, Ny);
        std::cout << "Left start: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleLeft[size_left - 1], Nx, Ny);
        std::cout << "Left end: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
    } else {
        std::cout << "Left size = 0" << std::endl;
    }

    size_t size_right = getSize_obstacleRight();
    if (size_right > 0) {
        std::cout << "Right: " << m_obstacleRight[0] << "|" << m_obstacleRight[size_right - 1] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleRight[0], Nx, Ny);
        std::cout << "Right start: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleRight[size_right - 1], Nx, Ny);
        std::cout << "Right end: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
    } else {
        std::cout << "Right size = 0" << std::endl;
    }
}

//======================================== Control ====================================
// ***************************************************************************************
/// \brief  Units test emergency solution
// ***************************************************************************************
void Obstacle::control() {
    std::string message;
    for (size_t i = 1; i < m_size_obstacleList; i++) {
        long int diff = static_cast<long int>(m_obstacleList[i]) - static_cast<long int>(m_obstacleList[i - 1]);
        if (diff < 0) {
            message += "sorting error at index " + std::to_string(i - 1) + "|" + std::to_string(i) + " with values " + std::to_string(m_obstacleList[i - 1]) + "|" + std::to_string(m_obstacleList[i]) +
                      "\n";
        }
    }
    if (!message.empty()) {
        message = "################ OBSTACLE CONTROL ################\n-- level " + std::to_string(m_level) + "\n" + message + "---------------- OBSTACLE CONTROL END ----------------";
        std::cout << message << std::endl;
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
/// \return real Calculated real grid coordiante
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
    if (m_k1 <= domain->get_index_z1(level)){
        m_size_obstacleFront = 0;
    }
    if (m_k2 >= domain->get_index_z2(level)){
        m_size_obstacleBack = 0;
    }
    if (m_j1 <= domain->get_index_y1(level)){
        m_size_obstacleBottom = 0;
    }
    if (m_j2 >= domain->get_index_y2(level)){
        m_size_obstacleTop = 0;
    }
    if (m_i1 <= domain->get_index_x1(level)){
        m_size_obstacleLeft = 0;
    }
    if (m_i2 >= domain->get_index_x2(level)){
        m_size_obstacleRight = 0;
    }
}

bool Obstacle::hasOverlap(size_t o1_coord1, size_t o1_coord2, size_t o2_coord1, size_t o2_coord2) {
    return ((o1_coord1 <= o2_coord1 && o2_coord1 <= o1_coord2) || (o1_coord1 <= o2_coord2 && o2_coord2 <= o1_coord2));
}

bool Obstacle::removeCellsFacingAnotherObstacle(Obstacle *o) {
    bool overlap = false;
    if (m_i1 - 1 == o->getCoordinates_i2()) {
        if (hasOverlap(m_j1, m_j2, o->getCoordinates_j1(), o->getCoordinates_j2()) && hasOverlap(m_k1, m_k2, o->getCoordinates_k1(), o->getCoordinates_k2())) {
            // another obstacle at the left side
            overlap = true;
            //remove cells from x = m_i1, y = max(m_j1,o_j1), z = max(m_k1,o_k1) til x = m_i1, y = min(m_j2,o_j2), z = min(m_k2,o_k2)
            size_t x1 = m_i1;
            size_t x2 = m_i1;
            size_t y1 = std::max(m_j1, o->getCoordinates_j1());
            size_t y2 = std::max(m_j2, o->getCoordinates_j2());
            size_t z1 = std::max(m_k1, o->getCoordinates_k1());
            size_t z2 = std::max(m_k2, o->getCoordinates_k2());
            size_t new_front_size = (x2 + 1 - x1) * (y2 + 1 - y2) * (z2 + 1 - z1);
           // delete[]
        }
    }

    if (m_i2 + 1 == o->getCoordinates_i1()) {
        if (hasOverlap(m_j1, m_j2, o->getCoordinates_j1(), o->getCoordinates_j2()) && hasOverlap(m_k1, m_k2, o->getCoordinates_k1(), o->getCoordinates_k2())) {
            // another obstacle at the right side
            overlap = true;
        }
    }

    if (m_j1 - 1 == o->getCoordinates_j2()) {
        if (hasOverlap(m_i1, m_i2, o->getCoordinates_i1(), o->getCoordinates_i2()) && hasOverlap(m_k1, m_k2, o->getCoordinates_k1(), o->getCoordinates_k2())) {
            // another obstacle at the bottom side
            overlap = true;
        }
    }

    if (m_j2 + 1 == o->getCoordinates_j1()) {
        if (hasOverlap(m_i1, m_i2, o->getCoordinates_i1(), o->getCoordinates_i2()) && hasOverlap(m_k1, m_k2, o->getCoordinates_k1(), o->getCoordinates_k2())) {
            // another obstacle at the top side
            overlap = true;
        }
    }

    if (m_k1 - 1 == o->getCoordinates_k2()) {
        if (hasOverlap(m_i1, m_i2, o->getCoordinates_i1(), o->getCoordinates_i2()) && hasOverlap(m_j1, m_j2, o->getCoordinates_j1(), o->getCoordinates_j2())) {
            // another obstacle at the front side
            overlap = true;
        }
    }

    if (m_k2 + 1 == o->getCoordinates_k1()) {
        if (hasOverlap(m_i1, m_i2, o->getCoordinates_i1(), o->getCoordinates_i2()) && hasOverlap(m_j1, m_j2, o->getCoordinates_j1(), o->getCoordinates_j2())) {
            // another obstacle at the back side
            overlap = true;
        }
    }
    return overlap;
}
