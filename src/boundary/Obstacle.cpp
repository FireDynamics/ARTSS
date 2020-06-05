/// \file 		Obstacle.h
/// \brief 		Data class of obstacle object
/// \date 		Oct 01, 2019
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include <vector>
#include "Obstacle.h"
#include "../Domain.h"
#include "../utility/Utility.h"

Obstacle::Obstacle(real x1, real x2, real y1, real y2, real z1, real z2) {
    Domain *domain = Domain::getInstance();

    real dx = domain->Getdx();
    real dy = domain->Getdy();
    real dz = domain->Getdz();

    real rdx = 1. / dx;
    real rdy = 1. / dy;
    real rdz = 1. / dz;

    real X1 = domain->GetX1();
    real Y1 = domain->GetY1();
    real Z1 = domain->GetZ1();

    real ox1 = matchGrid(x1, dx, X1);
    real ox2 = matchGrid(x2, dx, X1);
    real oy1 = matchGrid(y1, dy, Y1);
    real oy2 = matchGrid(y2, dy, Y1);
    real oz1 = matchGrid(z1, dz, Z1);
    real oz2 = matchGrid(z2, dz, Z1);
    real lox = fabs(ox2 - ox1);
    real loy = fabs(oy2 - oy1);
    real loz = fabs(oz2 - oz1);
    m_strideX = static_cast<size_t>(round(lox * rdx));
    m_strideY = static_cast<size_t>(round(loy * rdy));
    m_strideZ = static_cast<size_t>(round(loz * rdz));

    m_i1 = static_cast<size_t>((ox1 - X1) * rdx + 1);  // plus 1 for ghost cell
    m_j1 = static_cast<size_t>((oy1 - Y1) * rdy + 1);
    m_k1 = static_cast<size_t>((oz1 - Z1) * rdz + 1);

    init(0);
}


Obstacle::Obstacle(size_t coords_i1, size_t coords_j1, size_t coords_k1, size_t coords_i2, size_t coords_j2, size_t coords_k2, size_t level) {
    m_level = level;
    m_strideX = coords_i2 - coords_i1 + 1;
    m_strideY = coords_j2 - coords_j1 + 1;
    m_strideZ = coords_k2 - coords_k1 + 1;


    m_i1 = coords_i1;
    m_j1 = coords_j1;
    m_k1 = coords_k1;

    init(level);
}


//======================================== Init ====================================
// ***************************************************************************************
/// \brief  Initialize member variables (arrays)
/// \param  level Multigrid level
// ***************************************************************************************
void Obstacle::init(size_t level) {
#ifndef PROFILING
    m_logger = Utility::createLogger(typeid(this).name());
#endif
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->GetNx(level);
    size_t Ny = domain->GetNy(level);
    m_size_obstacleList = m_strideX * m_strideY * m_strideZ;
    m_obstacleList = new size_t[m_size_obstacleList];

    m_size_obstacleFront = m_strideY * m_strideX;
    m_size_obstacleBack = m_strideY * m_strideX;
    m_size_obstacleBottom = m_strideZ * m_strideX;
    m_size_obstacleTop = m_strideZ * m_strideX;
    m_size_obstacleLeft = m_strideZ * m_strideY;
    m_size_obstacleRight = m_strideZ * m_strideY;
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
    size_t i2 = getCoordinates_i2();
    size_t j2 = getCoordinates_j2();
    size_t k2 = getCoordinates_k2();

    size_t counter = 0;
    //fill obstacleList with corresponding indices
    for (size_t k = m_k1; k <= k2; ++k) {
        for (size_t j = m_j1; j <= j2; ++j) {
            for (size_t i = m_i1; i <= i2; ++i) {
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
        for (size_t j = 0; j < m_strideY; ++j) {
            for (size_t i = 0; i < m_strideX; ++i) {
                size_t index = i + m_strideX * j;
                size_t idx_front = IX(i, j, 0, m_strideX, m_strideY);
                *(m_obstacleFront + index) = m_obstacleList[idx_front];
            }
        }
    }
    if (m_size_obstacleBack > 0) {
        for (size_t j = 0; j < m_strideY; ++j) {
            for (size_t i = 0; i < m_strideX; ++i) {
                size_t index = i + m_strideX * j;
                size_t idx_back = IX(i, j, m_strideZ - 1, m_strideX, m_strideY);
                *(m_obstacleBack + index) = m_obstacleList[idx_back];
            }
        }
    }

    // TOP and BOTTOM of OBSTACLE
    // fill m_obstacleTop list with top indices of obstacle and oBottom list with bottom indices of obstacle
    if (m_size_obstacleBottom > 0) {
        for (size_t k = 0; k < m_strideZ; ++k) {
            for (size_t i = 0; i < m_strideX; ++i) {
                size_t index = i + m_strideX * k;
                size_t idx_bottom = IX(i, 0, k, m_strideX, m_strideY);
                *(m_obstacleBottom + index) = m_obstacleList[idx_bottom];
            }
        }
    }
    if (m_size_obstacleTop > 0) {
        for (size_t k = 0; k < m_strideZ; ++k) {
            for (size_t i = 0; i < m_strideX; ++i) {
                size_t index = i + m_strideX * k;
                size_t idx_top = IX(i, m_strideY - 1, k, m_strideX, m_strideY);
                *(m_obstacleTop + index) = m_obstacleList[idx_top];
            }
        }
    }

    // LEFT and RIGHT of OBSTACLE
    // fill oLeft list with left indices of obstacle and oRight list with right indices of obstacle
    if (m_size_obstacleLeft > 0) {
        for (size_t k = 0; k < m_strideZ; ++k) {
            for (size_t j = 0; j < m_strideY; ++j) {
                size_t index = j + m_strideY * k;
                size_t idx_left = IX(0, j, k, m_strideX, m_strideY);
                *(m_obstacleLeft + index) = m_obstacleList[idx_left];
            }
        }
    }
    if (m_size_obstacleRight > 0) {
        for (size_t k = 0; k < m_strideZ; ++k) {
            for (size_t j = 0; j < m_strideY; ++j) {
                size_t index = j + m_strideY * k;
                size_t idx_right = IX(m_strideX - 1, j, k, m_strideX, m_strideY);
                *(m_obstacleRight + index) = m_obstacleList[idx_right];
            }
        }
    }

    //// INNER of OBSTACLE
    //// fill oInner list with inner indices of obstacles
    //for (size_t k = 1; k < m_strideZ - 1; ++k) {
    //    for (size_t j = 1; j < m_strideY - 1; ++j) {
    //        for (size_t i = 1; i < m_strideX - 1; ++i) {
    //            size_t index = (i - 1) + (m_strideX - 2) * (j - 1) + (m_strideX - 2) * (m_strideY - 2) * (k - 1);

    //            size_t idx = IX(i, j, k, m_strideX, m_strideY);
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
    size_t i2 = getCoordinates_i2();
    size_t j2 = getCoordinates_j2();
    size_t k2 = getCoordinates_k2();
    m_logger->info("-- Obstacle");
    m_logger->info("\t strides (x y z): {} {} {}", m_strideX, m_strideY, m_strideZ);
    m_logger->info("\t size of slices  (Front|Back Bottom|Top Left|Right): {}|{} {}|{} {}|{}", m_size_obstacleFront, m_size_obstacleBack, m_size_obstacleBottom, m_size_obstacleTop, m_size_obstacleLeft, m_size_obstacleRight);
    m_logger->info("\t size of Obstacle: {}", m_size_obstacleList);
    m_logger->info("\t coords (x y z): ({}|{}) ({}|{}) ({}|{})", m_i1, i2, m_j1, j2, m_k1, k2);
}

//======================================== Print ====================================
// ***************************************************************************************
/// \brief  Print detailed obstacle infos
// ***************************************************************************************
void Obstacle::printDetails(){
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->GetNx(m_level);
    size_t Ny = domain->GetNy(m_level);

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
        std::cout << "Back: " << m_obstacleBack[0] << "|" << m_obstacleBack[size_back-1] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleBack[0], Nx, Ny);
        std::cout << "Back start: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleBack[size_back-1], Nx, Ny);
        std::cout << "Back end: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
    } else {
        std::cout << "Back size = 0" << std::endl;
    }

    size_t size_top = getSize_obstacleTop();
    if (size_top > 0) {
        std::cout << "Top: " << m_obstacleTop[0] << "|" << m_obstacleTop[size_top-1] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleTop[0], Nx, Ny);
        std::cout << "Top start: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleTop[size_top-1], Nx, Ny);
        std::cout << "Top end: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
    } else {
        std::cout << "Top size = 0" << std::endl;
    }

    size_t size_bottom = getSize_obstacleBottom();
    if (size_bottom > 0) {
        std::cout << "Bottom: " << m_obstacleBottom[0] << "|" << m_obstacleBottom[size_bottom-1] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleBottom[0], Nx, Ny);
        std::cout << "Bottom start: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleBottom[size_bottom - 1], Nx, Ny);
        std::cout << "Bottom end: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
    }else{
        std::cout << "Bottom size = 0" << std::endl;
    }

    size_t size_left = getSize_obstacleLeft();
    if (size_left > 0) {
        std::cout << "Left: " << m_obstacleLeft[0] << "|" << m_obstacleLeft[size_left-1] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleLeft[0], Nx, Ny);
        std::cout << "Left start: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleLeft[size_left-1], Nx, Ny);
        std::cout << "Left end: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
    } else {
        std::cout << "Left size = 0" << std::endl;
    }

    size_t size_right = getSize_obstacleRight();
    if (size_right > 0) {
        std::cout << "Right: " << m_obstacleRight[0] << "|" << m_obstacleRight[size_right-1] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleRight[0], Nx, Ny);
        std::cout << "Right start: " << coords[0] << "|" << coords[1] << "|" << coords[2] << std::endl;
        coords = Utility::coordinateFromLinearIndex(m_obstacleRight[size_right-1], Nx, Ny);
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
   // size_t all_cells = m_size_obstacleSliceX * 2 - 4 * m_strideX + m_size_obstacleSliceY * 2 - 4 * m_strideY + m_size_obstacleSliceZ * 2 - 4 * (m_strideZ - 2) + m_size_obstacleInner;
   // if (m_size_obstacleList != all_cells) {
   //     std::cout << "list size of obstacle does not match the size of its parts. Obstacle List: " << m_size_obstacleList << " all cells: " << all_cells << " inner: " << m_size_obstacleInner << std::endl;
   //     std::cout << "FRONT/BACK: " << m_size_obstacleSliceZ << " BOTTOM/TOP: " << m_size_obstacleSliceY << " LEFT/RIGHT: " << m_size_obstacleSliceX << std::endl;
   // }
}

//======================================== Is obstacle cell ====================================
// ***************************************************************************************
/// \brief  Check if cell is an obstacle cell
/// \param  i x-coordinate
/// \param  j y-coordinate
/// \param  k z-coordinate
/// \return  bool true if yes false if no
// ***************************************************************************************
bool Obstacle::isObstacleCell(size_t i, size_t j, size_t k) {
    size_t i2 = getCoordinates_i2();
    size_t j2 = getCoordinates_j2();
    size_t k2 = getCoordinates_k2();
    return m_i1 <= i && i <= i2 && m_j1 <= j && j <= j2 && m_k1 <= k && k <= k2;
}

size_t Obstacle::getCoordinates_i2() {
    return m_i1 + m_strideX - 1;
}

size_t Obstacle::getCoordinates_j2() {
    return m_j1 + m_strideY - 1;
}

size_t Obstacle::getCoordinates_k2() {
    return m_k1 + m_strideZ - 1;
}

//======================================== Match grid ====================================
// ***************************************************************************************
/// \brief  Snaps value to grid discretization
/// \param  obstacleCoordinate Coordinate of obstacle
/// \param  spacing dx/dy/dz
/// \param  startCoordinate X1/Y1/Z1
/// \return real Calculated real grid coordiante
// ***************************************************************************************
real Obstacle::matchGrid(double obstacleCoordinate, real spacing, real startCoordinate) {
    return std::round((-startCoordinate + obstacleCoordinate) / spacing) * spacing + startCoordinate;
}

//======================================== Remove cells at boundary ====================================
// ***************************************************************************************
/// \brief  Remove obstacle patch facing the boundary
/// \param  level Multigrid level
// ***************************************************************************************
void Obstacle::removeCellsAtBoundary(size_t level) {
    Domain *domain = Domain::getInstance();
    size_t i2 = getCoordinates_i2();
    size_t j2 = getCoordinates_j2();
    size_t k2 = getCoordinates_k2();

    if (m_k1 <= domain->GetIndexz1(level)){
        m_size_obstacleFront = 0;
    }
    if (k2 >= domain->GetIndexz2(level)){
        m_size_obstacleBack = 0;
    }
    if (m_j1 <= domain->GetIndexy1(level)){
        m_size_obstacleBottom = 0;
    }
    if (j2 >= domain->GetIndexy2(level)){
        m_size_obstacleTop = 0;
    }
    if (m_i1 <= domain->GetIndexx1(level)){
        m_size_obstacleLeft = 0;
    }
    if (i2 >= domain->GetIndexx2(level)){
        m_size_obstacleRight = 0;
    }
}
