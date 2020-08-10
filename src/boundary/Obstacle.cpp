/// \file       Obstacle.cpp
/// \brief      Data class of obstacle object
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Obstacle.h"


Obstacle::Obstacle(real x1, real x2, real y1, real y2, real z1, real z2) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
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

    real ox1 = matchGrid(x1, dx, X1);
    real ox2 = matchGrid(x2, dx, X1);
    real oy1 = matchGrid(y1, dy, Y1);
    real oy2 = matchGrid(y2, dy, Y1);
    real oz1 = matchGrid(z1, dz, Z1);
    real oz2 = matchGrid(z2, dz, Z1);

    m_i1 = static_cast<size_t>((ox1 - X1) * rdx + 1);  // plus 1 for ghost cell
    m_j1 = static_cast<size_t>((oy1 - Y1) * rdy + 1);
    m_k1 = static_cast<size_t>((oz1 - Z1) * rdz + 1);

    m_i2 = ((ox2 - X1) * rdx);
    m_j2 = ((oy2 - Y1) * rdy);
    m_k2 = ((oz2 - Z1) * rdz);

    init(0);
}


Obstacle::Obstacle(size_t coords_i1, size_t coords_j1, size_t coords_k1, size_t coords_i2, size_t coords_j2, size_t coords_k2, size_t level) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
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
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
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

    //// INNER of OBSTACLE
    //// fill oInner list with inner indices of obstacles
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
#ifdef BENCHMARKING
    return;
#else
    size_t strideX = getStrideX();
    size_t strideY = getStrideY();
    size_t strideZ = getStrideZ();

    m_logger->info("-- Obstacle");
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
void Obstacle::printDetails(){
#ifdef BENCHMARKING
    return;
#else
    Domain *domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(m_level);
    size_t Ny = domain->get_Ny(m_level);

    std::vector<size_t> coords;


    size_t size_front = getSize_obstacleFront();
    if (size_front > 0) {
        m_logger->info("Front: {} | {}",
                m_obstacleFront[0],
                m_obstacleFront[size_front - 1]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleFront[0], Nx, Ny);
        m_logger->info("Front start: {}|{}|{}",
                coords[0], coords[1], coords[2]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleFront[size_front - 1], Nx, Ny);
        m_logger->info("Front end: {}|{}|{}",
                coords[0], coords[1], coords[2]);
    } else {
        m_logger->info("Front size = 0");
    }

    size_t size_back = getSize_obstacleBack();
    if (size_back > 0) {
        m_logger->info("Back: {} | {}",
                m_obstacleBack[0],
                m_obstacleBack[size_back-1]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleBack[0], Nx, Ny);
        m_logger->info("Back start: {}|{}|{}",
                coords[0], coords[1], coords[2]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleBack[size_back-1], Nx, Ny);
        m_logger->info("Back end: {}|{}|{}",
                coords[0], coords[1], coords[2]);
    } else {
        m_logger->info("Back size = 0");
    }

    size_t size_top = getSize_obstacleTop();
    if (size_top > 0) {
        m_logger->info("Top: {} | {}",
                m_obstacleTop[0], m_obstacleTop[size_top-1]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleTop[0], Nx, Ny);
        m_logger->info("Top start: {}|{}|{}",
                coords[0], coords[1], coords[2]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleTop[size_top-1], Nx, Ny);
        m_logger->info("Top end: {}|{}|{}",
                coords[0], coords[1], coords[2]);
    } else {
        m_logger->info("Top size = 0");
    }

    size_t size_bottom = getSize_obstacleBottom();
    if (size_bottom > 0) {
        m_logger->info("Bottom: {} | {}",
                m_obstacleBottom[0], m_obstacleBottom[size_bottom-1]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleBottom[0], Nx, Ny);
        m_logger->info("Bottom start: {}|{}|{}",
                coords[0], coords[1], coords[2]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleBottom[size_bottom - 1], Nx, Ny);
        m_logger->info("Bottom end: {}|{}|{}",
                coords[0], coords[1], coords[2]);
    } else {
        m_logger->info("Bottom size = 0");
    }

    size_t size_left = getSize_obstacleLeft();
    if (size_left > 0) {
        m_logger->info("Left: {} | {}",
                m_obstacleLeft[0], m_obstacleLeft[size_left-1]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleLeft[0], Nx, Ny);
        m_logger->info("Left start: {}|{}|{}",
                coords[0], coords[1], coords[2]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleLeft[size_left-1], Nx, Ny);
        m_logger->info("Left end: {}|{}|{}",
                coords[0], coords[1], coords[2]);
    } else {
        m_logger->info("Left size = 0");
    }

    size_t size_right = getSize_obstacleRight();
    if (size_right > 0) {
        m_logger->info("Right: {} | {}",
                m_obstacleRight[0], m_obstacleRight[size_right-1]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleRight[0], Nx, Ny);
        m_logger->info("Right start: {}|{}|{}",
                coords[0], coords[1], coords[2]);

        coords = Utility::coordinateFromLinearIndex(m_obstacleRight[size_right-1], Nx, Ny);
        m_logger->info("Right end: {}|{}|{}",
                coords[0], coords[1], coords[2]);
    } else {
        m_logger->info("Right size = 0");
    }
#endif
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

    if (m_k1 <= domain->get_index_z1(level)){
        m_size_obstacleFront = 0;
    }
    if (k2 >= domain->get_index_z2(level)){
        m_size_obstacleBack = 0;
    }
    if (m_j1 <= domain->get_index_y1(level)){
        m_size_obstacleBottom = 0;
    }
    if (j2 >= domain->get_index_y2(level)){
        m_size_obstacleTop = 0;
    }
    if (m_i1 <= domain->get_index_x1(level)){
        m_size_obstacleLeft = 0;
    }
    if (i2 >= domain->get_index_x2(level)){
        m_size_obstacleRight = 0;
    }
}
