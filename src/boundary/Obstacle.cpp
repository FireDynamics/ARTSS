/// \file       Obstacle.cpp
/// \brief      Data class of obstacle object
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Obstacle.h"


Obstacle::Obstacle(real x1, real x2, real y1, real y2, real z1, real z2) :
    m_domain(*(Domain::getInstance())) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    real dx = m_domain.get_dx();
    real dy = m_domain.get_dy();
    real dz = m_domain.get_dz();

    real rdx = 1. / dx;
    real rdy = 1. / dy;
    real rdz = 1. / dz;

    real X1 = m_domain.get_X1();
    real Y1 = m_domain.get_Y1();
    real Z1 = m_domain.get_Z1();

    real ox1 = matchGrid(x1, dx, X1);
    real ox2 = matchGrid(x2, dx, X1);
    real oy1 = matchGrid(y1, dy, Y1);
    real oy2 = matchGrid(y2, dy, Y1);
    real oz1 = matchGrid(z1, dz, Z1);
    real oz2 = matchGrid(z2, dz, Z1);

    /* ox1 and ox2 as outer boundary
       #########################
       #       #       #       #
       #########################
       ^       ^               ^
      X1      ox1             ox2
                   ^       ^
                  m_i1    m_i2
    */

    m_i1 = static_cast<size_t>((ox1 - X1) * rdx + 1);  // plus 1 for ghost cell
    m_j1 = static_cast<size_t>((oy1 - Y1) * rdy + 1);
    m_k1 = static_cast<size_t>((oz1 - Z1) * rdz + 1);

    m_i2 = ((ox2 - X1) * rdx);
    m_j2 = ((oy2 - Y1) * rdy);
    m_k2 = ((oz2 - Z1) * rdz);

    init(0);
}


#ifndef BENCHMARKING
Obstacle::Obstacle(size_t coords_i1, size_t coords_j1, size_t coords_k1, size_t coords_i2, size_t coords_j2, size_t coords_k2, size_t level, std::shared_ptr<spdlog::logger> logger, const Domain &domain) :
    m_i1(coords_i1), m_j1(coords_j1), m_k1(coords_k1),
    m_i2(coords_i2), m_j2(coords_j2), m_k2(coords_k2),
    m_level(level),
    m_domain(domain), m_logger(logger) {
    init(m_level);
}
#endif


Obstacle::Obstacle(size_t coords_i1, size_t coords_j1, size_t coords_k1, size_t coords_i2, size_t coords_j2, size_t coords_k2, size_t level) :
    m_domain(*(Domain::getInstance())) {
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
    size_t Nx = m_domain.get_Nx(level);
    size_t Ny = m_domain.get_Ny(level);

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
#ifndef BENCHMARKING
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
#ifndef BENCHMARKING
    size_t Nx = m_domain.get_Nx(m_level);
    size_t Ny = m_domain.get_Ny(m_level);
    size_t strideX = getStrideX();
    size_t strideY = getStrideY();
    size_t strideZ = getStrideZ();
    size_t coords_i, coords_j, coords_k;

    m_logger->debug("############### OBSTACLE ###############");
    m_logger->debug("level: {}", m_level);
    m_logger->debug("strides (x y z): {} {} {}", strideX, strideY, strideZ);
    m_logger->debug("size of slices  (Front|Back Bottom|Top Left|Right): {}|{} {}|{} {}|{}", m_size_obstacleFront, m_size_obstacleBack, m_size_obstacleBottom, m_size_obstacleTop, m_size_obstacleLeft, m_size_obstacleRight);
    m_logger->debug("size of Obstacle: {}", m_size_obstacleList);
    m_logger->debug("coords (x y z): ({}|{}) ({}|{}) ({}|{})", m_i1, m_i2, m_j1, m_j2, m_k1, m_k2);

    size_t size_front = getSize_obstacleFront();
    if (size_front > 0) {
        m_logger->debug("Front: {} | {}",
                m_obstacleFront[0],
                m_obstacleFront[size_front - 1]);

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
        m_logger->debug("Back: {} | {}",
                m_obstacleBack[0],
                m_obstacleBack[size_back-1]);

        coords_k = getCoordinateK(m_obstacleBack[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleBack[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleBack[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Back start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacleBack[size_front - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleBack[size_front - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleBack[size_front - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Back end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Back size = 0");
    }

    size_t size_top = getSize_obstacleTop();
    if (size_top > 0) {
        m_logger->debug("Top: {} | {}",
                m_obstacleTop[0], m_obstacleTop[size_top-1]);

        coords_k = getCoordinateK(m_obstacleTop[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleTop[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleTop[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Top start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacleTop[size_front - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleTop[size_front - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleTop[size_front - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Top end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Top size = 0");
    }

    size_t size_bottom = getSize_obstacleBottom();
    if (size_bottom > 0) {
        m_logger->debug("Bottom: {} | {}",
                m_obstacleBottom[0], m_obstacleBottom[size_bottom-1]);

        coords_k = getCoordinateK(m_obstacleBottom[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleBottom[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleBottom[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Bottom start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacleBottom[size_front - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleBottom[size_front - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleBottom[size_front - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Bottom end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Bottom size = 0");
    }

    size_t size_left = getSize_obstacleLeft();
    if (size_left > 0) {
        m_logger->debug("Left: {} | {}",
                m_obstacleLeft[0], m_obstacleLeft[size_left-1]);

        coords_k = getCoordinateK(m_obstacleLeft[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleLeft[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleLeft[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Left start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacleLeft[size_front - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleLeft[size_front - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleLeft[size_front - 1], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Left end: {}|{}|{}", coords_i, coords_j, coords_k);
    } else {
        m_logger->debug("Left size = 0");
    }

    size_t size_right = getSize_obstacleRight();
    if (size_right > 0) {
        m_logger->debug("Right: {} | {}",
                m_obstacleRight[0], m_obstacleRight[size_right-1]);

        coords_k = getCoordinateK(m_obstacleRight[0], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleRight[0], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleRight[0], Nx, Ny, coords_j, coords_k);
        m_logger->debug("Right start: {}|{}|{}", coords_i, coords_j, coords_k);

        coords_k = getCoordinateK(m_obstacleRight[size_front - 1], Nx, Ny);
        coords_j = getCoordinateJ(m_obstacleRight[size_front - 1], Nx, Ny, coords_k);
        coords_i = getCoordinateI(m_obstacleRight[size_front - 1], Nx, Ny, coords_j, coords_k);
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
bool Obstacle::isObstacleCell(size_t i, size_t j, size_t k) const {
    size_t i2 = getCoordinates_i2();
    size_t j2 = getCoordinates_j2();
    size_t k2 = getCoordinates_k2();
    return m_i1 <= i && i <= i2 && m_j1 <= j && j <= j2 && m_k1 <= k && k <= k2;
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
    const auto point_i1 = getCoordinates_i1();
    const auto point_i2 = getCoordinates_i2();
    const auto point_j1 = getCoordinates_j1();
    const auto point_j2 = getCoordinates_j2();
    const auto point_k1 = getCoordinates_k1();
    const auto point_k2 = getCoordinates_k2();
    real indeces[6] = {static_cast<real>(point_i1), static_cast<real>(point_i2),
                       static_cast<real>(point_j1), static_cast<real>(point_j2),
                       static_cast<real>(point_k1), static_cast<real>(point_k2)};

    if (isObstacleCell(i, j, k))
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
/// \param  obstacleCoordinate Coordinate of obstacle
/// \param  spacing dx/dy/dz
/// \param  startCoordinate X1/Y1/Z1
/// \return real Calculated real grid coordinate
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
    size_t i2 = getCoordinates_i2();
    size_t j2 = getCoordinates_j2();
    size_t k2 = getCoordinates_k2();

    if (m_k1 <= m_domain.get_index_z1(level)){
        m_size_obstacleFront = 0;
    }
    if (k2 >= m_domain.get_index_z2(level)){
        m_size_obstacleBack = 0;
    }
    if (m_j1 <= m_domain.get_index_y1(level)){
        m_size_obstacleBottom = 0;
    }
    if (j2 >= m_domain.get_index_y2(level)){
        m_size_obstacleTop = 0;
    }
    if (m_i1 <= m_domain.get_index_x1(level)){
        m_size_obstacleLeft = 0;
    }
    if (i2 >= m_domain.get_index_x2(level)){
        m_size_obstacleRight = 0;
    }
}
