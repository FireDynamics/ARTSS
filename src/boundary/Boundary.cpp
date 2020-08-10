/// \file 		Boundary.cpp
/// \brief 		Data class of boundary object
/// \date 		Oct 01, 2019
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>

#include "Boundary.h"
#include "../Field.h"
#include "../Domain.h"
#include "../utility/GlobalMacrosTypes.h"

Boundary::Boundary(size_t level) {
    m_level = level;
    init(0);
    innerCells();

#ifndef BENCHMARKING
    // print(0);
    control(0);
#endif
}

Boundary::Boundary(Obstacle **obstacleList, size_t numberOfObstacles, size_t size_obstacles, size_t level) {
    m_level = level;
    init(size_obstacles);
    innerCells(obstacleList, numberOfObstacles);

#ifndef BENCHMARKING
    // print(size_obstacles);
    control(size_obstacles);
#endif
}


//======================================== Init ====================================
// ***************************************************************************************
/// \brief  Initialize member variables (arrays)
/// \param  size_obstacles Amount of obstacle cells
// ***************************************************************************************
void Boundary::init(size_t size_obstacles){
    auto domain = Domain::getInstance();

    const size_t nx = domain->get_nx(m_level);
    const size_t ny = domain->get_ny(m_level);
    const size_t nz = domain->get_nz(m_level);

    m_size_boundaryList = 2 * nx * ny + 2 * (nz - 2) * (ny - 2) + 2 * (nz - 2) * nx;
    m_boundaryList = new size_t[m_size_boundaryList];

    m_size_boundaryFront = ny * nx;
    m_size_boundaryBack = ny * nx;
    m_boundaryFront = new size_t[m_size_boundaryFront];
    m_boundaryBack = new size_t[m_size_boundaryBack];

    m_size_boundaryTop = nz * nx;
    m_size_boundaryBottom = nz * nx;
    m_boundaryBottom = new size_t[m_size_boundaryBottom];
    m_boundaryTop = new size_t[m_size_boundaryTop];

    m_size_boundaryLeft = nz * ny;
    m_size_boundaryRight = nz * ny;
    m_boundaryLeft = new size_t[m_size_boundaryLeft];
    m_boundaryRight = new size_t[m_size_boundaryRight];

    m_size_innerList = (nx - 2) * (ny - 2) * (nz - 2) - size_obstacles;
    m_innerList = new size_t[m_size_innerList];

    boundaryCells();
}


//======================================== Print ====================================
// ***************************************************************************************
/// \brief  Prints boundary infos
/// \param  size_obstacles Amount of obstacle cells
// ***************************************************************************************
void Boundary::print(size_t size_obstacles) {
#ifdef BENCHMARKING
    return;
#else
    auto m_logger = Utility::create_logger(typeid(Boundary).name());
    m_logger->info("################ BOUNDARY ################");
    m_logger->info("list size of bList: {}", m_size_boundaryList);
    m_logger->info("Boundary starts with {} and ends with {}", *(m_boundaryList + 0),
                                                             *(m_boundaryList + m_size_boundaryList - 1));
    m_logger->info("list size of size_z: {}", m_size_boundaryFront);
    m_logger->info("Front starts with {} and ends with {}", *(m_boundaryFront + 0),
                                                          *(m_boundaryList + m_size_boundaryList - 1));
    m_logger->info("Back starts with {} and ends with {}", *(m_boundaryBack + 0),
                                                         *(m_boundaryList + m_size_boundaryList - 1));
    m_logger->info("list size of size_y: ", m_size_boundaryBottom);
    m_logger->info("Bottom starts with {} and ends with {}", *(m_boundaryBottom + 0),
                                                           *(m_boundaryBottom + m_size_boundaryBottom - 1));
    m_logger->info("Top starts with {} and ends with {}", *(m_boundaryTop + 0),
                                                        *(m_boundaryTop + m_size_boundaryTop - 1));
    m_logger->info("list size of size_x: ", m_size_boundaryLeft);
    m_logger->info("Left starts with {} and ends with {}", *(m_boundaryLeft + 0),
                                                         *(m_boundaryLeft + m_size_boundaryLeft - 1));
    m_logger->info("Right starts with {} and ends with {}", *(m_boundaryRight + 0),
                                                          *(m_boundaryRight + m_size_boundaryRight - 1));
    m_logger->info("list size of innerList: {} obstacle size: {}", m_size_innerList,
                                                                 size_obstacles);
    m_logger->info("Inner starts with {} and ends with {}", *(m_innerList + 0),
                                                          *(m_innerList + m_size_innerList - 1));
    m_logger->info("--------------- END BOUNDARY ---------------");
#endif
}

//======================================== Control ====================================
// ***************************************************************************************
/// \brief  Units test emergency solution
/// \param  size_obstacles Amount of obstacle cells
// ***************************************************************************************
void Boundary::control(size_t size_obstacles) {
    // TODO(n16h7): clean up
    std::string message;
    Domain* domain = Domain::getInstance();
    size_t nx = domain->get_nx(m_level);
    size_t ny = domain->get_ny(m_level);
    size_t nz = domain->get_nz(m_level);
    size_t Nx = domain->get_Nx(m_level);
    size_t Ny = domain->get_Ny(m_level);
    size_t all_cells =  m_size_boundaryFront + m_size_boundaryBack + m_size_boundaryBottom + m_size_boundaryTop + m_size_boundaryLeft + m_size_boundaryRight;
    size_t duplicates = 4 * nx + 4 * (ny - 2) + 4 * (nz - 2) + 8;
    if ( m_size_boundaryList != all_cells - duplicates){
        message = message + "list size of all boundary cells does not fit with sum of it parts. Boundary List: " + std::to_string(m_size_boundaryList) + " sum: " + std::to_string(all_cells) + " duplicates: " + std::to_string(duplicates) + "\n";
        message = message + "Front: " + std::to_string(m_size_boundaryFront) + " Back: " + std::to_string(m_size_boundaryBack) + " Bottom: " + std::to_string(m_size_boundaryBottom) + " Top: " + std::to_string(m_size_boundaryTop) + " Left: " + std::to_string(m_size_boundaryLeft) + " Right: " + std::to_string(m_size_boundaryRight) + "\n";
    }
    if ( m_size_boundaryList + m_size_innerList + size_obstacles != nx*ny*nz){
        message = message + "list size of all domain cells is not equal with domain size. Boundary List: " + std::to_string(m_size_boundaryList) + " Inner List: " + std::to_string(m_size_innerList) + " Domain Size: " + std::to_string(domain->get_size(m_level)) + " Obstacle size: " + std::to_string(size_obstacles) + "\n";
    }
    size_t innerCells = (nz - 2) * (ny - 2) * (nx - 2);
    if (m_size_innerList != innerCells - size_obstacles){
        message = message + "list size of inner cell is not equal with domain inner size minus size of obstacles. Inner List: " + std::to_string(m_size_innerList) + " Domain inner size: " + std::to_string(innerCells) + " obstacle size: " + std::to_string(size_obstacles) + "\n";
    }
    size_t startIndex = IX((domain->get_index_x1(m_level) - 1), (domain->get_index_y1(m_level) - 1), (domain->get_index_z1(m_level) - 1), Nx, Ny);
    size_t endIndex = IX((domain->get_index_x2(m_level) + 1), (domain->get_index_y2(m_level) + 1), (domain->get_index_z2(m_level) + 1), Nx, Ny);
    if (*(m_boundaryList) != startIndex || *(m_boundaryList + m_size_boundaryList - 1) != endIndex){
        message = message + "first or last index of boundary list not correct (" + std::to_string(startIndex) + "|" + std::to_string(*(m_boundaryList)) + "(" + std::to_string(endIndex) + "|" + std::to_string(*(m_boundaryList + m_size_boundaryList - 1)) + "\n";
    }
    size_t front2 = IX(domain->get_index_x2(m_level) + 1, domain->get_index_y2(m_level) + 1, domain->get_index_z1(m_level) - 1, Nx, Ny);
    if (*(m_boundaryFront) != startIndex || *(m_boundaryFront + m_size_boundaryFront - 1) != front2){
        message = message + "first or last index of boundary Front not correct (" + std::to_string(startIndex) + "|" + std::to_string(*(m_boundaryFront)) + "(" + std::to_string(front2) + "|" + std::to_string(*(m_boundaryFront + m_size_boundaryFront - 1)) + "\n";
    }
    size_t back1 = IX(domain->get_index_x1(m_level) - 1, domain->get_index_y1(m_level) - 1 , domain->get_index_z2(m_level) + 1, Nx, Ny);
    if (*(m_boundaryBack) != back1 || *(m_boundaryBack + m_size_boundaryBack - 1) != endIndex){
        message = message + "first or last index of boundary Back not correct (" + std::to_string(back1) + "|" + std::to_string(*(m_boundaryBack)) + "(" + std::to_string(endIndex) + "|" + std::to_string(*(m_boundaryBack + m_size_boundaryBack - 1)) + "\n";
    }
    size_t bottom2 = IX(domain->get_index_x2(m_level) + 1, domain->get_index_y1(m_level) - 1, domain->get_index_z2(m_level) + 1, Nx, Ny);
    if (*(m_boundaryBottom) != startIndex || *(m_boundaryBottom + m_size_boundaryBottom - 1) != bottom2){
        message = message + "first or last index of boundary Bottom not correct (" + std::to_string(startIndex) + "|" + std::to_string(*(m_boundaryBottom)) + "(" + std::to_string(bottom2) + "|" + std::to_string(*(m_boundaryBottom + m_size_boundaryBottom - 1)) + "\n";
    }
    size_t top1 = IX(domain->get_index_x1(m_level) - 1, domain->get_index_y2(m_level) + 1, domain->get_index_z1(m_level) - 1, Nx, Ny);
    if (*(m_boundaryTop) !=  top1 || *(m_boundaryTop + m_size_boundaryTop - 1) != endIndex){
        message = message + "first or last index of boundary Back not correct (" + std::to_string(top1) + "|" + std::to_string(*(m_boundaryBack)) + "(" + std::to_string(endIndex) + "|" + std::to_string(*(m_boundaryBack + m_size_boundaryBack - 1)) + "\n";
    }
    size_t left2 = IX(domain->get_index_x1(m_level) - 1, domain->get_index_y2(m_level) + 1, domain->get_index_z2(m_level) + 1, Nx, Ny);
    if (*(m_boundaryLeft) != startIndex || *(m_boundaryLeft + m_size_boundaryLeft - 1) != left2){
        message = message + "first or last index of boundary Left not correct (" + std::to_string(startIndex) + "|" + std::to_string(*(m_boundaryLeft)) + "(" + std::to_string(left2) + "|" + std::to_string(*(m_boundaryLeft + m_size_boundaryLeft - 1)) + "\n";
    }
    size_t right1 = IX(domain->get_index_x2(m_level) + 1, domain->get_index_y1(m_level) - 1, domain->get_index_z1(m_level) - 1, Nx, Ny);
    if (*(m_boundaryRight) != right1 || *(m_boundaryRight + m_size_boundaryRight - 1) != endIndex){
        message = message + "first or last index of boundary Right not correct (" + std::to_string(right1) + "|" + std::to_string(*(m_boundaryRight)) + "(" + std::to_string(endIndex) + "|" + std::to_string(*(m_boundaryRight + m_size_boundaryRight - 1)) + "\n";
    }

    for (size_t i = 1; i < m_size_boundaryList; i++) {
        int diff = static_cast<int>(m_boundaryList[i] - m_boundaryList[i - 1]);
        if (diff < 0) {
            message = message + "sorting error at index " + std::to_string(i - 1) + "|" + std::to_string(i) + " with values " + std::to_string(m_boundaryList[i - 1]) + "|" + std::to_string(m_boundaryList[i]) + "\n";
        }
    }
    if(!message.empty()) {
        message = "############### BOUNDARY CONTROL ###############\n-- level " + std::to_string(m_level) + "\n" + message + "--------------- END BOUNDARY CONTROL ---------------";
#ifndef BENCHMARKING
        auto m_logger = Utility::create_logger(typeid(Boundary).name());
        m_logger->warn(message);
#endif
    }
}

Boundary::~Boundary() {
    delete[] m_innerList;
    delete[] m_boundaryList;
    delete[] m_boundaryFront;
    delete[] m_boundaryBack;
    delete[] m_boundaryTop;
    delete[] m_boundaryBottom;
    delete[] m_boundaryLeft;
    delete[] m_boundaryRight;
}

//======================================== Boundary cells ====================================
// ***************************************************************************************
/// \brief  Creates lists of indices of boundary cells
// ***************************************************************************************
void Boundary::boundaryCells() {
    auto domain = Domain::getInstance();

    const size_t Nx = domain->get_Nx(m_level);
    const size_t Ny = domain->get_Ny(m_level);

    // DETAILED and CONCATENATED LISTS
// BOUNDARY

    //TODO: boundaries for physical domain. New method for computational domain -> redefine XML usage of boundaries

    //start indices for computational domain minus 1 for ghost cells
    size_t x1 = domain->get_index_x1(m_level) - 1;
    size_t y1 = domain->get_index_y1(m_level) - 1;
    size_t z1 = domain->get_index_z1(m_level) - 1;

    //end indices for computational domain plus 1 for ghost cells
    size_t x2 = domain->get_index_x2(m_level) + 1;
    size_t y2 = domain->get_index_y2(m_level) + 1;
    size_t z2 = domain->get_index_z2(m_level) + 1;

    // fill boundaryList with boundary indices of computational domain (sorted)
    size_t counter = 0;
    //go through computational domain in z slices
    //front
    for (size_t j = y1; j <= y2; ++j) {
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, j, z1, Nx, Ny);
            *(m_boundaryList + counter) = idx;
            counter++;
        }
    }
    //left, right, bottom, top
    for (size_t k = z1 + 1; k < z2; ++k) {
        //bottom stride
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, y1, k, Nx, Ny);
            *(m_boundaryList + counter) = idx;
            counter++;
        }
        //cell on the left and on the right
        for (size_t j = y1 + 1; j < y2; ++j) {
            size_t idx = IX(x1, j, k, Nx, Ny);
            *(m_boundaryList + counter) = idx;
            counter++;

            idx = IX(x2, j, k, Nx, Ny);
            *(m_boundaryList + counter) = idx;
            counter++;
        }
        //top stride
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, y2, k, Nx, Ny);
            *(m_boundaryList + counter) = idx;
            counter++;
        }
    }
    //back
    for (size_t j = y1; j <= y2; ++j) {
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, j, z2, Nx, Ny);
            *(m_boundaryList + counter) = idx;
            counter++;
        }
    }

    counter = 0;
    // FRONT and BACK
    for (size_t j = y1; j <= y2; ++j) {
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, j, z1, Nx, Ny);
            *(m_boundaryFront + counter) = idx;

            idx = IX(i, j, z2, Nx, Ny);
            *(m_boundaryBack + counter) = idx;
            counter++;
        }
    }

    counter = 0;
    // TOP and BOTTOM
    for (size_t k = z1; k <= z2; ++k) {
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, y1, k, Nx, Ny);
            *(m_boundaryBottom + counter) = idx;

            idx = IX(i, y2, k, Nx, Ny);
            *(m_boundaryTop + counter) = idx;
            counter++;
        }
    }

    counter = 0;
    // LEFT and RIGHT
    for (size_t k = z1; k <= z2; ++k) {
        for (size_t j = y1; j <= y2; ++j) {
            size_t idx = IX(x1, j, k, Nx, Ny);
            *(m_boundaryLeft + counter) = idx;

            idx = IX(x2, j, k, Nx, Ny);
            *(m_boundaryRight + counter) = idx;
            counter++;
        }
    }
}

//======================================== Inner cells ====================================
// ***************************************************************************************
/// \brief  Creates lists of indices of inner cells
/// \param  obstacleList List of all obstacles of each multigrid level
/// \param  numberOfObstacles Amount of obstacles
// ***************************************************************************************
void Boundary::innerCells(Obstacle **obstacleList, size_t numberOfObstacles) {
    Domain *domain = Domain::getInstance();
    size_t k1 = domain->get_index_z1(m_level);
    size_t j1 = domain->get_index_y1(m_level);
    size_t i1 = domain->get_index_x1(m_level);
    size_t k2 = domain->get_index_z2(m_level);
    size_t j2 = domain->get_index_y2(m_level);
    size_t i2 = domain->get_index_x2(m_level);

    size_t Nx = domain->get_Nx(m_level);
    size_t Ny = domain->get_Ny(m_level);

    size_t counter = 0;
    for (size_t k = k1; k <= k2; ++k) {
        for (size_t j = j1; j <= j2; ++j) {
            for (size_t i = i1; i <= i2; ++i) {
                bool isInnerCell = true;
                //check if cell is part of an obstacle
                for (size_t o = 0; o < numberOfObstacles && isInnerCell; o++) {
                    if (obstacleList[o]->isObstacleCell(i, j, k)) {
                        isInnerCell = false;
                    }
                }
                if (isInnerCell) {
                    size_t idx = IX(i, j, k, Nx, Ny);
                    *(m_innerList + counter) = idx;
                    counter++;
                }
            }
        }
    }
}

//======================================== Inner cells ====================================
// ***************************************************************************************
/// \brief  Creates lists of indices of inner cells without obstacles
// ***************************************************************************************
void Boundary::innerCells() {
    Domain *domain = Domain::getInstance();
    size_t k1 = domain->get_index_z1(m_level);
    size_t j1 = domain->get_index_y1(m_level);
    size_t i1 = domain->get_index_x1(m_level);
    size_t k2 = domain->get_index_z2(m_level);
    size_t j2 = domain->get_index_y2(m_level);
    size_t i2 = domain->get_index_x2(m_level);

    size_t Nx = domain->get_Nx(m_level);
    size_t Ny = domain->get_Ny(m_level);

    size_t counter = 0;
    for (size_t k = k1; k <= k2; ++k) {
        for (size_t j = j1; j <= j2; ++j) {
            for (size_t i = i1; i <= i2; ++i) {
                size_t idx = IX(i, j, k, Nx, Ny);
                *(m_innerList + counter) = idx;
                counter++;
            }
        }
    }
}

//======================================== Update lists ====================================
// ***************************************************************************************
/// \brief  Updates lists of indices
/// \param  obstacleList List of all obstacles of each multigrid level
/// \param  numberOfObstacles Number of obstacles
// ***************************************************************************************
void Boundary::updateLists(Obstacle** obstacleList, size_t numberOfObstacles, size_t size_obstacles) {
    //TODO update for GPU -- delete old data + enter new data
    //TODO write to joined list directly?
    clearLists();
    init(size_obstacles);
    innerCells(obstacleList, numberOfObstacles);
}

//======================================== Update lists ====================================
// ***************************************************************************************
/// \brief  Updates lists of indices
// ***************************************************************************************
void Boundary::updateLists() {
    //TODO update for GPU -- delete old data + enter new data
    //TODO write to joined list directly?
    clearLists();
    init(0);
    innerCells();
}

//======================================== Update lists ====================================
// ***************************************************************************************
/// \brief  removes all allocated arrays
// ***************************************************************************************
void Boundary::clearLists(){
    delete[] m_innerList;
    delete[] m_boundaryList;
    delete[] m_boundaryFront;
    delete[] m_boundaryBack;
    delete[] m_boundaryTop;
    delete[] m_boundaryBottom;
    delete[] m_boundaryLeft;
    delete[] m_boundaryRight;
}

