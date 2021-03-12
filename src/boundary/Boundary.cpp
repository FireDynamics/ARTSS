/// \file       Boundary.cpp
/// \brief      Data class of boundary object
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Boundary.h"
#include <cmath>
#include <string>
#include "../field/Field.h"
#include "../Domain.h"
#include "../utility/GlobalMacrosTypes.h"

Boundary::Boundary(size_t level) {
    m_level = level;
    init(0);
    inner_cells();

#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(Boundary).name());
    print(0);
    control(0);
#endif
}

Boundary::Boundary(Obstacle **obstacle_list, size_t number_of_obstacles, size_t size_obstacles, size_t level) {
    m_level = level;
    init(size_obstacles);
    inner_cells(obstacle_list, number_of_obstacles);

#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(Boundary).name());
    print(size_obstacles);
    control(size_obstacles);
#endif
}


//======================================== Init ====================================================
// *************************************************************************************************
/// \brief  Initialize member variables (arrays)
/// \param  size_obstacles Amount of obstacle cells
// *************************************************************************************************
void Boundary::init(size_t size_obstacles) {
    auto domain = Domain::getInstance();

    const size_t nx = domain->get_nx(m_level);
    const size_t ny = domain->get_ny(m_level);
    const size_t nz = domain->get_nz(m_level);

    m_size_boundary_list = 2 * nx * ny + 2 * (nz - 2) * (ny - 2) + 2 * (nz - 2) * nx;
    m_boundary_list = new size_t[m_size_boundary_list];

    m_size_boundary_front = ny * nx;
    m_size_boundary_back = ny * nx;
    m_boundary_front = new size_t[m_size_boundary_front];
    m_boundary_back = new size_t[m_size_boundary_back];

    m_size_boundary_top = nz * nx;
    m_size_boundary_bottom = nz * nx;
    m_boundary_bottom = new size_t[m_size_boundary_bottom];
    m_boundary_top = new size_t[m_size_boundary_top];

    m_size_boundary_left = nz * ny;
    m_size_boundary_right = nz * ny;
    m_boundary_left = new size_t[m_size_boundary_left];
    m_boundary_right = new size_t[m_size_boundary_right];

    m_size_innerList = (nx - 2) * (ny - 2) * (nz - 2) - size_obstacles;
    m_inner_list = new size_t[m_size_innerList];

    boundary_cells();
}


//======================================== Print ===================================================
// *************************************************************************************************
/// \brief  Prints boundary infos
/// \param  size_obstacles Amount of obstacle cells
// *************************************************************************************************
void Boundary::print(size_t size_obstacles) {
#ifndef BENCHMARKING
    m_logger->debug("################ BOUNDARY ################");
    m_logger->debug("list size of bList: {}", m_size_boundary_list);
    m_logger->debug("Boundary starts with {} and ends with {}",
                    *(m_boundary_list + 0), *(m_boundary_list + m_size_boundary_list - 1));
    m_logger->debug("list size of size_z: {}", m_size_boundary_front);
    m_logger->debug("Front starts with {} and ends with {}",
                    *(m_boundary_front + 0), *(m_boundary_list + m_size_boundary_list - 1));
    m_logger->debug("Back starts with {} and ends with {}",
                    *(m_boundary_back + 0), *(m_boundary_list + m_size_boundary_list - 1));
    m_logger->debug("list size of size_y: ", m_size_boundary_bottom);
    m_logger->debug("Bottom starts with {} and ends with {}",
                    *(m_boundary_bottom + 0), *(m_boundary_bottom + m_size_boundary_bottom - 1));
    m_logger->debug("Top starts with {} and ends with {}",
                    *(m_boundary_top + 0), *(m_boundary_top + m_size_boundary_top - 1));
    m_logger->debug("list size of size_x: ", m_size_boundary_left);
    m_logger->debug("Left starts with {} and ends with {}",
                    *(m_boundary_left + 0), *(m_boundary_left + m_size_boundary_left - 1));
    m_logger->debug("Right starts with {} and ends with {}",
                    *(m_boundary_right + 0), *(m_boundary_right + m_size_boundary_right - 1));
    m_logger->debug("list size of innerList: {} obstacle size: {}", m_size_innerList, size_obstacles);
    m_logger->debug("Inner starts with {} and ends with {}",
                    *(m_inner_list + 0), *(m_inner_list + m_size_innerList - 1));
    m_logger->debug("--------------- END BOUNDARY ---------------");
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
    Domain *domain = Domain::getInstance();
    size_t nx = domain->get_nx(m_level);
    size_t ny = domain->get_ny(m_level);
    size_t nz = domain->get_nz(m_level);
    size_t Nx = domain->get_Nx(m_level);
    size_t Ny = domain->get_Ny(m_level);
    size_t all_cells = m_size_boundary_front + m_size_boundary_back
                       + m_size_boundary_bottom + m_size_boundary_top
                       + m_size_boundary_left + m_size_boundary_right;
    size_t duplicates = 4 * nx + 4 * (ny - 2) + 4 * (nz - 2) + 8;
    if (m_size_boundary_list != all_cells - duplicates) {
        message = message + "list size of all boundary cells does not fit with sum of it parts. Boundary List: "
                  + std::to_string(m_size_boundary_list) + " sum: " + std::to_string(all_cells)
                  + " duplicates: " + std::to_string(duplicates) + "\n";
        message = message + "Front: " + std::to_string(m_size_boundary_front)
                          + " Back: " + std::to_string(m_size_boundary_back)
                          + " Bottom: " + std::to_string(m_size_boundary_bottom)
                          + " Top: " + std::to_string(m_size_boundary_top)
                          + " Left: " + std::to_string(m_size_boundary_left)
                          + " Right: " + std::to_string(m_size_boundary_right) + "\n";
    }
    if (m_size_boundary_list + m_size_innerList + size_obstacles != nx * ny * nz) {
        message = message + "list size of all domain cells is not equal with domain size."
                            "Boundary List: " + std::to_string(m_size_boundary_list)
                            + " Inner List: " + std::to_string(m_size_innerList)
                            + " Domain Size: " + std::to_string(domain->get_size(m_level))
                            + " Obstacle size: " + std::to_string(size_obstacles) + "\n";
    }
    size_t innerCells = (nz - 2) * (ny - 2) * (nx - 2);
    if (m_size_innerList != innerCells - size_obstacles) {
        message = message + "list size of inner cell is not equal with domain inner size minus size of obstacles."
                            " Inner List: " + std::to_string(m_size_innerList)
                            + " Domain inner size: " + std::to_string(innerCells)
                            + " Obstacle size: " + std::to_string(size_obstacles) + "\n";
    }
    size_t startIndex = IX((domain->get_index_x1(m_level) - 1), (domain->get_index_y1(m_level) - 1), (domain->get_index_z1(m_level) - 1), Nx, Ny);
    size_t endIndex = IX((domain->get_index_x2(m_level) + 1), (domain->get_index_y2(m_level) + 1), (domain->get_index_z2(m_level) + 1), Nx, Ny);
    if (*(m_boundary_list) != startIndex || *(m_boundary_list + m_size_boundary_list - 1) != endIndex) {
        message = message + "first or last index of boundary list not correct ("
                  + std::to_string(startIndex) + "|" + std::to_string(*(m_boundary_list)) + ")("
                  + std::to_string(endIndex) + "|" + std::to_string(*(m_boundary_list + m_size_boundary_list - 1)) + ")\n";
    }
    size_t front2 = IX(domain->get_index_x2(m_level) + 1, domain->get_index_y2(m_level) + 1, domain->get_index_z1(m_level) - 1, Nx, Ny);
    if (*(m_boundary_front) != startIndex || *(m_boundary_front + m_size_boundary_front - 1) != front2) {
        message = message + "first or last index of boundary Front not correct ("
                  + std::to_string(startIndex) + "|" + std::to_string(*(m_boundary_front)) + ")("
                  + std::to_string(front2) + "|" + std::to_string(*(m_boundary_front + m_size_boundary_front - 1)) + ")\n";
    }
    size_t back1 = IX(domain->get_index_x1(m_level) - 1, domain->get_index_y1(m_level) - 1, domain->get_index_z2(m_level) + 1, Nx, Ny);
    if (*(m_boundary_back) != back1 || *(m_boundary_back + m_size_boundary_back - 1) != endIndex) {
        message = message + "first or last index of boundary Back not correct ("
                  + std::to_string(back1) + "|" + std::to_string(*(m_boundary_back)) + ")("
                  + std::to_string(endIndex) + "|" + std::to_string(*(m_boundary_back + m_size_boundary_back - 1)) + ")\n";
    }
    size_t bottom2 = IX(domain->get_index_x2(m_level) + 1, domain->get_index_y1(m_level) - 1, domain->get_index_z2(m_level) + 1, Nx, Ny);
    if (*(m_boundary_bottom) != startIndex || *(m_boundary_bottom + m_size_boundary_bottom - 1) != bottom2) {
        message = message + "first or last index of boundary Bottom not correct ("
                  + std::to_string(startIndex) + "|" + std::to_string(*(m_boundary_bottom)) + ")("
                  + std::to_string(bottom2) + "|" + std::to_string(*(m_boundary_bottom + m_size_boundary_bottom - 1)) + ")\n";
    }
    size_t top1 = IX(domain->get_index_x1(m_level) - 1, domain->get_index_y2(m_level) + 1, domain->get_index_z1(m_level) - 1, Nx, Ny);
    if (*(m_boundary_top) != top1 || *(m_boundary_top + m_size_boundary_top - 1) != endIndex) {
        message = message + "first or last index of boundary Top not correct ("
                  + std::to_string(top1) + "|" + std::to_string(*(m_boundary_top)) + ")("
                  + std::to_string(endIndex) + "|" + std::to_string(*(m_boundary_top + m_size_boundary_top - 1)) + ")\n";
    }
    size_t left2 = IX(domain->get_index_x1(m_level) - 1, domain->get_index_y2(m_level) + 1, domain->get_index_z2(m_level) + 1, Nx, Ny);
    if (*(m_boundary_left) != startIndex || *(m_boundary_left + m_size_boundary_left - 1) != left2) {
        message = message + "first or last index of boundary Left not correct ("
                  + std::to_string(startIndex) + "|" + std::to_string(*(m_boundary_left)) + ")("
                  + std::to_string(left2) + "|" + std::to_string(*(m_boundary_left + m_size_boundary_left - 1)) + ")\n";
    }
    size_t right1 = IX(domain->get_index_x2(m_level) + 1, domain->get_index_y1(m_level) - 1, domain->get_index_z1(m_level) - 1, Nx, Ny);
    if (*(m_boundary_right) != right1 || *(m_boundary_right + m_size_boundary_right - 1) != endIndex) {
        message = message + "first or last index of boundary Right not correct ("
                  + std::to_string(right1) + "|" + std::to_string(*(m_boundary_right)) + ")("
                  + std::to_string(endIndex) + "|" + std::to_string(*(m_boundary_right + m_size_boundary_right - 1)) + ")\n";
    }

    for (size_t i = 1; i < m_size_boundary_list; i++) {
        int diff = static_cast<int>(m_boundary_list[i] - m_boundary_list[i - 1]);
        if (diff < 0) {
            message = message + "sorting error at index " + std::to_string(i - 1) + "|" + std::to_string(i)
                      + " with values " + std::to_string(m_boundary_list[i - 1]) + "|" + std::to_string(m_boundary_list[i]) + "\n";
        }
    }
    if (!message.empty()) {
        message = "############### BOUNDARY CONTROL ###############\n-- level "
                  + std::to_string(m_level) + "\n" + message
                  + "--------------- END BOUNDARY CONTROL ---------------";
#ifndef BENCHMARKING
        m_logger->warn(message);
#endif
    }
}

Boundary::~Boundary() {
    delete[] m_inner_list;
    delete[] m_boundary_list;
    delete[] m_boundary_front;
    delete[] m_boundary_back;
    delete[] m_boundary_top;
    delete[] m_boundary_bottom;
    delete[] m_boundary_left;
    delete[] m_boundary_right;
}

//======================================== Boundary cells ==========================================
// *************************************************************************************************
/// \brief  Creates lists of indices of boundary cells
// *************************************************************************************************
void Boundary::boundary_cells() {
    auto domain = Domain::getInstance();

    const size_t Nx = domain->get_Nx(m_level);
    const size_t Ny = domain->get_Ny(m_level);

    // DETAILED and CONCATENATED LISTS
    // BOUNDARY

    // TODO(issue 86): boundaries for physical domain.
    // TODO(issue 86): New method for computational domain -> redefine XML usage of boundaries

    // start indices for computational domain minus 1 for ghost cells
    size_t x1 = domain->get_index_x1(m_level) - 1;
    size_t y1 = domain->get_index_y1(m_level) - 1;
    size_t z1 = domain->get_index_z1(m_level) - 1;

    // end indices for computational domain plus 1 for ghost cells
    size_t x2 = domain->get_index_x2(m_level) + 1;
    size_t y2 = domain->get_index_y2(m_level) + 1;
    size_t z2 = domain->get_index_z2(m_level) + 1;

    // fill boundaryList with boundary indices of computational domain (sorted)
    size_t counter = 0;
    // go through computational domain in z slices
    // front
    for (size_t j = y1; j <= y2; ++j) {
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, j, z1, Nx, Ny);
            *(m_boundary_list + counter) = idx;
            counter++;
        }
    }
    // left, right, bottom, top
    for (size_t k = z1 + 1; k < z2; ++k) {
        // bottom stride
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, y1, k, Nx, Ny);
            *(m_boundary_list + counter) = idx;
            counter++;
        }
        // cell on the left and on the right
        for (size_t j = y1 + 1; j < y2; ++j) {
            size_t idx = IX(x1, j, k, Nx, Ny);
            *(m_boundary_list + counter) = idx;
            counter++;

            idx = IX(x2, j, k, Nx, Ny);
            *(m_boundary_list + counter) = idx;
            counter++;
        }
        // top stride
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, y2, k, Nx, Ny);
            *(m_boundary_list + counter) = idx;
            counter++;
        }
    }
    // back
    for (size_t j = y1; j <= y2; ++j) {
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, j, z2, Nx, Ny);
            *(m_boundary_list + counter) = idx;
            counter++;
        }
    }

    counter = 0;
    // FRONT and BACK
    for (size_t j = y1; j <= y2; ++j) {
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, j, z1, Nx, Ny);
            *(m_boundary_front + counter) = idx;

            idx = IX(i, j, z2, Nx, Ny);
            *(m_boundary_back + counter) = idx;
            counter++;
        }
    }

    counter = 0;
    // TOP and BOTTOM
    for (size_t k = z1; k <= z2; ++k) {
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, y1, k, Nx, Ny);
            *(m_boundary_bottom + counter) = idx;

            idx = IX(i, y2, k, Nx, Ny);
            *(m_boundary_top + counter) = idx;
            counter++;
        }
    }

    counter = 0;
    // LEFT and RIGHT
    for (size_t k = z1; k <= z2; ++k) {
        for (size_t j = y1; j <= y2; ++j) {
            size_t idx = IX(x1, j, k, Nx, Ny);
            *(m_boundary_left + counter) = idx;

            idx = IX(x2, j, k, Nx, Ny);
            *(m_boundary_right + counter) = idx;
            counter++;
        }
    }
}

//======================================== Inner cells =============================================
// *************************************************************************************************
/// \brief  Creates lists of indices of inner cells
/// \param  obstacle_list List of all obstacles of each multigrid level
/// \param  number_of_obstacles Amount of obstacles
// *************************************************************************************************
void Boundary::inner_cells(Obstacle **obstacle_list, size_t number_of_obstacles) {
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
    size_t counter_obstacle = 0;
    for (size_t k = k1; k <= k2; ++k) {
        for (size_t j = j1; j <= j2; ++j) {
            for (size_t i = i1; i <= i2; ++i) {
                bool is_inner_cell = true;
                // check if cell is part of an obstacle
                for (size_t o = 0; o < number_of_obstacles && is_inner_cell; o++) {
                    if (obstacle_list[o]->is_obstacle_cell(i, j, k)) {
                        is_inner_cell = false;
                        counter_obstacle++;
                    }
                }
                if (is_inner_cell) {
                    size_t idx = IX(i, j, k, Nx, Ny);
                    *(m_inner_list + counter) = idx;
                    counter++;
                }
            }
        }
    }
}

//======================================== Inner cells =============================================
// *************************************************************************************************
/// \brief  Creates lists of indices of inner cells without obstacles
// *************************************************************************************************
void Boundary::inner_cells() {
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
                *(m_inner_list + counter) = idx;
                counter++;
            }
        }
    }
}

//======================================== Update lists ============================================
// *************************************************************************************************
/// \brief  Updates lists of indices
/// \param  obstacle_list List of all obstacles of each multigrid level
/// \param  number_of_obstacles Number of obstacles
// *************************************************************************************************
void Boundary::update_lists(Obstacle **obstacle_list, size_t number_of_obstacles, size_t size_obstacles) {
    clear_lists();
    init(size_obstacles);
    inner_cells(obstacle_list, number_of_obstacles);
}

//======================================== Update lists ============================================
// *************************************************************************************************
/// \brief  Updates lists of indices
// *************************************************************************************************
void Boundary::update_lists() {
    clear_lists();
    init(0);
    inner_cells();
}

//======================================== Update lists ============================================
// *************************************************************************************************
/// \brief  removes all allocated arrays
// *************************************************************************************************
void Boundary::clear_lists() {
    delete[] m_inner_list;
    delete[] m_boundary_list;
    delete[] m_boundary_front;
    delete[] m_boundary_back;
    delete[] m_boundary_top;
    delete[] m_boundary_bottom;
    delete[] m_boundary_left;
    delete[] m_boundary_right;
}
