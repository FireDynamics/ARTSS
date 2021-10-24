/// \file       Boundary.cpp
/// \brief      Data class of boundary object
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Domain.h"
#include <cmath>
#include <string>
#include "../field/Field.h"
#include "../DomainData.h"
#include "../utility/GlobalMacrosTypes.h"

Domain::Domain(size_t multigrid_level) :
        m_multigrid_level(multigrid_level), m_size_boundary() {
    init(0);
    inner_cells();

    Utility::merge_sort(m_inner_list, m_boundary_list, m_size_inner_list, m_size_boundary_list, m_domain_list);
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(Domain).name());
    print(0);
    control(0);
#endif
}

Domain::Domain(
        Obstacle **obstacle_list,
        size_t number_of_obstacles,
        size_t size_obstacles,
        size_t multigrid_level)  :
        m_multigrid_level(multigrid_level), m_size_boundary() {
    init(size_obstacles);
    inner_cells(obstacle_list, number_of_obstacles);

    Utility::merge_sort(m_inner_list, m_boundary_list, m_size_inner_list, m_size_boundary_list, m_domain_list);
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(Domain).name());
    print(size_obstacles);
    control(size_obstacles);
#endif
}


//======================================== Init ====================================================
// *************************************************************************************************
/// \brief  Initialize member variables (arrays)
/// \param  size_obstacles Amount of obstacle cells
// *************************************************************************************************
void Domain::init(size_t size_obstacles) {
    auto domain = DomainData::getInstance();

    const size_t Nx = domain->get_Nx(m_multigrid_level);
    const size_t Ny = domain->get_Ny(m_multigrid_level);
    const size_t Nz = domain->get_Nz(m_multigrid_level);

    m_boundary_patch_divided = new size_t*[number_of_patches];
    m_size_boundary.m_patches[FRONT] = Ny * Nx;
    m_size_boundary[BACK] = Ny * Nx;
    m_boundary_patch_divided[FRONT] = new size_t[m_size_boundary[FRONT]];
    m_boundary_patch_divided[BACK] = new size_t[m_size_boundary[BACK]];

    m_size_boundary[TOP] = Nz * Nx;
    m_size_boundary[BOTTOM] = Nz * Nx;
    m_boundary_patch_divided[BOTTOM] = new size_t[m_size_boundary[BOTTOM]];
    m_boundary_patch_divided[TOP] = new size_t[m_size_boundary[TOP]];

    m_size_boundary[LEFT] = Nz * Ny;
    m_size_boundary[RIGHT] = Nz * Ny;
    m_boundary_patch_divided[LEFT] = new size_t[m_size_boundary[LEFT]];
    m_boundary_patch_divided[RIGHT] = new size_t[m_size_boundary[RIGHT]];

    m_size_inner_list = domain->get_nx(m_multigrid_level) * domain->get_ny(m_multigrid_level) * domain->get_nz(m_multigrid_level) - size_obstacles;
    m_inner_list = new size_t[m_size_inner_list];

    m_size_boundary_list = 2 * Nx * Ny + 2 * (Nz - 2) * (Ny - 2) + 2 * (Nz - 2) * Nx;
    m_boundary_list = new size_t[m_size_boundary_list];

    m_size_domain_list = m_size_inner_list + m_size_boundary_list;
    m_domain_list = new size_t[m_size_domain_list];

    boundary_cells();
}


//======================================== Print ===================================================
// *************************************************************************************************
/// \brief  Prints boundary infos
/// \param  size_obstacles Amount of obstacle cells
// *************************************************************************************************
void Domain::print(size_t size_obstacles) {
#ifndef BENCHMARKING
    m_logger->debug("################ DOMAIN ################");
    m_logger->debug("level: {}", m_multigrid_level);
    m_logger->debug("list size of bList: {}", m_size_domain_list);
    m_logger->debug("Domain starts with {} and ends with {}",
                    *(m_domain_list + 0), *(m_domain_list + m_size_domain_list - 1));
    m_logger->debug("list size of size_z: {}", m_size_boundary[FRONT]);
    m_logger->debug("Front starts with {} and ends with {}",
                    *(m_boundary_patch_divided[FRONT] + 0), *(m_domain_list + m_size_domain_list - 1));
    m_logger->debug("Back starts with {} and ends with {}",
                    *(m_boundary_patch_divided[BACK] + 0), *(m_domain_list + m_size_domain_list - 1));
    m_logger->debug("list size of size_y: ", m_size_boundary[BOTTOM]);
    m_logger->debug("Bottom starts with {} and ends with {}",
                    *(m_boundary_patch_divided[BOTTOM] + 0), *(m_boundary_patch_divided[BOTTOM] + m_size_boundary[BOTTOM] - 1));
    m_logger->debug("Top starts with {} and ends with {}",
                    *(m_boundary_patch_divided[TOP] + 0), *(m_boundary_patch_divided[TOP] + m_size_boundary[TOP] - 1));
    m_logger->debug("list size of size_x: ", m_size_boundary[LEFT]);
    m_logger->debug("Left starts with {} and ends with {}",
                    *(m_boundary_patch_divided[LEFT] + 0), *(m_boundary_patch_divided[LEFT] + m_size_boundary[LEFT] - 1));
    m_logger->debug("Right starts with {} and ends with {}",
                    *(m_boundary_patch_divided[RIGHT] + 0), *(m_boundary_patch_divided[RIGHT] + m_size_boundary[RIGHT] - 1));
    m_logger->debug("list size of innerList: {} obstacle size: {}", m_size_inner_list, size_obstacles);
    m_logger->debug("Inner starts with {} and ends with {}",
                    *(m_inner_list + 0), *(m_inner_list + m_size_inner_list - 1));
    m_logger->debug("--------------- END DOMAIN ---------------");
#endif
}

//======================================== Control =================================================
// *************************************************************************************************
/// \brief  Units test emergency solution
/// \param  size_obstacles Amount of obstacle cells
// *************************************************************************************************
void Domain::control(size_t size_obstacles) {
    // TODO(n16h7): clean up
    std::string message;
    DomainData *domain = DomainData::getInstance();
    size_t size = domain->get_size(m_multigrid_level);

    size_t nx = domain->get_nx(m_multigrid_level);
    size_t ny = domain->get_ny(m_multigrid_level);
    size_t nz = domain->get_nz(m_multigrid_level);
    size_t Nx = domain->get_Nx(m_multigrid_level);
    size_t Ny = domain->get_Ny(m_multigrid_level);
    size_t all_cells = m_size_inner_list
                     + m_size_boundary[FRONT] + m_size_boundary[BACK]
                     + m_size_boundary[BOTTOM] + m_size_boundary[TOP]
                     + m_size_boundary[LEFT] + m_size_boundary[RIGHT];
    size_t duplicates = 4 * Nx + 4 * Ny + 4 * nz;
    if (m_size_domain_list != all_cells - duplicates) {
        message = message + "list size of all domain cells does not fit with sum of it parts. Domain List: "
                  + std::to_string(m_size_domain_list) + " sum: " + std::to_string(all_cells)
                  + " duplicates: " + std::to_string(duplicates) + "\n";
        message = message + "Front: " + std::to_string(m_size_boundary[FRONT])
                          + " Back: " + std::to_string(m_size_boundary[BACK])
                          + " Bottom: " + std::to_string(m_size_boundary[BOTTOM])
                          + " Top: " + std::to_string(m_size_boundary[TOP])
                          + " Left: " + std::to_string(m_size_boundary[LEFT])
                          + " Right: " + std::to_string(m_size_boundary[RIGHT]) + "\n";
    }
    if (m_size_domain_list + size_obstacles != size) {
        message = message + "list size of all domain cells is not equal with domain size."
                            " Domain List: " + std::to_string(m_size_domain_list)
                            + " Domain Size: " + std::to_string(domain->get_size(m_multigrid_level))
                            + " Obstacle size: " + std::to_string(size_obstacles) + "\n";
    }
    size_t innerCells = nz * ny * nx;
    if (m_size_inner_list != innerCells - size_obstacles) {
        message = message + "list size of inner cell is not equal with domain inner size minus size of obstacles."
                            " Inner List: " + std::to_string(m_size_inner_list)
                            + " Domain inner size: " + std::to_string(innerCells)
                            + " Obstacle size: " + std::to_string(size_obstacles) + "\n";
    }
    size_t startIndex = IX(
            (domain->get_index_x1(m_multigrid_level) - 1),
            (domain->get_index_y1(m_multigrid_level) - 1),
            (domain->get_index_z1(m_multigrid_level) - 1),
            Nx, Ny);
    size_t endIndex = IX(
            (domain->get_index_x2(m_multigrid_level) + 1),
            (domain->get_index_y2(m_multigrid_level) + 1),
            (domain->get_index_z2(m_multigrid_level) + 1),
            Nx, Ny);
    if (*(m_domain_list) != startIndex || *(m_domain_list + m_size_domain_list - 1) != endIndex) {
        message = message + "first or last index of domain list not correct ("
                  + std::to_string(startIndex) + "|"
                  + std::to_string(*(m_domain_list)) + ")("
                  + std::to_string(endIndex) + "|"
                  + std::to_string(*(m_domain_list + m_size_domain_list - 1)) + ")\n";
    }
    size_t front2 = IX(
            domain->get_index_x2(m_multigrid_level) + 1,
            domain->get_index_y2(m_multigrid_level) + 1,
            domain->get_index_z1(m_multigrid_level) - 1,
            Nx, Ny);
    if (m_boundary_patch_divided[FRONT][0] != startIndex || m_boundary_patch_divided[FRONT][m_size_boundary[FRONT] - 1] != front2) {
        message = message + "first or last index of boundary Front not correct ("
                  + std::to_string(startIndex) + "|"
                  + std::to_string(*(m_boundary_patch_divided[FRONT])) + ")("
                  + std::to_string(front2) + "|"
                  + std::to_string(*(m_boundary_patch_divided[FRONT] + m_size_boundary[FRONT] - 1)) + ")\n";
    }
    size_t back1 = IX(
            domain->get_index_x1(m_multigrid_level) - 1,
            domain->get_index_y1(m_multigrid_level) - 1,
            domain->get_index_z2(m_multigrid_level) + 1,
            Nx, Ny);
    if (*(m_boundary_patch_divided[BACK]) != back1 || *(m_boundary_patch_divided[BACK] + m_size_boundary[BACK] - 1) != endIndex) {
        message = message + "first or last index of boundary Back not correct ("
                  + std::to_string(back1) + "|"
                  + std::to_string(*(m_boundary_patch_divided[BACK])) + ")("
                  + std::to_string(endIndex) + "|"
                  + std::to_string(*(m_boundary_patch_divided[BACK] + m_size_boundary[BACK] - 1)) + ")\n";
    }
    size_t bottom2 = IX(
            domain->get_index_x2(m_multigrid_level) + 1,
            domain->get_index_y1(m_multigrid_level) - 1,
            domain->get_index_z2(m_multigrid_level) + 1,
            Nx, Ny);
    if (*(m_boundary_patch_divided[BOTTOM]) != startIndex || *(m_boundary_patch_divided[BOTTOM] + m_size_boundary[BOTTOM] - 1) != bottom2) {
        message = message + "first or last index of boundary Bottom not correct ("
                  + std::to_string(startIndex) + "|"
                  + std::to_string(*(m_boundary_patch_divided[BOTTOM])) + ")("
                  + std::to_string(bottom2) + "|"
                  + std::to_string(*(m_boundary_patch_divided[BOTTOM] + m_size_boundary[BOTTOM] - 1)) + ")\n";
    }
    size_t top1 = IX(
            domain->get_index_x1(m_multigrid_level) - 1,
            domain->get_index_y2(m_multigrid_level) + 1,
            domain->get_index_z1(m_multigrid_level) - 1,
            Nx, Ny);
    if (*(m_boundary_patch_divided[TOP]) != top1 || *(m_boundary_patch_divided[TOP] + m_size_boundary[TOP] - 1) != endIndex) {
        message = message + "first or last index of boundary Top not correct ("
                  + std::to_string(top1) + "|"
                  + std::to_string(*(m_boundary_patch_divided[TOP])) + ")("
                  + std::to_string(endIndex) + "|"
                  + std::to_string(*(m_boundary_patch_divided[TOP] + m_size_boundary[TOP] - 1)) + ")\n";
    }
    size_t left2 = IX(
            domain->get_index_x1(m_multigrid_level) - 1,
            domain->get_index_y2(m_multigrid_level) + 1,
            domain->get_index_z2(m_multigrid_level) + 1,
            Nx, Ny);
    if (*(m_boundary_patch_divided[LEFT]) != startIndex || *(m_boundary_patch_divided[LEFT] + m_size_boundary[LEFT] - 1) != left2) {
        message = message + "first or last index of boundary Left not correct ("
                  + std::to_string(startIndex) + "|"
                  + std::to_string(*(m_boundary_patch_divided[LEFT])) + ")("
                  + std::to_string(left2) + "|"
                  + std::to_string(*(m_boundary_patch_divided[LEFT] + m_size_boundary[LEFT] - 1)) + ")\n";
    }
    size_t right1 = IX(
            domain->get_index_x2(m_multigrid_level) + 1,
            domain->get_index_y1(m_multigrid_level) - 1,
            domain->get_index_z1(m_multigrid_level) - 1,
            Nx, Ny);
    if (*(m_boundary_patch_divided[RIGHT]) != right1 || *(m_boundary_patch_divided[RIGHT] + m_size_boundary[RIGHT] - 1) != endIndex) {
        message = message + "first or last index of boundary Right not correct ("
                  + std::to_string(right1) + "|"
                  + std::to_string(*(m_boundary_patch_divided[RIGHT])) + ")("
                  + std::to_string(endIndex) + "|"
                  + std::to_string(*(m_boundary_patch_divided[RIGHT] + m_size_boundary[RIGHT] - 1)) + ")\n";
    }

    for (size_t i = 1; i < m_size_domain_list; i++) {
        int diff = static_cast<int>(m_domain_list[i] - m_domain_list[i - 1]);
        if (diff < 0) {
            message = message + "sorting error at index "
                      + std::to_string(i - 1) + "|"
                      + std::to_string(i) + " with values "
                      + std::to_string(m_domain_list[i - 1]) + "|"
                      + std::to_string(m_domain_list[i]) + "\n";
        }
    }
    if (!message.empty()) {
        message = "############### DOMAIN CONTROL ###############\n-- level "
                  + std::to_string(m_multigrid_level) + "\n" + message
                  + "--------------- END DOMAIN CONTROL ---------------";
#ifndef BENCHMARKING
        m_logger->warn(message);
#endif
    }
}

Domain::~Domain() {
    delete[] m_inner_list;
    delete[] m_domain_list;
    delete[] m_boundary_patch_divided[FRONT];
    delete[] m_boundary_patch_divided[BACK];
    delete[] m_boundary_patch_divided[TOP];
    delete[] m_boundary_patch_divided[BOTTOM];
    delete[] m_boundary_patch_divided[LEFT];
    delete[] m_boundary_patch_divided[RIGHT];
}

//======================================== Boundary cells ==========================================
// *************************************************************************************************
/// \brief  Creates lists of indices of boundary cells
// *************************************************************************************************
void Domain::boundary_cells() {
    auto domain = DomainData::getInstance();

    const size_t Nx = domain->get_Nx(m_multigrid_level);
    const size_t Ny = domain->get_Ny(m_multigrid_level);

    // DETAILED and CONCATENATED LISTS
    // BOUNDARY

    // TODO(issue 86): boundaries for physical domain.
    // TODO(issue 86): New method for computational domain -> redefine XML usage of boundaries

    // start indices for computational domain minus 1 for ghost cells
    size_t x1 = domain->get_index_x1(m_multigrid_level) - 1;
    size_t y1 = domain->get_index_y1(m_multigrid_level) - 1;
    size_t z1 = domain->get_index_z1(m_multigrid_level) - 1;

    // end indices for computational domain plus 1 for ghost cells
    size_t x2 = domain->get_index_x2(m_multigrid_level) + 1;
    size_t y2 = domain->get_index_y2(m_multigrid_level) + 1;
    size_t z2 = domain->get_index_z2(m_multigrid_level) + 1;

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
            *(m_boundary_patch_divided[FRONT] + counter) = idx;

            idx = IX(i, j, z2, Nx, Ny);
            *(m_boundary_patch_divided[BACK] + counter) = idx;
            counter++;
        }
    }

    counter = 0;
    // TOP and BOTTOM
    for (size_t k = z1; k <= z2; ++k) {
        for (size_t i = x1; i <= x2; ++i) {
            size_t idx = IX(i, y1, k, Nx, Ny);
            *(m_boundary_patch_divided[BOTTOM] + counter) = idx;

            idx = IX(i, y2, k, Nx, Ny);
            *(m_boundary_patch_divided[TOP] + counter) = idx;
            counter++;
        }
    }

    counter = 0;
    // LEFT and RIGHT
    for (size_t k = z1; k <= z2; ++k) {
        for (size_t j = y1; j <= y2; ++j) {
            size_t idx = IX(x1, j, k, Nx, Ny);
            *(m_boundary_patch_divided[LEFT] + counter) = idx;

            idx = IX(x2, j, k, Nx, Ny);
            *(m_boundary_patch_divided[RIGHT] + counter) = idx;
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
void Domain::inner_cells(Obstacle **obstacle_list, size_t number_of_obstacles) {
    DomainData *domain = DomainData::getInstance();
    size_t k1 = domain->get_index_z1(m_multigrid_level);
    size_t j1 = domain->get_index_y1(m_multigrid_level);
    size_t i1 = domain->get_index_x1(m_multigrid_level);
    size_t k2 = domain->get_index_z2(m_multigrid_level);
    size_t j2 = domain->get_index_y2(m_multigrid_level);
    size_t i2 = domain->get_index_x2(m_multigrid_level);

    size_t Nx = domain->get_Nx(m_multigrid_level);
    size_t Ny = domain->get_Ny(m_multigrid_level);

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
void Domain::inner_cells() {
    DomainData *domain = DomainData::getInstance();
    size_t k1 = domain->get_index_z1(m_multigrid_level);
    size_t j1 = domain->get_index_y1(m_multigrid_level);
    size_t i1 = domain->get_index_x1(m_multigrid_level);
    size_t k2 = domain->get_index_z2(m_multigrid_level);
    size_t j2 = domain->get_index_y2(m_multigrid_level);
    size_t i2 = domain->get_index_x2(m_multigrid_level);

    size_t Nx = domain->get_Nx(m_multigrid_level);
    size_t Ny = domain->get_Ny(m_multigrid_level);

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
void Domain::update_lists(Obstacle **obstacle_list, size_t number_of_obstacles, size_t size_obstacles) {
    clear_lists();
    init(size_obstacles);
    inner_cells(obstacle_list, number_of_obstacles);
}

//======================================== Update lists ============================================
// *************************************************************************************************
/// \brief  Updates lists of indices
// *************************************************************************************************
void Domain::update_lists() {
    clear_lists();
    init(0);
    inner_cells();
}

//======================================== Update lists ============================================
// *************************************************************************************************
/// \brief  removes all allocated arrays
// *************************************************************************************************
void Domain::clear_lists() {
    delete[] m_inner_list;
    delete[] m_domain_list;
    delete[] m_boundary_patch_divided[FRONT];
    delete[] m_boundary_patch_divided[BACK];
    delete[] m_boundary_patch_divided[BOTTOM];
    delete[] m_boundary_patch_divided[TOP];
    delete[] m_boundary_patch_divided[LEFT];
    delete[] m_boundary_patch_divided[RIGHT];
}
