/// \file       Boundary.cpp
/// \brief      Data class of boundary object
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Domain.h"
#include <cmath>
#include <string>
#include "../DomainData.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Algorithm.h"

Domain::Domain(
        size_t *obstacle_list,
        size_t size_obstacle_list,
        size_t **surface_list,
        PatchObject &size_surface_list,
        size_t multigrid_level) :
        m_multigrid_level(multigrid_level), m_size_boundary() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(Domain).name());
#endif
    init(size_obstacle_list, size_surface_list);
    inner_cells(obstacle_list, size_obstacle_list);
    boundary_cells(surface_list, size_surface_list);

    Algorithm::merge_sort(m_inner_list, m_boundary_list, m_size_inner_list, m_size_boundary_list, m_domain_list);
#ifndef BENCHMARKING
    print(size_obstacle_list, size_surface_list);
    control(size_obstacle_list, size_surface_list);
#endif
}


//======================================== Init ====================================================
// *************************************************************************************************
/// \brief  Initialize member variables (arrays)
/// \param  size_obstacles Amount of obstacle cells
// *************************************************************************************************
void Domain::init(size_t size_obstacle_list, PatchObject &size_surface_list) {
    auto domain_data = DomainData::getInstance();

    const size_t Nx = domain_data->get_Nx(m_multigrid_level);
    const size_t Ny = domain_data->get_Ny(m_multigrid_level);
    const size_t Nz = domain_data->get_Nz(m_multigrid_level);

    m_boundary_patch_divided = new size_t *[number_of_patches];
    m_size_boundary[FRONT] = Ny * Nx - size_surface_list[FRONT];
    m_size_boundary[BACK] = Ny * Nx - size_surface_list[BACK];
    m_size_boundary[TOP] = Nz * Nx - size_surface_list[TOP];
    m_size_boundary[BOTTOM] = Nz * Nx - size_surface_list[BOTTOM];
    m_size_boundary[LEFT] = Nz * Ny - size_surface_list[LEFT];
    m_size_boundary[RIGHT] = Nz * Ny - size_surface_list[RIGHT];

    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_boundary_patch_divided[patch] = new size_t[m_size_boundary[patch]];
    }

    m_size_inner_list = domain_data->get_nx(m_multigrid_level) * domain_data->get_ny(m_multigrid_level) * domain_data->get_nz(m_multigrid_level) - size_obstacle_list;
    m_inner_list = new size_t[m_size_inner_list];

    // surfaces cannot be at the edge/corner cells
    m_size_boundary_list = 2 * Nx * Ny + 2 * (Nz - 2) * (Ny - 2) + 2 * (Nz - 2) * Nx - size_surface_list.get_sum();
    m_boundary_list = new size_t[m_size_boundary_list];

    m_size_domain_list = m_size_inner_list + m_size_boundary_list;
    m_domain_list = new size_t[m_size_domain_list];
}


//======================================== Print ===================================================
// *************************************************************************************************
/// \brief  Prints boundary infos
/// \param  size_obstacle_list Amount of obstacle cells
/// \param  size_surface_list amount of surface cells divided by patch
// *************************************************************************************************
void Domain::print(size_t size_obstacle_list, PatchObject &size_surface_list) {
#ifndef BENCHMARKING
    m_logger->debug("################ DOMAIN ################");
    m_logger->debug("level: {}", m_multigrid_level);
    m_logger->debug("list size of domain list: {}", m_size_domain_list);
    m_logger->debug("Domain starts with {} and ends with {}",
                    m_domain_list[0], m_domain_list[m_size_domain_list - 1]);
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        std::string patch_name = PatchObject::get_patch_name(patch);
        m_logger->debug("list size of domain boundary {}: {} surface size: {}",
                        patch_name,
                        m_size_boundary[patch],
                        size_surface_list[patch]);
        m_logger->debug("{} starts with {} and ends with {}", patch_name,
                        m_boundary_patch_divided[FRONT][0],
                        m_boundary_patch_divided[FRONT][m_size_boundary[patch] - 1]);
    }
    m_logger->debug("list size of inner list: {} obstacle size: {}",
                    m_size_inner_list, size_obstacle_list);
    m_logger->debug("Inner starts with {} and ends with {}",
                    m_inner_list[0], m_inner_list[m_size_inner_list - 1]);
    m_logger->debug("--------------- END DOMAIN ---------------");
#endif
}

//======================================== Control =================================================
// *************************************************************************************************
/// \brief  Units test emergency solution
/// \param  size_obstacles Amount of obstacle cells
// *************************************************************************************************
void Domain::control(size_t size_obstacle_list, PatchObject &size_surface_list) {
#ifndef BENCHMARKING
    // TODO(n16h7): clean up
    std::string message;
    DomainData *domain = DomainData::getInstance();
    size_t size = domain->get_size(m_multigrid_level);

    size_t nx = domain->get_nx(m_multigrid_level);
    size_t ny = domain->get_ny(m_multigrid_level);
    size_t nz = domain->get_nz(m_multigrid_level);
    size_t Nx = domain->get_Nx(m_multigrid_level);
    size_t Ny = domain->get_Ny(m_multigrid_level);
    size_t all_cells = m_size_inner_list + m_size_boundary.get_sum();
    size_t duplicates = 4 * Nx + 4 * Ny + 4 * nz;
    if (m_size_domain_list != all_cells - duplicates) {
        message += fmt::format("list size of all domain cells does not fit with sum of it parts."
                               "Domain List: {} sum: {} duplicates: {}",
                               m_size_domain_list, all_cells, duplicates);
        for (size_t patch = 0; patch < number_of_patches; patch++) {
            message += fmt::format("{}: {}\n",
                                   PatchObject::get_patch_name(patch), m_size_boundary[patch]);
        }
    }
    if (m_size_domain_list + size_obstacle_list + size_surface_list.get_sum() != size) {
        message += fmt::format("list size of all domain cells is not equal with domain size."
                            " Domain List: {} Domain Size: {} Obstacle size: {} Surface size: {}",
                            m_size_domain_list, domain->get_size(m_multigrid_level),
                            size_obstacle_list, size_surface_list.get_sum());
    }
    size_t innerCells = nz * ny * nx;
    if (m_size_inner_list != innerCells - size_obstacle_list) {
        message += fmt::format("list size of inner cell does not equal domain inner size minus size of obstacles."
                            " Inner List: {} Domain inner size: {}, size obstacles: {}",
                            m_size_inner_list, innerCells, size_obstacle_list);
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
        m_logger->warn(message);
    }
#endif
}

Domain::~Domain() {
    delete[] m_inner_list;
    delete[] m_domain_list;
    delete[] m_boundary_list;
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        delete[] m_boundary_patch_divided[patch];
    }
    delete[] m_boundary_patch_divided;
}

//======================================== Boundary cells ==========================================
// *************************************************************************************************
/// \brief  Creates lists of indices of boundary cells
// *************************************************************************************************
void Domain::boundary_cells(size_t **surface_list, PatchObject &size_surface_list) {
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

    // fill boundary patches with boundary indices of computational domain (sorted)
    { // FRONT and BACK
        size_t counter_front = 0;
        size_t counter_back = 0;
        size_t counter_surface_front = 0;
        size_t counter_surface_back = 0;
        for (size_t j = y1; j <= y2; ++j) {
            for (size_t i = x1; i <= x2; ++i) {
                size_t idx = IX(i, j, z1, Nx, Ny);
                if (counter_surface_front >= size_surface_list[FRONT] || surface_list[FRONT][counter_surface_front] != idx) {
                    m_boundary_patch_divided[FRONT][counter_front++] = idx;
                } else {
                    counter_surface_front++;
                }

                idx = IX(i, j, z2, Nx, Ny);
                if (counter_surface_back >= size_surface_list[BACK] || surface_list[BACK][counter_surface_back] != idx) {
                    m_boundary_patch_divided[BACK][counter_back++] = idx;
                } else {
                    counter_surface_back++;
                }
            }
        }
    }
    { // TOP and BOTTOM
        size_t counter_bottom = 0;
        size_t counter_top = 0;
        size_t counter_surface_bottom = 0;
        size_t counter_surface_top = 0;
        for (size_t k = z1; k <= z2; ++k) {
            for (size_t i = x1; i <= x2; ++i) {
                size_t idx = IX(i, y1, k, Nx, Ny);
                if (counter_surface_bottom >= size_surface_list[BOTTOM] || surface_list[BOTTOM][counter_surface_bottom] != idx) {
                    m_boundary_patch_divided[BOTTOM][counter_bottom++] = idx;
                } else {
                    counter_surface_bottom++;
                }

                idx = IX(i, y2, k, Nx, Ny);
                if (counter_surface_top >= size_surface_list[TOP] || surface_list[TOP][counter_surface_top] != idx) {
                    m_boundary_patch_divided[TOP][counter_top++] = idx;
                } else {
                    counter_surface_top++;
                }
            }
        }
    }
    { // LEFT and RIGHT
        size_t counter_left = 0;
        size_t counter_right = 0;
        size_t counter_surface_left = 0;
        size_t counter_surface_right = 0;
        for (size_t k = z1; k <= z2; ++k) {
            for (size_t j = y1; j <= y2; ++j) {
                size_t idx = IX(x1, j, k, Nx, Ny);
                if (counter_surface_left >= size_surface_list[LEFT] || surface_list[LEFT][counter_surface_left] != idx) {
                    m_boundary_patch_divided[LEFT][counter_left++] = idx;
                } else {
                    counter_surface_left++;
                }

                idx = IX(x2, j, k, Nx, Ny);
                if (counter_surface_right >= size_surface_list[RIGHT] || surface_list[RIGHT][counter_surface_right] != idx) {
                    m_boundary_patch_divided[RIGHT][counter_right++] = idx;
                } else {
                    counter_surface_right++;
                }
            }
        }
    }
    // create boundary list via merging the six boundary patches
    std::vector<size_t> boundary_cells;
    boundary_cells.assign(m_boundary_patch_divided[0], m_boundary_patch_divided[0] + m_size_boundary[0]);
    for (size_t patch = 1; patch < number_of_patches; patch++) {
        //TODO(cvm): is that even legal?
        boundary_cells = Algorithm::merge_sort_with_duplicates(boundary_cells.data(), boundary_cells.size(), m_boundary_patch_divided[patch], m_size_boundary[patch]);
    }
    std::copy(boundary_cells.begin(), boundary_cells.end(), m_boundary_list);
}

//======================================== Inner cells =============================================
// *************************************************************************************************
/// \brief  Creates lists of indices of inner cells
/// \param  obstacle_list List of all obstacles of each multigrid level
/// \param  number_of_obstacles Amount of obstacles
// *************************************************************************************************
void Domain::inner_cells(const size_t *obstacle_list, size_t size_obstacle_list) {
    DomainData *domain = DomainData::getInstance();
    size_t k1 = domain->get_index_z1(m_multigrid_level);
    size_t j1 = domain->get_index_y1(m_multigrid_level);
    size_t i1 = domain->get_index_x1(m_multigrid_level);
    size_t k2 = domain->get_index_z2(m_multigrid_level);
    size_t j2 = domain->get_index_y2(m_multigrid_level);
    size_t i2 = domain->get_index_x2(m_multigrid_level);

    size_t Nx = domain->get_Nx(m_multigrid_level);
    size_t Ny = domain->get_Ny(m_multigrid_level);

    size_t counter_inner_cells = 0;
    size_t counter_obstacle_list = 0;
    for (size_t k = k1; k <= k2; ++k) {
        for (size_t j = j1; j <= j2; ++j) {
            for (size_t i = i1; i <= i2; ++i) {
                size_t idx = IX(i, j, k, Nx, Ny);
                if (counter_obstacle_list >= size_obstacle_list || obstacle_list[counter_obstacle_list] != idx) {
                    m_inner_list[counter_inner_cells++] = idx;
                } else {
                    counter_obstacle_list++;
                }
            }
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("inner list size: {} counter: {} obstacle counter: {}",
                    m_size_inner_list, counter_inner_cells, counter_obstacle_list);
#endif
}

//======================================== Update lists ============================================
// *************************************************************************************************
/// \brief  Updates lists of indices
/// \param  obstacle_list List of all obstacles of each multigrid level
/// \param  number_of_obstacles Number of obstacles
// *************************************************************************************************
void Domain::update_lists(size_t *obstacle_list, size_t size_obstacle_list,
                          size_t **surface_list, PatchObject &size_surface_list) {
    clear_lists();
    init(size_obstacle_list, size_surface_list);
    inner_cells(obstacle_list, size_obstacle_list);
    boundary_cells(surface_list, size_surface_list);
}

//======================================== Update lists ============================================
// *************************************************************************************************
/// \brief  removes all allocated arrays
// *************************************************************************************************
void Domain::clear_lists() {
    delete[] m_inner_list;
    delete[] m_domain_list;
    delete[] m_boundary_list;
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        delete[] m_boundary_patch_divided[patch];
    }
    delete[] m_boundary_patch_divided;
}
