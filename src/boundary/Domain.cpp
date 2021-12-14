/// \file       Boundary.cpp
/// \brief      Data class of boundary object
/// \date       Oct 01, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Domain.h"
#include <cmath>
#include <string>
#include "../DomainData.h"
#include "../utility/GlobalMacrosTypes.h"
#include "../utility/Algorithm.h"

Domain::Domain(size_t *obstacle_list, size_t size_obstacle_list,
               size_t **surface_list, PatchObject &size_surface_list,
               size_t multigrid_level) :
        m_multigrid_level(multigrid_level), m_size_boundary() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    init(size_obstacle_list, size_surface_list);
    inner_cells(obstacle_list, size_obstacle_list);
    boundary_cells(surface_list, size_surface_list);
    joined_list();

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

    m_boundary_patch_divided = new size_t *[number_of_patches];
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        size_t axis1 = (patch / 2 + 1) % number_of_axes;  // patch left (0) -> axis 1 (Y)
        size_t axis2 = (patch / 2 + 2) % number_of_axes;  // patch left (0) -> axis 2 (Z)
        m_size_boundary[patch] =
                domain_data->get_number_of_cells(CoordinateAxis(axis1), m_multigrid_level)
              * domain_data->get_number_of_cells(CoordinateAxis(axis2), m_multigrid_level)
              - size_surface_list[patch];

        m_boundary_patch_divided[patch] = new size_t[m_size_boundary[patch]];
    }

    m_size_inner_list = domain_data->get_nx(m_multigrid_level) * domain_data->get_ny(m_multigrid_level) * domain_data->get_nz(m_multigrid_level) - size_obstacle_list;
    m_inner_list = new size_t[m_size_inner_list];
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
        std::string patch_name = Mapping::get_patch_name(Patch(patch));
        m_logger->debug("list size of domain boundary {}: {} surface size: {}",
                        patch_name,
                        m_size_boundary[patch],
                        size_surface_list[patch]);
        m_logger->debug("{} starts with {} and ends with {}",
                        patch_name,
                        m_boundary_patch_divided[patch][0],
                        m_boundary_patch_divided[patch][m_size_boundary[patch] - 1]);
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
                               " Domain List: {} sum: {} duplicates: {}\n",
                               m_size_domain_list, all_cells, duplicates);
        for (size_t patch = 0; patch < number_of_patches; patch++) {
            message += fmt::format("  {}: {}\n",
                                   Mapping::get_patch_name(Patch(patch)), m_size_boundary[patch]);
        }
    }
    if (m_size_domain_list + size_obstacle_list + size_surface_list.get_sum() != size) {
        message += fmt::format("list size of all domain cells is not equal with domain size."
                            " Domain List: {} Domain Size: {} Obstacle size: {} Surface size: {}\n",
                            m_size_domain_list, domain->get_size(m_multigrid_level),
                            size_obstacle_list, size_surface_list.get_sum());
    }
    size_t innerCells = nz * ny * nx;
    if (m_size_inner_list != innerCells - size_obstacle_list) {
        message += fmt::format("list size of inner cell does not equal domain inner size minus size of obstacles."
                            " Inner List: {} Domain inner size: {}, size obstacles: {}\n",
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

    for (size_t i = 1; i < m_size_inner_list; i++) {
        int diff = static_cast<int>(m_inner_list[i] - m_inner_list[i - 1]);
        if (diff < 0) {
            message = message + "inner list: sorting error at index "
                      + std::to_string(i - 1) + "|"
                      + std::to_string(i) + " with values "
                      + std::to_string(m_inner_list[i - 1]) + "|"
                      + std::to_string(m_inner_list[i]) + "\n";
        }
    }
    for (size_t i = 1; i < m_size_boundary_list; i++) {
        int diff = static_cast<int>(m_boundary_list[i] - m_boundary_list[i - 1]);
        if (diff < 0) {
            message = message + "boundary list: sorting error at index "
                      + std::to_string(i - 1) + "|"
                      + std::to_string(i) + " with values "
                      + std::to_string(m_boundary_list[i - 1]) + "|"
                      + std::to_string(m_boundary_list[i]) + "\n";
        }
    }
    for (size_t i = 1; i < m_size_domain_list; i++) {
        int diff = static_cast<int>(m_domain_list[i] - m_domain_list[i - 1]);
        if (diff < 0) {
            message = message + "domain list: sorting error at index "
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
    auto domain_data = DomainData::getInstance();

    const size_t Nx = domain_data->get_Nx(m_multigrid_level);
    const size_t Ny = domain_data->get_Ny(m_multigrid_level);

    // DETAILED and CONCATENATED LISTS
    // BOUNDARY

    // TODO(issue 86): boundaries for physical domain.
    // TODO(issue 86): New method for computational domain -> redefine XML usage of boundaries

    auto tmp_start = new Coordinate<size_t>();
    auto tmp_end = new Coordinate<size_t>();
    // comments for X axis (0)
    for (size_t axis = 0; axis < number_of_axes; axis++) {
        CoordinateAxis axis1;
        CoordinateAxis axis2;
        if (CoordinateAxis(axis) == CoordinateAxis::X) {
            axis1 = CoordinateAxis::Y;
            axis2 = CoordinateAxis::Z;
        } else if (CoordinateAxis(axis) == CoordinateAxis::Y) {
            axis1 = CoordinateAxis::X;
            axis2 = CoordinateAxis::Z;
        } else {
            axis1 = CoordinateAxis::X;
            axis2 = CoordinateAxis::Y;
        }

        auto patch_start = Patch(axis * 2);  // Patch Left (0)
        auto patch_end = Patch(axis * 2 + 1);  // Patch Right (1)

        // start indices for computational domain minus 1 for ghost cells
        size_t z1 = domain_data->get_start_index_CD(axis2, m_multigrid_level) - 1;
        size_t y1 = domain_data->get_start_index_CD(axis1, m_multigrid_level) - 1;
        (*tmp_start)[axis] = domain_data->get_start_index_CD(CoordinateAxis(axis), m_multigrid_level) - 1;

        // end indices for computational domain plus 1 for ghost cells
        size_t z2 = domain_data->get_end_index_CD(axis2, m_multigrid_level) + 1;
        size_t y2 = domain_data->get_end_index_CD(axis1, m_multigrid_level) + 1;
        (*tmp_end)[axis] = domain_data->get_end_index_CD(CoordinateAxis(axis), m_multigrid_level) + 1;

        size_t counter_start = 0;
        size_t counter_end = 0;
        size_t counter_surface_start = 0;
        size_t counter_surface_end = 0;
        for (size_t k = z1; k <= z2; k++) {
            (*tmp_start)[axis2] = k;
            (*tmp_end)[axis2] = k;
            for (size_t j = y1; j <= y2; j++) {
                (*tmp_start)[axis1] = j;
                size_t idx = tmp_start->get_index(Nx, Ny);
                if (counter_surface_start >= size_surface_list[patch_start] || surface_list[patch_start][counter_surface_start] != idx) {
                    m_boundary_patch_divided[patch_start][counter_start++] = idx;
                } else {
                    counter_surface_start++;
                }

                (*tmp_end)[axis1] = j;
                idx = tmp_end->get_index(Nx, Ny);
                if (counter_surface_end >= size_surface_list[patch_end] || surface_list[patch_end][counter_surface_end] != idx) {
                    m_boundary_patch_divided[patch_end][counter_end++] = idx;
                } else {
                    counter_surface_end++;
                }
            }
        }
#ifndef BENCHMARKING
        m_logger->debug("for patch {} boundary list size: {} counter: {} surface list size: {} surface counter: {}",
                        Mapping::get_patch_name(patch_start),
                        m_size_boundary[patch_start], counter_start,
                        size_surface_list[patch_start], counter_surface_start);
        m_logger->debug("for patch {} boundary list size: {} counter: {} surface list size: {} surface counter: {}",
                        Mapping::get_patch_name(patch_end),
                        m_size_boundary[patch_end], counter_end,
                        size_surface_list[patch_end], counter_surface_end);
#endif
    }
    delete tmp_start;
    delete tmp_end;
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
                if (counter_obstacle_list >= size_obstacle_list) {
                    m_inner_list[counter_inner_cells++] = idx;
                } else if (obstacle_list[counter_obstacle_list] == idx){
                    counter_obstacle_list++;
                } else {
                    m_inner_list[counter_inner_cells++] = idx;
                }
            }
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("inner list size: {} counter: {} obstacle list size: {} obstacle counter: {}",
                    m_size_inner_list, counter_inner_cells, size_obstacle_list, counter_obstacle_list);
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
    //TODO (issue 86) it may be not necessary to update every list
    clear_lists();
    init(size_obstacle_list, size_surface_list);
    inner_cells(obstacle_list, size_obstacle_list);
    boundary_cells(surface_list, size_surface_list);
    joined_list();
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

void Domain::joined_list() {
    // create boundary list via merging the six boundary patches
    std::vector<size_t> boundary_cells;
    boundary_cells.assign(m_boundary_patch_divided[0], m_boundary_patch_divided[0] + m_size_boundary[0]);
    for (size_t patch = 1; patch < number_of_patches; patch++) {
        //TODO(cvm): is that even legal?
        boundary_cells = Algorithm::merge_sort_with_duplicates(boundary_cells.data(), boundary_cells.size(), m_boundary_patch_divided[patch], m_size_boundary[patch]);
    }
    m_size_boundary_list = boundary_cells.size();
    m_boundary_list = new size_t[m_size_boundary_list];
    std::copy(boundary_cells.begin(), boundary_cells.end(), m_boundary_list);

    m_size_domain_list = m_size_inner_list + m_size_boundary_list;
    m_domain_list = new size_t[m_size_domain_list];
    Algorithm::merge_sort(m_inner_list, m_boundary_list,
                          m_size_inner_list, m_size_boundary_list,
                          m_domain_list);
}
