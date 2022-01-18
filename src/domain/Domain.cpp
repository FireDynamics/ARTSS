/// \file       Domain.cpp
/// \brief      stores index lists of domain
/// \date       Oct 01, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Domain.h"
#include <cmath>
#include <string>
#include "DomainData.h"
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
#endif
}


//======================================== Init ====================================================
// *************************************************************************************************
/// \brief  Initialize member variables (arrays)
/// \param  size_obstacles Amount of obstacle cells
// *************************************************************************************************
void Domain::init(size_t size_obstacle_list, PatchObject &size_surface_list) {
    auto domain_data = DomainData::getInstance();

    m_boundary_patch_divided.resize(number_of_patches);
    for (Patch patch: {Patch::LEFT, Patch::RIGHT, Patch::BOTTOM, Patch::TOP, Patch::FRONT, Patch::BACK}) {
        auto axes = Mapping::get_axes(patch);
        m_size_boundary[patch] =
                (domain_data->get_number_of_inner_cells(axes[0], m_multigrid_level) + 2) *
                (domain_data->get_number_of_inner_cells(axes[1], m_multigrid_level) + 2) -
                size_surface_list[patch];

        m_boundary_patch_divided[patch].resize(m_size_boundary[patch]);
    }

    m_inner_list.resize(
            domain_data->get_number_of_inner_cells(X, m_multigrid_level) *
            domain_data->get_number_of_inner_cells(Y, m_multigrid_level) *
            domain_data->get_number_of_inner_cells(Z, m_multigrid_level) -
            size_obstacle_list);
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
    m_logger->debug("list size of domain list: {}", m_domain_list.size());
    m_logger->debug("Domain starts with {} and ends with {}",
                    m_domain_list[0], m_domain_list[m_domain_list.size() - 1]);
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
                    m_inner_list.size(), size_obstacle_list);
    m_logger->debug("Inner starts with {} and ends with {}",
                    m_inner_list[0], m_inner_list[m_inner_list.size() - 1]);
    m_logger->debug("--------------- END DOMAIN ---------------");
#endif
}

//======================================== Boundary cells ==========================================
// *************************************************************************************************
/// \brief  Creates lists of indices of boundary cells
// *************************************************************************************************
void Domain::boundary_cells(size_t **surface_list, PatchObject &size_surface_list) {
    auto domain_data = DomainData::getInstance();

    const size_t Nx = domain_data->get_number_of_cells(X, m_multigrid_level);
    const size_t Ny = domain_data->get_number_of_cells(Y, m_multigrid_level);

    // DETAILED and CONCATENATED LISTS
    // BOUNDARY

    // TODO(issue 86): boundaries for physical domain.
    // TODO(issue 86): New method for computational domain -> redefine XML usage of boundaries

    Coordinate<size_t> tmp_start;
    Coordinate<size_t> tmp_end;
    // comments for X axis (0)
    for (CoordinateAxis axis: {X, Y, Z}) {
        CoordinateAxis axis1;
        CoordinateAxis axis2;
        if (axis == CoordinateAxis::X) {
            axis1 = CoordinateAxis::Y;
            axis2 = CoordinateAxis::Z;
        } else if (axis == CoordinateAxis::Y) {
            axis1 = CoordinateAxis::X;
            axis2 = CoordinateAxis::Z;
        } else {
            axis1 = CoordinateAxis::X;
            axis2 = CoordinateAxis::Y;
        }

        auto patch_start = Mapping::to_patch(axis, true);  // Patch Left (0)
        auto patch_end = Mapping::to_patch(axis, false);  // Patch Right (1)

        // start indices for computational domain minus 1 for ghost cells
        size_t z1 = domain_data->get_start_index_CD(axis2, m_multigrid_level) - 1;
        size_t y1 = domain_data->get_start_index_CD(axis1, m_multigrid_level) - 1;
        tmp_start[axis] = domain_data->get_start_index_CD(CoordinateAxis(axis), m_multigrid_level) - 1;

        // end indices for computational domain plus 1 for ghost cells
        size_t z2 = domain_data->get_end_index_CD(axis2, m_multigrid_level) + 1;
        size_t y2 = domain_data->get_end_index_CD(axis1, m_multigrid_level) + 1;
        tmp_end[axis] = domain_data->get_end_index_CD(CoordinateAxis(axis), m_multigrid_level) + 1;

        size_t counter_start = 0;
        size_t counter_end = 0;
        size_t counter_surface_start = 0;
        size_t counter_surface_end = 0;
        for (size_t k = z1; k <= z2; k++) {
            tmp_start[axis2] = k;
            tmp_end[axis2] = k;
            for (size_t j = y1; j <= y2; j++) {
                tmp_start[axis1] = j;
                size_t idx = tmp_start.get_index(Nx, Ny);
                if (counter_surface_start >= size_surface_list[patch_start] || surface_list[patch_start][counter_surface_start] != idx) {
                    m_boundary_patch_divided[patch_start][counter_start++] = idx;
                } else {
                    counter_surface_start++;
                }

                tmp_end[axis1] = j;
                idx = tmp_end.get_index(Nx, Ny);
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
}

//======================================== Inner cells =============================================
// *************************************************************************************************
/// \brief  Creates lists of indices of inner cells
/// \param  obstacle_list List of all obstacles of each multigrid level
/// \param  number_of_obstacles Amount of obstacles
// *************************************************************************************************
void Domain::inner_cells(const size_t *obstacle_list, size_t size_obstacle_list) {
    auto domain_data = DomainData::getInstance();
    size_t k1 = domain_data->get_start_index_CD(Z, m_multigrid_level);
    size_t j1 = domain_data->get_start_index_CD(Y, m_multigrid_level);
    size_t i1 = domain_data->get_start_index_CD(X, m_multigrid_level);
    size_t k2 = domain_data->get_end_index_CD(Z, m_multigrid_level);
    size_t j2 = domain_data->get_end_index_CD(Y, m_multigrid_level);
    size_t i2 = domain_data->get_end_index_CD(X, m_multigrid_level);

    size_t Nx = domain_data->get_number_of_cells(X, m_multigrid_level);
    size_t Ny = domain_data->get_number_of_cells(Y, m_multigrid_level);

    size_t counter_inner_cells = 0;
    size_t counter_obstacle_list = 0;
    for (size_t k = k1; k <= k2; ++k) {
        for (size_t j = j1; j <= j2; ++j) {
            for (size_t i = i1; i <= i2; ++i) {
                size_t idx = IX(i, j, k, Nx, Ny);
                if (counter_obstacle_list < size_obstacle_list && obstacle_list[counter_obstacle_list] == idx) {
                    counter_obstacle_list++;
                } else {
                    m_inner_list[counter_inner_cells++] = idx;
                }
            }
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("inner list size: {} counter: {} obstacle list size: {} obstacle counter: {}",
                    m_inner_list.size(), counter_inner_cells, size_obstacle_list, counter_obstacle_list);
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
    m_inner_list.clear();
    m_domain_list.clear();
    m_boundary_list.clear();
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_boundary_patch_divided[patch].clear();
    }
}

void Domain::joined_list() {
    // create boundary list via merging the six boundary patches
    std::vector<size_t> boundary_cells;
    boundary_cells.assign(m_boundary_patch_divided[0].begin(), m_boundary_patch_divided[0].end());
    for (size_t patch = 1; patch < number_of_patches; patch++) {
        //TODO(cvm): is that even legal?
        boundary_cells = Algorithm::merge_sort_with_duplicates(boundary_cells.data(), boundary_cells.size(), m_boundary_patch_divided[patch].data(), m_size_boundary[patch]);
    }
    m_boundary_list = boundary_cells;
    m_domain_list = Algorithm::merge_sort(m_inner_list.data(), m_boundary_list.data(),
                                          m_inner_list.size(), m_boundary_list.size());
}
