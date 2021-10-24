/// \file       Multigrid.h
/// \brief      Creates all lists needed for multigrid
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Multigrid.h"
#include "PatchObject.h"
#include <string>
#include <vector>

Multigrid::Multigrid(BoundaryDataController *bdc_boundary, size_t multigrid_levels) :
        m_multigrid_levels(multigrid_levels),
        m_jl_domain_inner_list(multigrid_levels),
        m_jl_domain_boundary_list(multigrid_levels),
        m_jl_obstacle_list(multigrid_levels),
        m_bdc_boundary(bdc_boundary) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
#ifdef GPU_DEBUG
    m_gpu_logger = Utility::create_gpu_logger(typeid(this).name());
#endif
    init_domain();
    add_MG_lists();
    send_lists_to_GPU();
#ifndef BENCHMARKING
    print();
    control();
#endif
}

Multigrid::Multigrid(
        size_t number_of_surfaces, Surface **surface_list,
        size_t number_of_obstacles, Obstacle **obstacle_list,
        BoundaryDataController *bdc_boundary,
        BoundaryDataController **bdc_obstacles,
        size_t multigrid_levels) :
        m_multigrid_levels(multigrid_levels),
        m_number_of_surface_objects(number_of_surfaces),
        m_number_of_obstacle_objects(number_of_obstacles),
        m_jl_domain_inner_list(multigrid_levels),
        m_jl_domain_boundary_list(multigrid_levels),
        m_jl_obstacle_list(multigrid_levels),
        m_bdc_boundary(bdc_boundary),
        m_bdc_obstacle(bdc_obstacles) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
#ifdef GPU_DEBUG
    m_gpu_logger = Utility::create_gpu_logger(typeid(this).name());
#endif
    init_domain();

    if (m_number_of_surface_objects > 0) {
        init_surfaces(surface_list);
    }

    if (m_number_of_obstacle_objects > 0) {
        init_obstacles(obstacle_list);
    }

    add_MG_lists();
    send_lists_to_GPU();
#ifndef BENCHMARKING
    print();
    control();
#endif
}

//======================================== Init ====================================================
// *************************************************************************************************
/// \brief  Initialize member variables for domain
// *************************************************************************************************
void Multigrid::init_domain() {
    // list of domain boundary for each level
    m_MG_domain_object_list = new Domain *[m_multigrid_levels + 1];

    m_jl_domain_boundary_list_patch_divided = new SimpleJoinedList*[number_of_patches];
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_jl_domain_boundary_list_patch_divided[patch]  = new SimpleJoinedList(m_multigrid_levels);
    }
}

void Multigrid::init_obstacles(Obstacle **obstacle_list){
    // list of obstacles for each level
    m_MG_obstacle_object_list = new Obstacle **[m_multigrid_levels + 1];
    m_MG_obstacle_object_list[0] = obstacle_list;  // level 0

    m_jl_obstacle_boundary_list_patch_divided = new ObstacleJoinedList*[number_of_patches];
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_jl_obstacle_boundary_list_patch_divided[patch]  = new ObstacleJoinedList(m_multigrid_levels, m_number_of_obstacle_objects);
    }
}

void Multigrid::init_surfaces(Surface **surface_list) {
    // list of surfaces for each level
    m_MG_surface_object_list = new Surface **[m_multigrid_levels + 1];
    m_MG_surface_object_list[0] = surface_list;  // level 0

    // surface indices divided by level
    m_MG_surface_index_list = new size_t *[m_multigrid_levels + 1];

    // start index of each surface in level joined list
    m_size_MG_surface_index_list_level = new size_t[m_multigrid_levels * m_number_of_surface_objects + 1];
    // start index of first level in joined list = 0
    *(m_size_MG_surface_index_list_level) = 0;
}

Multigrid::~Multigrid() {
    //TODO(issue 130)
    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        if (m_number_of_surface_objects > 0) {
            Surface **surface_level = *(m_MG_surface_object_list + level);
            for (size_t surface = 0; surface < m_number_of_surface_objects; surface++) {
                delete (*(surface_level + surface));
            }
            delete[] surface_level;
            delete (*(m_MG_surface_index_list + level));
        }
        if (m_number_of_obstacle_objects > 0) {
            Obstacle **obstacle_level = *(m_MG_obstacle_object_list);
            for (size_t obstacle = 0; obstacle < m_number_of_obstacle_objects; obstacle++) {
                delete (*(obstacle_level + obstacle));
            }
            delete[] obstacle_level;
        }
        delete (*(m_MG_domain_object_list + level));
    }
    delete[] m_MG_surface_object_list;
    delete[] m_MG_domain_object_list;
    delete[] m_MG_obstacle_object_list;

    if (m_number_of_surface_objects > 0) {
        size_t size_surface_list = get_length_of_surface_index_list_joined();
#pragma acc exit data delete(m_data_MG_surface_list_level_joined[:size_surface_list])
        delete[] m_MG_surface_object_list;
        delete[] m_MG_surface_index_list;
        delete[] m_size_MG_surface_index_list_level;
        delete[] m_data_MG_surface_list_level_joined;
    }
    if (m_number_of_obstacle_objects > 0) {
        delete[] m_MG_obstacle_object_list;
    }

    for (size_t patch = 0; patch < number_of_patches; patch++) {
        delete m_jl_obstacle_boundary_list_patch_divided[patch];
        delete m_jl_domain_boundary_list_patch_divided[patch];
    }
}

//======================================== Control =================================================
// *************************************************************************************************
/// \brief  Units test emergency solution
// *************************************************************************************************
void Multigrid::control() {
#ifndef BENCHMARKING
    std::string message;
    auto domain = DomainData::getInstance();
    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        size_t nx = domain->get_nx(level);
        size_t ny = domain->get_ny(level);
        size_t nz = domain->get_nz(level);

        if (!message.empty()) {
            message += "For Level " + std::to_string(level) + "\n";
            message += "size control\n";
        }
        size_t original_len = (static_cast<Domain *>(m_MG_domain_object_list[level]))->get_size_inner_list();
        size_t calculated_len = m_jl_domain_inner_list.get_last_index(level) - m_jl_domain_inner_list.get_first_index(level) + 1;
        size_t saved_len = m_jl_domain_inner_list.get_slice_size(level);
        if (calculated_len != original_len) {
            size_t control = domain->get_nx(level) * domain->get_ny(level) * domain->get_nz(level);
            message += fmt::format("length calculated/stored of domain inner cells does not equals its original size size:"
                                   "original: {} saved: {} calculated: {} control: {}\n",
                                   original_len, saved_len, calculated_len, control);
        }
        original_len = (static_cast<Domain *>(m_MG_domain_object_list[level]))->get_size_domain_list();
        calculated_len = m_jl_domain_boundary_list.get_last_index(level) - m_jl_domain_boundary_list.get_first_index(level) + 1;
        saved_len = m_jl_domain_boundary_list.get_slice_size(level);
        if (calculated_len != original_len) {
            size_t control = (domain->get_Nx(level) * domain->get_Ny(level) * 2)
                             + (domain->get_Nx(level) * (domain->get_Nz(level) - 2)) * 2
                             + ((domain->get_Ny(level) - 2) * (domain->get_Nz(level) - 2)) * 2;
            message += fmt::format("length calculated/stored of domain boundary cells does not equals its original size size:"
                                   "original: {} saved: {} calculated: {} control: {}\n",
                                   original_len, saved_len, calculated_len, control);
        }
        for (size_t patch = 0; patch < number_of_patches; patch++) {
            original_len = (static_cast<Domain *>(m_MG_domain_object_list[level]))->get_size_boundary_list()[patch];
            saved_len = m_jl_domain_boundary_list_patch_divided[patch]->get_slice_size(level);
            calculated_len = m_jl_domain_boundary_list_patch_divided[patch]->get_first_index(level) - m_jl_domain_boundary_list_patch_divided[patch]->get_slice_size(level + 1);
            if (calculated_len != original_len || saved_len != original_len) {
                size_t control = domain->get_ny(level) * domain->get_nz(level);
                message += fmt::format("length calculated/stored of domain boundary '{}' does not equals its original size size:"
                                       "original: {} saved: {} calculated: {} control: {}\n",
                                       BoundaryData::get_patch_name(static_cast<Patch>(patch)),
                                       original_len, saved_len, calculated_len, control);
            }
        }

        size_t csize_inner = domain->get_nx(level) * domain->get_ny(level) * domain->get_nz(level)
                             - m_jl_obstacle_list.get_slice_size(level);
        size_t bsize_inner = get_slice_size_domain_inner_cells_level_joined(level);
        if (csize_inner != bsize_inner) {
            message += "get_size_domain_inner_cells(level) does not equal (nx-2)*(ny-2)*(nz-2) "
                       + std::to_string(bsize_inner) + "|" + std::to_string(csize_inner) + "\n";
        }
        size_t cindex_inner_start = get_start_index_domain_inner_cells_level_joined(level);
        size_t cindex_inner_end = get_end_index_domain_inner_cells_level_joined(level);
        if (cindex_inner_end - cindex_inner_start + 1 != bsize_inner) {
            message += "get_size_domain_inner_cells(level) does not equal the difference between start and end "
                       + std::to_string(cindex_inner_start) + "|" + std::to_string(cindex_inner_end) + "\n";
        }
        size_t bsize_boundary = domain->get_size(level) - nx * ny * nz;
        size_t csize_boundary = m_jl_domain_boundary_list.get_slice_size(level);
        if (csize_boundary != bsize_boundary) {
            message += "get_size_domain_list(level) does not equal size-(nx-2)*(ny-2)*(nz-2) "
                       + std::to_string(bsize_boundary) + "|" + std::to_string(csize_boundary) + "\n";
        }

        if (!message.empty()) {
            message = "For Level " + std::to_string(level) + "\nsize control\n" + message;
        }
    }
    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        Domain *b = *(m_MG_domain_object_list + level);
        b->control(m_jl_obstacle_list.get_slice_size(level));
    }
    {
        size_t bsize_inner = 0;
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            bsize_inner += get_slice_size_domain_inner_cells_level_joined(level);
        }
        size_t csize_inner = get_size_domain_inner_cells_level_joined();
        if (bsize_inner != csize_inner) {
            message += "get_size_domain_inner_cells_level_joined does not equal the sum of each inner list "
                       + std::to_string(bsize_inner) + "|" + std::to_string(csize_inner) + "\n";
        }
    }

    if (!message.empty()) {
        message = "################ MULTIGRID CONTROL ################\n" + message
                  + "---------------- MULTIGRID CONTROL END ----------------";
        m_logger->warn(message);
    }
#endif
}

//======================================== Print ===================================================
// *************************************************************************************************
/// \brief  Print multigrid infos
// *************************************************************************************************
void Multigrid::print() {
#ifndef BENCHMARKING
    m_logger->debug("################ MULTIGRID ################");
    m_logger->debug("Number of Obstacles: {}, Number of Surfaces: {}",
                    m_number_of_obstacle_objects, m_number_of_surface_objects);
    m_logger->debug("Multigrid levels: {}", m_multigrid_levels);

    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        m_logger->debug("For Level {} domain inner list starts at index {} and ends with index {}", level,
                        m_jl_domain_inner_list.get_first_index(level),
                        m_jl_domain_inner_list.get_last_index(level));
        m_logger->debug("and the corresponding indices at this position: {}|{}",
                        m_jl_domain_inner_list[m_jl_domain_inner_list.get_first_index(level)],
                        m_jl_domain_inner_list[m_jl_domain_inner_list.get_last_index(level)]);
    }
    m_logger->debug("Total number of domain inner cells: {}", m_jl_domain_inner_list.get_size());

    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        m_logger->debug("For Level {} domain boundary list starts at index {} and ends with index {}", level,
                        m_jl_domain_boundary_list.get_first_index(level),
                        m_jl_domain_boundary_list.get_last_index(level));
        m_logger->debug("and the corresponding indices at this position: {} | {}",
                        m_jl_domain_boundary_list.get_data()[m_jl_domain_boundary_list.get_first_index(level)],
                        m_jl_domain_boundary_list.get_data()[m_jl_domain_boundary_list.get_last_index(level)]);
    }
    m_logger->debug("Total number of domain boundary cells: {}", m_jl_domain_boundary_list.get_size());

    for (size_t patch = 0; patch < number_of_patches; patch++) {
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            size_t start_index = m_jl_domain_boundary_list_patch_divided[patch]->get_first_index(level);
            size_t end_index = m_jl_domain_boundary_list_patch_divided[patch]->get_last_index(level);
            m_logger->debug("For Level {} domain boundary '{}' starts at index {} and ends with index {}",
                            level, BoundaryData::get_patch_name(patch),
                            start_index, end_index);
            m_logger->debug(" and the corresponding indices at this position {} | {}",
                            m_jl_domain_boundary_list_patch_divided[patch]->get_data()[start_index],
                            m_jl_domain_boundary_list_patch_divided[patch]->get_data()[end_index]);
        }
        m_logger->debug("Total number of '{}' boundary cells: {}",
                        BoundaryData::get_patch_name(patch),
                        m_jl_domain_boundary_list_patch_divided[patch]->get_size());
    }

    for (size_t level = 0; level < m_multigrid_levels; level++) {
        m_logger->debug("For Level {} obstacle_list has {} elements",
                        level, m_jl_obstacle_list.get_slice_size(level));
    }

    for (size_t patch = 0; patch < number_of_patches; patch++) {
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            m_logger->debug("For Level {} obstacle '{}' starts at index {} and ends with index {} with length {}",
                            level,
                            BoundaryData::get_patch_name(patch),
                            m_jl_obstacle_boundary_list_patch_divided[patch]->get_first_index(level, 0),
                            m_jl_obstacle_boundary_list_patch_divided[patch]->get_last_index(level, m_number_of_obstacle_objects - 1),
                            m_jl_obstacle_boundary_list_patch_divided[patch]->get_slice_size(level));
        }
        m_logger->debug("Total number of '{}' obstacle cells: {}",
                        BoundaryData::get_patch_name(patch),
                        m_jl_obstacle_boundary_list_patch_divided[patch]->get_size());
    }

    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        m_logger->debug("For Level {} surface_list starts at index {} and ends with index {}",
                        level,
                        get_first_index_of_surface_index_list(level),
                        get_last_index_of_surface_index_list(level));
    }
    m_logger->debug("Total length of surface_list: {}", get_length_of_surface_index_list_joined());
    m_logger->debug("---------------- MULTIGRID END ----------------");
#endif
}

// ================================= Add MG lists ==================================================
// *************************************************************************************************
/// \brief  adds lists (outer, inner, surfaces, obstacles) in case of grid restriction (dominant)
// *************************************************************************************************
void Multigrid::add_MG_lists() {
    size_t sum_domain_inner_cells = 0;
    auto * sum_domain_boundary = new PatchObject();
    size_t sum_obstacle_cells = 0;
    auto *obstacle_cells = new size_t[m_multigrid_levels];
    auto sum_obstacle_boundary = new PatchObject();

    {  // create domain object of level 0
        Domain *domain;
        if (m_number_of_obstacle_objects > 0) {
            sum_obstacle_cells += calc_obstacles(m_MG_obstacle_object_list[0], sum_obstacle_boundary);
            obstacle_cells[0] = sum_obstacle_cells;
            domain = new Domain(m_MG_obstacle_object_list[0], m_number_of_obstacle_objects,
                                m_jl_obstacle_list.get_slice_size(0));
        } else {
            domain = new Domain();
        }
        sum_domain_inner_cells += domain->get_size_inner_list();
        //m_jl_inner_list.add_data(0, domain->get_size_inner_list(), domain->get_inner_list());
        PatchObject &boundary_size = domain->get_size_boundary_list();
        *sum_domain_boundary += boundary_size;

        // save boundary object in multigrid list
        m_MG_domain_object_list[0] = domain;
    }

    if (m_number_of_surface_objects > 0) {
        calc_surfaces(*(m_MG_surface_object_list));
    }

    // create boundary object, surfaces and obstacles for each multigrid level
    size_t **tmp_store_obstacle = new size_t*[m_multigrid_levels];
    for (size_t level = 1; level < m_multigrid_levels + 1; level++) {
        surface_dominant_restriction(level);
        size_t sum = obstacle_dominant_restriction(level, sum_domain_boundary, tmp_store_obstacle);
        sum_obstacle_cells += sum;
        obstacle_cells[level] = sum;

        Domain *domain;
        if (m_number_of_obstacle_objects > 0) {
            domain = new Domain(m_MG_obstacle_object_list[level], m_number_of_obstacle_objects,
                                m_jl_obstacle_list.get_slice_size(level), level);
        } else {
            domain = new Domain(level);
        }
        sum_domain_inner_cells += domain->get_size_inner_list();
        *sum_domain_boundary += domain->get_size_boundary_list();

        m_MG_domain_object_list[level] = domain;
    }
    m_jl_domain_inner_list.set_size(sum_domain_inner_cells);
    m_jl_obstacle_list.set_size(sum_obstacle_cells);

    for (size_t patch = 0 ; patch < number_of_patches; patch++) {
        m_jl_domain_boundary_list_patch_divided[patch]->set_size((*sum_domain_boundary)[patch]);
        m_jl_obstacle_boundary_list_patch_divided[patch]->set_size((*sum_obstacle_boundary)[patch]);
    }

    m_jl_obstacle_list.set_size(sum_obstacle_cells);
    for (size_t level = 0; level < m_multigrid_levels; level++) {
        m_jl_obstacle_list.add_data(level, obstacle_cells[level], tmp_store_obstacle[level]);
    }
}

// ================================= Calc obstacles ================================================
// *************************************************************************************************
/// \brief  create obstacles of level 0, no check if obstacles are overlapping
/// \param obstacle_object_list List of obstacle objects
// *************************************************************************************************
size_t Multigrid::calc_obstacles(Obstacle **obstacle_object_list, PatchObject *patch_sum) const {
    size_t sum_obstacle_list = 0;
    for (size_t id = 0; id < m_number_of_obstacle_objects; id++) {
        *patch_sum += obstacle_object_list[id]->get_size_boundary_list();
        sum_obstacle_list += obstacle_object_list[id]->get_size_obstacle_list();
    }
    return sum_obstacle_list;
}

// ================================= Calc surfaces =================================================
// *************************************************************************************************
/// \brief  Create surfaces of level 0
/// \param surface_list List of surface objects
// *************************************************************************************************
void Multigrid::calc_surfaces(Surface **surface_object_list) {
    if (m_number_of_surface_objects > 0) {
        size_t *surface_index_list = new size_t[get_first_index_of_surface_index_list(1)];
        size_t counter_s = 0;
        for (size_t s = 0; s < m_number_of_surface_objects; s++) {
            Surface *surface_tmp = surface_object_list[s];
            for (size_t i = 0; i < surface_tmp->getSize_surfaceList(); i++) {
                *(surface_index_list + counter_s) = surface_tmp->getSurfaceList()[i];
                counter_s++;
            }
            *(m_size_MG_surface_index_list_level + s + 1) = counter_s;
        }
        *(m_MG_surface_index_list) = surface_index_list;
    }
}

// ================================= Surface dominant restriction ==================================
// *************************************************************************************************
/// \brief  Calculates coarse level of multigrid of surfaces (dominant restriction)
/// \param level Multigrid level
// *************************************************************************************************
void Multigrid::surface_dominant_restriction(size_t level) {
// SURFACES
    if (m_number_of_surface_objects > 0) {
        // add index to m_surface_list if any of l-1 indices building the l index was a surface
        // dominant restriction)
        DomainData *domain = DomainData::getInstance();
        size_t Nx = domain->get_Nx(level);
        size_t Ny = domain->get_Ny(level);

        Surface **surface_list_fine = *(m_MG_surface_object_list + (level - 1));
        Surface **surface_list_coarse = new Surface *[m_number_of_surface_objects];
        *(m_MG_surface_object_list + level) = surface_list_coarse;
        // loop through surfaces
        for (size_t surfaceID = 0; surfaceID < m_number_of_surface_objects; ++surfaceID) {
            Surface *surface_fine = surface_list_fine[surfaceID];

            size_t stride_x_fine = surface_fine->getStrideX();
            size_t stride_y_fine = surface_fine->getStrideY();
            size_t stride_z_fine = surface_fine->getStrideZ();
            size_t startIndex_fine = surface_fine->getSurfaceList()[0];
            size_t k_fine = getCoordinateK(startIndex_fine, Nx, Ny);
            size_t j_fine = getCoordinateJ(startIndex_fine, Nx, Ny, k_fine);
            size_t i_fine = getCoordinateI(startIndex_fine, Nx, Ny, j_fine, k_fine);

            size_t stride_x_coarse, stride_y_coarse, stride_z_coarse;
            if (i_fine % 2 == 0) {
                stride_x_coarse = stride_x_fine / 2 + 1;
            } else {
                stride_x_coarse = (stride_x_fine + 1) / 2;
            }
            if (j_fine % 2 == 0) {
                stride_y_coarse = stride_y_fine / 2;
            } else {
                stride_y_coarse = (stride_y_fine + 1) / 2;
            }
            if (k_fine % 2 == 0) {
                stride_z_coarse = stride_z_fine / 2 + 1;
            } else {
                stride_z_coarse = (stride_z_fine + 1) / 2;
            }
            size_t startIndex_coarse = IX(i_fine / 2, j_fine / 2, k_fine / 2, Nx, Ny);
#ifndef BENCHMARKING
            m_logger->debug("startIndex multigrid surface: {} {}|{}",
                            startIndex_fine,
                            startIndex_coarse,
                            startIndex_fine / 2);
#endif

            auto *surface_coarse = new Surface(surfaceID, startIndex_coarse,
                                                  stride_x_coarse, stride_y_coarse, stride_z_coarse, level);
            *(surface_list_coarse + surfaceID) = surface_coarse;
            size_t index = level * m_number_of_surface_objects + surfaceID + 1;
            *(m_size_MG_surface_index_list_level + index) = *(m_size_MG_surface_index_list_level + index - 1)
                                                            + surface_coarse->getSize_surfaceList();
#ifndef BENCHMARKING
            m_logger->debug("control multigrid surface index: {} {}",
                            *(m_size_MG_surface_index_list_level + index),
                            surface_coarse->getSize_surfaceList());
#endif
        }  // end surface id loop
    }
}

// ================================= Obstacle dominant restriction =================================
// *************************************************************************************************
/// \brief  Calculates coarse level of multigrid of obstacles (dominant restriction)
/// \param level Multigrid level
/// \return Obstacle** list of created obstacles
// *************************************************************************************************
size_t Multigrid::obstacle_dominant_restriction(size_t level, PatchObject *sum_patches, size_t **tmp_store_obstacle) {
    // OBSTACLES
    // add index to m_obstacle_lists if any of l-1 indices building the l index was an obstacle
    // (dominant restriction)
    if (m_number_of_obstacle_objects == 0) {
        return 0;
    }
    DomainData *domain_data = DomainData::getInstance();
    Obstacle **obstacle_list_fine = m_MG_obstacle_object_list[level - 1];
    auto obstacle_list_coarse = new Obstacle *[m_number_of_obstacle_objects];
    m_MG_obstacle_object_list[level] = obstacle_list_coarse;
    for (size_t id = 0; id < m_number_of_obstacle_objects; id++) {
        Obstacle *obstacle_fine = obstacle_list_fine[id];
        size_t i1_fine = obstacle_fine->get_coordinates_i1();
        size_t j1_fine = obstacle_fine->get_coordinates_j1();
        size_t k1_fine = obstacle_fine->get_coordinates_k1();
        size_t i2_fine = obstacle_fine->get_coordinates_i2();
        size_t j2_fine = obstacle_fine->get_coordinates_j2();
        size_t k2_fine = obstacle_fine->get_coordinates_k2();

        size_t i1_coarse = (i1_fine + 1) / 2;
        size_t j1_coarse = (j1_fine + 1) / 2;
        size_t k1_coarse = (k1_fine + 1) / 2;
        size_t i2_coarse = (i2_fine + 1) / 2;
        size_t j2_coarse = (j2_fine + 1) / 2;
        size_t k2_coarse = (k2_fine + 1) / 2;

#ifndef BENCHMARKING
        if (i2_fine - i1_fine + 1 < domain_data->get_nx(level - 1)
            && i2_coarse - i1_coarse + 1 >= domain_data->get_nx(level)) {
            m_logger->warn("Be cautious! Obstacle '{}' fills up inner cells in x-direction at level {}", obstacle_fine->get_name(), level);
        }
        if (j2_fine - j1_fine + 1 < domain_data->get_ny(level - 1)
            && j2_coarse - j1_coarse + 1 >= domain_data->get_ny(level)) {
            m_logger->warn("Be cautious! Obstacle '{}' fills up inner cells in y-direction at level {}", obstacle_fine->get_name(), level);
        }
        if (k2_fine - k1_fine + 1 < domain_data->get_nz(level - 1)
            && k2_coarse - k1_coarse + 1 >= domain_data->get_nz(level)) {
            m_logger->warn("Be cautious! Obstacle '{}' fills up inner cells in z-direction at level {}", obstacle_fine->get_name(), level);
        }

        for (size_t c = 0; c < id; c++) {
            if (obstacle_list_coarse[c]->has_overlap(i1_coarse, i2_coarse, j1_coarse, j2_coarse, k1_coarse, k2_coarse)) {
                m_logger->debug("overlapping of obstacle {} with obstacle {} on level {}",
                                obstacle_list_coarse[c]->get_name(), obstacle_fine->get_name(), level);
            }
        }
#endif

        auto obstacle_coarse = new Obstacle(i1_coarse, j1_coarse, k1_coarse,
                                            i2_coarse, j2_coarse, k2_coarse,
                                            level, obstacle_fine->get_name());
#ifndef BENCHMARKING
        if (obstacle_coarse->get_stride_x() <= 1) {
            m_logger->warn("Obstacle '{}' is too small with size 1 in x-direction at level {}. "
                           "Consider less multigrid level, a higher resolution at the finest grid "
                           "or expanding the obstacle. Otherwise only the right boundary condition "
                           "will be applied.", obstacle_fine->get_name(), level);
        }
        if (obstacle_coarse->get_stride_y() <= 1) {
            m_logger->warn("Obstacle '{}' is too small with size 1 in y-direction at level {}. "
                           "Consider less multigrid level, a higher resolution at the finest grid "
                           "or expanding the obstacle. Otherwise only the top boundary condition "
                           "will be applied.", obstacle_fine->get_name(), level);
        }
        if (obstacle_coarse->get_stride_z() <= 1) {
            m_logger->warn("Obstacle '{}' is too small with size 1 in z-direction at level {}. "
                           "Consider less multigrid level, a higher resolution at the finest grid "
                           "or expanding the obstacle. Otherwise only the back boundary condition "
                           "will be applied.", obstacle_fine->get_name(), level);
        }
#endif
        obstacle_list_coarse[id] = obstacle_coarse;
        *sum_patches += obstacle_coarse->get_size_boundary_list();
    }

    size_t *list = obstacle_list_coarse[0]->get_obstacle_list();
    size_t size = obstacle_list_coarse[0]->get_size_obstacle_list();
    std::vector<size_t> data;
    data.assign(list, list + size);

    for (size_t o = 1; o < m_number_of_obstacle_objects; o++) {
        Obstacle *obstacle = obstacle_list_coarse[o];
        data = Utility::mergeSortedListsToUniqueList(list, size, obstacle->get_obstacle_list(),
                                                     obstacle->get_size_obstacle_list());
        list = data.data();
        size = data.size();
    }

    auto obstacle_list_tmp = new size_t[size];
    std::copy(&data[0], &data[size], obstacle_list_tmp);
    tmp_store_obstacle[level] = obstacle_list_tmp;
    return size;
}

// ================================= Send lists to GPU =============================================
// *************************************************************************************************
/// \brief  create joined list and send them to GPU
// *************************************************************************************************
void Multigrid::send_lists_to_GPU() {
    send_surface_lists_to_GPU();
    create_domain_lists_for_GPU();
    send_obstacle_lists_to_GPU();
}

// ================================= Send boundary lists to GPU ====================================
// *************************************************************************************************
/// \brief  create boundary joined list and send them to GPU
// *************************************************************************************************
void Multigrid::create_domain_lists_for_GPU() {
    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        Domain *domain = m_MG_domain_object_list[level];
        m_jl_domain_inner_list.add_data(level, domain->get_size_inner_list(), domain->get_inner_list());

        size_t **all_boundaries = domain->get_boundary_list();
        PatchObject all_boundary_sizes = domain->get_size_boundary_list();
        for (size_t patch = 0; patch < number_of_patches; patch++) {
            m_jl_domain_boundary_list_patch_divided[patch]->add_data(level, all_boundary_sizes[patch], all_boundaries[patch]);
        }
    }

    m_jl_domain_inner_list.copyin();
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_jl_domain_boundary_list_patch_divided[patch]->copyin();
    }
}

// ================================= Send surface list to GPU ======================================
// *************************************************************************************************
/// \brief  create surface joined list and send it to GPU
// *************************************************************************************************
void Multigrid::send_surface_lists_to_GPU() {
    size_t counter_surface_list = 0;

    size_t size_surface_list = 0;
    if (m_number_of_surface_objects > 0) {
        size_surface_list = get_length_of_surface_index_list_joined();
        m_data_MG_surface_list_level_joined = new size_t[size_surface_list];
#ifndef BENCHMARKING
        m_logger->debug("control send_surface_lists_to_GPU size surface {}", size_surface_list);
#endif
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            Surface **surface_list = *(m_MG_surface_object_list + level);
            for (size_t s = 0; s < m_number_of_surface_objects; s++) {
                Surface *surface = *(surface_list + s);
                for (size_t i = 0; i < surface->getSize_surfaceList(); i++) {
                    *(m_data_MG_surface_list_level_joined + counter_surface_list) = surface->getSurfaceList()[i];
                    counter_surface_list++;
                }
            }
        }
#ifndef BENCHMARKING
        m_logger->debug("control sendMGListsToGPU surface {} | {}",
                        counter_surface_list, size_surface_list);
#endif
#ifdef GPU_DEBUG
    m_gpu_logger->info("copyin surface list joined (m_data_MG_surface_list_level_joined) with pointer {} and size {}", static_cast<void *>(m_data_MG_surface_list_level_joined), size_surface_list);
#endif
#pragma acc enter data copyin(m_data_MG_surface_list_level_joined[:size_surface_list])
    }
}

// ================================= Send obstacle lists to GPU ====================================
// *************************************************************************************************
/// \brief  create obstacle joined list and send them to GPU
// *************************************************************************************************
void Multigrid::send_obstacle_lists_to_GPU() {
    if (m_number_of_obstacle_objects > 0) {
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            Obstacle **obstacle_list = m_MG_obstacle_object_list[level];
            for (size_t id = 0; id < m_number_of_obstacle_objects; id++) {
                Obstacle *obstacle = obstacle_list[id];
                size_t **all_obstacle_boundary_cells = obstacle->get_boundary_list();
                PatchObject all_obstacle_size = obstacle->get_size_boundary_list();
                for (size_t patch = 0; patch < number_of_patches; patch++) {
                    m_jl_obstacle_boundary_list_patch_divided[patch]->add_data(level, id, all_obstacle_size[patch], all_obstacle_boundary_cells[patch]);
                }
            }
        }
        m_jl_obstacle_list.copyin();
        for (size_t patch = 0; patch < number_of_patches; patch++) {
            m_jl_obstacle_boundary_list_patch_divided[patch]->copyin();
        }
    }
}

// ================================= Apply boundary condition ======================================
// *************************************************************************************************
/// \brief  Apply boundary condition for obstacles, surfaces and domain
/// \param  field Field
/// \param  sync synchronous kernel launching (true, default: false)
// *************************************************************************************************
void Multigrid::apply_boundary_condition(Field &field, bool sync) {
    real *d = field.data;
    size_t level = field.get_level();
    FieldType f = field.get_type();
    if (m_number_of_surface_objects > 0) {
        Surface **surface_list = *(m_MG_surface_object_list + level);
        for (size_t id = 0; id < m_number_of_surface_objects; ++id) {
            (static_cast<Surface *>(*(surface_list + id)))->applyBoundaryConditions(d, f, level, sync);
        }
    }
    if (m_number_of_obstacle_objects > 0) {
        for (size_t id = 0; id < m_number_of_obstacle_objects; ++id) {
            (static_cast<BoundaryDataController *> (m_bdc_obstacle[id]))->apply_boundary_condition_obstacle(field, m_jl_obstacle_boundary_list_patch_divided, f, id, sync);
        }
    }
    m_bdc_boundary->apply_boundary_condition(field, m_jl_domain_boundary_list_patch_divided, sync);
}

//======================================== Update lists ====================================
// ***************************************************************************************
/// \brief  Updates lists of indices
// ***************************************************************************************
void Multigrid::update_lists() {
    //TODO(issue 178) rework update function especially with obstacles
    remove_boundary_lists_from_GPU();

    if (m_number_of_obstacle_objects > 0) {
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            Domain *domain = *(m_MG_domain_object_list + level);
            domain->update_lists(*(m_MG_obstacle_object_list + level),
                                 m_number_of_obstacle_objects,
                                 m_jl_obstacle_list.get_slice_size(level));

            //m_size_MG_inner_index_list_level[level + 1] = m_size_MG_inner_index_list_level[level]
            //                                              + domain->get_size_inner_list();
            //m_size_MG_boundary_index_list_level[level + 1] = m_size_MG_boundary_index_list_level[level]
            //                                                 + domain->get_size_domain_list();
            //m_size_MG_boundary_slice_z_level[level + 1] = m_size_MG_boundary_slice_z_level[level]
            //                                              + domain->get_size_boundary_list()[FRONT];
            //m_size_MG_boundary_slice_y_level[level + 1] = m_size_MG_boundary_slice_y_level[level]
            //                                              + domain->get_size_boundary_list()[TOP];
            //m_size_MG_boundary_slice_x_level[level + 1] = m_size_MG_boundary_slice_x_level[level]
            //                                              + domain->get_size_boundary_list()[LEFT];
        }
    } else {
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            Domain *domain = m_MG_domain_object_list[level];
            domain->update_lists();
            //m_size_MG_inner_index_list_level[level + 1] = m_size_MG_inner_index_list_level[level]
            //                                              + domain->get_size_inner_list();
            //m_size_MG_boundary_index_list_level[level + 1] = m_size_MG_boundary_index_list_level[level]
            //                                                 + domain->get_size_domain_list();
            //m_size_MG_boundary_slice_z_level[level + 1] = m_size_MG_boundary_slice_z_level[level]
            //                                              + domain->get_size_boundary_list()[FRONT];
            //m_size_MG_boundary_slice_y_level[level + 1] = m_size_MG_boundary_slice_y_level[level]
            //                                              + domain->get_size_boundary_list()[TOP];
            //m_size_MG_boundary_slice_x_level[level + 1] = m_size_MG_boundary_slice_x_level[level]
            //                                              + domain->get_size_boundary_list()[LEFT];
        }
    }
    create_domain_lists_for_GPU();
}

void Multigrid::remove_boundary_lists_from_GPU() {
    //TODO(issue 130)
}

// surface_list
size_t Multigrid::get_length_of_surface_index_list_joined() const {
    size_t len = 0;
    if (m_number_of_surface_objects > 0) {
        len = get_first_index_of_surface_index_list(m_multigrid_levels + 1) + 1;
    }
    return len;
}

size_t Multigrid::get_first_index_of_surface_index_list(size_t level) const {
    size_t index = 0;
    if (m_number_of_surface_objects > 0) {
        index = *(m_size_MG_surface_index_list_level + level * m_number_of_surface_objects);
    }
    return index;
}

size_t Multigrid::get_last_index_of_surface_index_list(size_t level) const {
    size_t index = 0;
    if (m_number_of_surface_objects > 0) {
        index = *(m_size_MG_surface_index_list_level + (level + 1) * m_number_of_surface_objects) - 1;
    }
    return index;
}

bool Multigrid::is_obstacle_cell(const size_t level, const size_t index) {
    Obstacle **obstacle_list = m_MG_obstacle_object_list[level];
    for (size_t id = 0; id < m_number_of_obstacle_objects; id++) {
        Obstacle *obstacle = obstacle_list[id];
        if (obstacle->is_obstacle_cell(index)) {
            return true;
        }
    }
    return false;
}

bool Multigrid::is_obstacle_cell(const size_t level,
                                 const size_t i, const size_t j, const size_t k) {
    const size_t Nx = DomainData::getInstance()->get_Nx(level);
    const size_t Ny = DomainData::getInstance()->get_Ny(level);
    return is_obstacle_cell(level, IX(i, j, k, Nx, Ny));
}
