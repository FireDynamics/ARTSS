/// \file       Multigrid.cpp
/// \brief      Creates all lists needed for multigrid
/// \date       Oct 01, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Multigrid.h"
#include "PatchObject.h"
#include "../utility/Algorithm.h"
#include <string>
#include <vector>

Multigrid::Multigrid(size_t number_of_surfaces, Surface **surface_list,
                     size_t number_of_obstacles, Obstacle **obstacle_list,
                     BoundaryDataController *bdc_boundary,
                     BoundaryDataController **bdc_obstacles,
                     BoundaryDataController **bdc_surfaces,
                     size_t multigrid_levels) :
        m_multigrid_levels(multigrid_levels),
        m_number_of_surface_objects(number_of_surfaces),
        m_number_of_obstacle_objects(number_of_obstacles),
        m_jl_domain_list(multigrid_levels),
        m_jl_domain_inner_list(multigrid_levels),
        m_jl_obstacle_list(multigrid_levels),
        m_bdc_boundary(bdc_boundary),
        m_bdc_obstacle(bdc_obstacles),
        m_bdc_surface(bdc_surfaces) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
    m_logger->debug("starting multigrid");
#endif
#ifdef GPU_DEBUG
    m_gpu_logger = Utility::create_gpu_logger(typeid(this).name());
#endif
    // init domain
    // list of domain objects for each level
    m_MG_domain_object_list = new Domain *[m_multigrid_levels + 1];
    m_jl_domain_boundary_list_patch_divided = new SingleJoinedList*[number_of_patches];
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_jl_domain_boundary_list_patch_divided[patch]  = new SingleJoinedList(m_multigrid_levels);
    }

    // init obstacle
    // list of obstacle objects for each level
    m_MG_obstacle_object_list = new Obstacle **[m_multigrid_levels + 1];
    m_MG_obstacle_object_list[0] = obstacle_list;  // level 0

    m_jl_obstacle_boundary_list_patch_divided = new MultipleJoinedList*[number_of_patches];
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_jl_obstacle_boundary_list_patch_divided[patch]  = new MultipleJoinedList(m_multigrid_levels, m_number_of_obstacle_objects);
    }

    m_MG_surface_object_list = new Surface **[m_multigrid_levels + 1];
    m_MG_surface_object_list[0] = surface_list;  // level 0

    m_jl_surface_list_patch_divided = new MultipleJoinedList*[number_of_patches];
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_jl_surface_list_patch_divided[patch] = new MultipleJoinedList(m_multigrid_levels, m_number_of_surface_objects);
    }

    if (m_number_of_obstacle_objects > 0) {
        create_multigrid_obstacle_lists();
        send_obstacle_lists_to_GPU();
    }
    if (m_number_of_surface_objects > 0) {
        create_multigrid_surface_lists();
        send_surface_lists_to_GPU();
    }
    create_multigrid_domain_lists();

    send_domain_lists_to_GPU();
#ifndef BENCHMARKING
    print();
    control();
    m_logger->debug("end multigrid");
#endif
}

Multigrid::~Multigrid() {
    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        Surface **surface_level = m_MG_surface_object_list[level];
        for (size_t surface = 0; surface < m_number_of_surface_objects; surface++) {
            delete surface_level[surface];
        }
        delete[] surface_level;
        Obstacle **obstacle_level = m_MG_obstacle_object_list[level];
        for (size_t obstacle = 0; obstacle < m_number_of_obstacle_objects; obstacle++) {
            delete obstacle_level[obstacle];
        }
        delete[] obstacle_level;
        delete m_MG_domain_object_list[level];
    }
    delete[] m_MG_surface_object_list;
    delete[] m_MG_domain_object_list;
    delete[] m_MG_obstacle_object_list;


    for (size_t patch = 0; patch < number_of_patches; patch++) {
        delete m_jl_surface_list_patch_divided[patch];
        delete m_jl_obstacle_boundary_list_patch_divided[patch];
        delete m_jl_domain_boundary_list_patch_divided[patch];
    }

    delete m_jl_domain_boundary_list_patch_divided;
    delete m_jl_obstacle_boundary_list_patch_divided;
    delete m_jl_surface_list_patch_divided;

    for (size_t id = 0; id < m_number_of_obstacle_objects; id++) {
        delete m_bdc_obstacle[id];
    }
    delete m_bdc_boundary;
    delete m_bdc_obstacle;
}

//======================================== Control =================================================
// *************************************************************************************************
/// \brief  Units test emergency solution
// *************************************************************************************************
void Multigrid::control() {
#ifndef BENCHMARKING
    std::string control_message;
    auto domain_data = DomainData::getInstance();
    for (size_t level = 0; level < m_multigrid_levels; level++) {
        std::string message;
        size_t original_len = (static_cast<Domain *>(m_MG_domain_object_list[level]))->get_size_inner_list();
        size_t calculated_len = m_jl_domain_inner_list.get_last_index(level) - m_jl_domain_inner_list.get_first_index(level) + 1;
        size_t saved_len = m_jl_domain_inner_list.get_slice_size(level);
        if (calculated_len != original_len) {
            size_t control = domain_data->get_nx(level) * domain_data->get_ny(level) * domain_data->get_nz(level);
            message += fmt::format("length calculated/stored of domain inner cells does not equals its original size size: "
                                   "original: {} saved: {} calculated: {} control: {}\n",
                                   original_len, saved_len, calculated_len, control);
        }
        original_len = (static_cast<Domain *>(m_MG_domain_object_list[level]))->get_size_domain_list();
        calculated_len = m_jl_domain_list.get_last_index(level) - m_jl_domain_list.get_first_index(level) + 1;
        saved_len = m_jl_domain_list.get_slice_size(level);
        if (calculated_len != original_len) {
            size_t control = (domain_data->get_Nx(level) * domain_data->get_Ny(level) * 2)
                             + (domain_data->get_Nx(level) * (domain_data->get_Nz(level) - 2)) * 2
                             + ((domain_data->get_Ny(level) - 2) * (domain_data->get_Nz(level) - 2)) * 2;
            message += fmt::format("length calculated/stored of domain boundary cells does not equals its original size size: "
                                   "original: {} saved: {} calculated: {} control: {}\n",
                                   original_len, saved_len, calculated_len, control);
        }
        PatchObject *boundary_size = (static_cast<Domain *>(m_MG_domain_object_list[level]))->get_size_boundary_list();
        for (size_t patch = 0; patch < number_of_patches; patch++) {
            original_len = (*boundary_size)[patch];
            saved_len = m_jl_domain_boundary_list_patch_divided[patch]->get_slice_size(level);
            calculated_len = m_jl_domain_boundary_list_patch_divided[patch]->get_last_index(level) - m_jl_domain_boundary_list_patch_divided[patch]->get_first_index(level) + 1;
            if (calculated_len != original_len || saved_len != original_len) {
                size_t control;
                if (patch == Patch::FRONT || patch == Patch::BACK) {
                    control = domain_data->get_Nx(level) * domain_data->get_Ny(level);
                } else if (patch == Patch::BOTTOM || patch == Patch::TOP) {
                    control = domain_data->get_Nx(level) * domain_data->get_Nz(level);
                } else if (patch == Patch::LEFT || patch == Patch::RIGHT) {
                    control = domain_data->get_Ny(level) * domain_data->get_Nz(level);
                }
                message += fmt::format("length calculated/stored of domain boundary '{}' does not equals its original size size: "
                                       "original: {} saved: {} calculated: {} control: {}\n",
                                       PatchObject::get_patch_name(static_cast<Patch>(patch)),
                                       original_len, saved_len, calculated_len, control);
            }
        }

        size_t csize_inner = domain_data->get_nx(level) * domain_data->get_ny(level) * domain_data->get_nz(level)
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
        if (!message.empty()) {
            control_message += "For Level " + std::to_string(level) + "\nsize control\n" + message;
        }
    }
    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        Domain *domain = *(m_MG_domain_object_list + level);
        auto *sum_surfaces = new PatchObject();
        for (size_t patch = 0; patch < number_of_patches; patch++) {
            sum_surfaces->add_value(patch, m_jl_surface_list_patch_divided[patch]->get_slice_size(level));
        }
        domain->control(m_jl_obstacle_list.get_slice_size(level), *sum_surfaces);
        delete sum_surfaces;
    }
    {
        size_t bsize_inner = 0;
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            bsize_inner += get_slice_size_domain_inner_cells_level_joined(level);
        }
        size_t csize_inner = get_size_domain_inner_cells_level_joined();
        if (bsize_inner != csize_inner) {
            control_message += "get_size_domain_inner_cells_level_joined does not equal the sum of each inner list "
                       + std::to_string(bsize_inner) + "|" + std::to_string(csize_inner) + "\n";
        }
    }

    if (!control_message.empty()) {
        control_message = "################ MULTIGRID CONTROL ################\n" + control_message
                  + "---------------- MULTIGRID CONTROL END ----------------";
        m_logger->warn(control_message);
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
        m_logger->debug("For Level {} domain list starts at index {} and ends with index {}", level,
                        m_jl_domain_list.get_first_index(level),
                        m_jl_domain_list.get_last_index(level));
        m_logger->debug("and the corresponding indices at this position: {} | {}",
                        m_jl_domain_list.get_data()[m_jl_domain_list.get_first_index(level)],
                        m_jl_domain_list.get_data()[m_jl_domain_list.get_last_index(level)]);
    }
    m_logger->debug("Total number of domain cells: {}", m_jl_domain_list.get_size());

    for (size_t patch = 0; patch < number_of_patches; patch++) {
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            size_t start_index = m_jl_domain_boundary_list_patch_divided[patch]->get_first_index(level);
            size_t end_index = m_jl_domain_boundary_list_patch_divided[patch]->get_last_index(level);
            m_logger->debug("For Level {} domain boundary '{}' starts at index {} and ends with index {}",
                            level, PatchObject::get_patch_name(patch),
                            start_index, end_index);
            m_logger->debug(" and the corresponding indices at this position {} | {}",
                            m_jl_domain_boundary_list_patch_divided[patch]->get_data()[start_index],
                            m_jl_domain_boundary_list_patch_divided[patch]->get_data()[end_index]);
        }
        m_logger->debug("Total number of '{}' boundary cells: {}",
                        PatchObject::get_patch_name(patch),
                        m_jl_domain_boundary_list_patch_divided[patch]->get_size());
    }

    if (m_number_of_obstacle_objects > 0) {
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            m_logger->debug("For Level {} obstacle_list has {} elements",
                            level, m_jl_obstacle_list.get_slice_size(level));
        }

        for (size_t patch = 0; patch < number_of_patches; patch++) {
            for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
                m_logger->debug("For Level {} obstacle '{}' starts at index {} and ends with index {} with length {}",
                                level,
                                PatchObject::get_patch_name(patch),
                                m_jl_obstacle_boundary_list_patch_divided[patch]->get_first_index(level, 0),
                                m_jl_obstacle_boundary_list_patch_divided[patch]->get_last_index(level, m_number_of_obstacle_objects - 1),
                                m_jl_obstacle_boundary_list_patch_divided[patch]->get_slice_size(level));
            }
            m_logger->debug("Total number of '{}' obstacle cells: {}",
                            PatchObject::get_patch_name(patch),
                            m_jl_obstacle_boundary_list_patch_divided[patch]->get_size());
        }
    }

    if (m_number_of_surface_objects > 0) {
        size_t sum_surface_cells = 0;
        for (size_t patch = 0; patch < number_of_patches; patch++) {
            for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
                m_logger->debug("For Level {} surface_list starts at index {} and ends with index {}",
                                level,
                                m_jl_surface_list_patch_divided[patch]->get_first_index(level, 0),
                                m_jl_surface_list_patch_divided[patch]->get_last_index(level, m_number_of_surface_objects - 1));
            }
            sum_surface_cells += m_jl_surface_list_patch_divided[patch]->get_size();
        }
        m_logger->debug("Total length of surface_list: {}", sum_surface_cells);
    }
    m_logger->debug("---------------- MULTIGRID END ----------------");
#endif
}


//======================================== Print ===================================================
// *************************************************************************************************
/// \brief  create obstacle for all multigrid level. Obstacle objects are stored in
/// m_MG_obstacle_object_list and the joined list for the GPU are stored in
/// m_jl_obstacle_list and m_jl_obstacle_list_patch_divided
// *************************************************************************************************
void Multigrid::create_multigrid_obstacle_lists() {
#ifndef BENCHMARKING
    m_logger->debug("create_multigrid_obstacle_list");
#endif
    size_t sum_obstacle_cells = 0;  // total sum of obstacle cells for allocation m_jl_obstacle_list
    auto sum_obstacle_boundary = new PatchObject();  // total sum of obstacle patches for allocation m_jl_obstacle_boundary_list_patch_divided
    auto sum_obstacle_cells_level_divided = new size_t[m_multigrid_levels + 1];  // sum for obstacle cells for each level (total sum equals sum_obstacle_cells)
    std::fill(sum_obstacle_cells_level_divided, sum_obstacle_cells_level_divided + m_multigrid_levels +1, 0);

    if (m_number_of_obstacle_objects > 0) {
        auto **tmp_store_obstacle = new size_t *[m_multigrid_levels + 1];
        {  // count obstacle cells level 0
            size_t level = 0;
            for (size_t id = 0; id < m_number_of_obstacle_objects; id++) {
                sum_obstacle_cells_level_divided[level] += m_MG_obstacle_object_list[level][id]->get_size_obstacle_list();
                PatchObject *obstacle_boundary = m_MG_obstacle_object_list[level][id]->get_size_boundary_list();
                *sum_obstacle_boundary += *obstacle_boundary;
                sum_obstacle_cells += sum_obstacle_cells_level_divided[level];
            }

            // merge obstacles
            size_t *list = m_MG_obstacle_object_list[level][0]->get_obstacle_list();
            size_t size = m_MG_obstacle_object_list[level][0]->get_size_obstacle_list();
            for (size_t o = 1; o < m_number_of_obstacle_objects; o++) {
                //TODO (cvm) better idea for merging all obstacles ?
                Obstacle *obstacle = m_MG_obstacle_object_list[level][o];
                size_t new_size = size + obstacle->get_size_obstacle_list();
                auto obstacle_list_tmp = new size_t[new_size];
                Algorithm::merge_sort(list, obstacle->get_obstacle_list(),
                                      size, obstacle->get_size_obstacle_list(),
                                      obstacle_list_tmp);
                delete[] list;
                list = new size_t[new_size];
                std::copy(obstacle_list_tmp, obstacle_list_tmp + new_size, list);
                size = new_size;
                delete[] obstacle_list_tmp;
            }
            tmp_store_obstacle[level] = list;
#ifndef BENCHMARKING
            m_logger->debug("create_multigrid_obstacle_list ! calculated size {} original size {}", size, sum_obstacle_cells_level_divided[level]);
#endif
        }

        // create obstacles for each multigrid level
        for (size_t level = 1; level < m_multigrid_levels + 1; level++) {
            size_t sum = obstacle_dominant_restriction(level, sum_obstacle_boundary, tmp_store_obstacle);
            sum_obstacle_cells_level_divided[level] = sum;
            sum_obstacle_cells += sum_obstacle_cells_level_divided[level];
        }

        for (size_t patch = 0; patch < number_of_patches; patch++) {
            m_jl_obstacle_boundary_list_patch_divided[patch]->set_size((*sum_obstacle_boundary)[patch]);
        }
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            Obstacle **obstacle_list = m_MG_obstacle_object_list[level];
            for (size_t id = 0; id < m_number_of_obstacle_objects; id++) {
                Obstacle *obstacle = obstacle_list[id];
                size_t **obstacle_boundary_cells = obstacle->get_boundary_list();
                PatchObject *obstacle_size = obstacle->get_size_boundary_list();
                for (size_t patch = 0; patch < number_of_patches; patch++) {
#ifndef BENCHMARKING
                    m_logger->debug("add obstacle boundary data to MJL for obstacle '{}' patch '{}' level={}, number of obstacle boundary cells={}",
                                    obstacle->get_name(), PatchObject::get_patch_name(patch), level,
                                    (*obstacle_size)[patch]);
#endif
                    m_jl_obstacle_boundary_list_patch_divided[patch]->add_data(level, id, (*obstacle_size)[patch], obstacle_boundary_cells[patch]);
                }
            }
        }

        m_jl_obstacle_list.set_size(sum_obstacle_cells);
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
#ifndef BENCHMARKING
            m_logger->debug("add obstacle data to SJL, level={} number of obstacle cells={}", level,
                            sum_obstacle_cells_level_divided[level]);
#endif
            m_jl_obstacle_list.add_data(level, sum_obstacle_cells_level_divided[level], tmp_store_obstacle[level]);
        }

        delete sum_obstacle_boundary;
        delete[] sum_obstacle_cells_level_divided;
        for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
            delete[] tmp_store_obstacle[level];
        }
        delete[] tmp_store_obstacle;
    }
}


void Multigrid::create_multigrid_surface_lists() {
#ifndef BENCHMARKING
    m_logger->debug("create_multigrid_surface_list");
#endif
    size_t sum_surface_cells = 0;
    auto sum_surface_patch_divided = new PatchObject();  // total sum of surfaces for allocation m_jl_surface_list_patch_divided
    {  // count surface cells level 0
        size_t level = 0;
        for (size_t id = 0; id < m_number_of_surface_objects; id++) {
            Patch patch = m_MG_surface_object_list[level][id]->get_patch();
            sum_surface_patch_divided->add_value(patch, m_MG_surface_object_list[level][id]->get_size_surface_list());
        }
        sum_surface_cells = sum_surface_patch_divided->get_sum();
    }

    // create surfaces for each multigrid level
    for (size_t level = 1; level < m_multigrid_levels + 1; level++) {
        sum_surface_cells += surface_dominant_restriction(level, sum_surface_patch_divided);
    }

    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_jl_surface_list_patch_divided[patch]->set_size((*sum_surface_patch_divided)[patch]);
    }
    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        for (size_t id = 0; id < m_number_of_surface_objects; level++) {
            Patch patch = m_MG_surface_object_list[level][id]->get_patch();
#ifndef BENCHMARKING
            m_logger->debug("add surface data ({}) to MJL with patch={} level={} pointer={}",
                            m_MG_surface_object_list[level][id]->get_name(),
                            PatchObject::get_patch_name(patch), level,
                            m_MG_surface_object_list[level][id]->get_size_surface_list());
#endif
            // all other patches have size 0
            m_jl_surface_list_patch_divided[patch]->add_data(
                    level, id,
                    m_MG_surface_object_list[level][id]->get_size_surface_list(),
                    m_MG_surface_object_list[level][id]->get_surface_list());
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("surface sums: {}|{}", sum_surface_cells, sum_surface_patch_divided->get_sum());
#endif
    delete sum_surface_patch_divided;
}

void Multigrid::create_multigrid_domain_lists() {
#ifndef BENCHMARKING
    m_logger->debug("create_multigrid_domain_list");
#endif
    size_t sum_domain_inner_cells = 0;  // total sum of domain inner cells for allocating m_jl_domain_inner_list
    size_t sum_domain_cells = 0;  // total sum of domain cells for allocating m_jl_domain_list
    auto sum_domain_boundary = new PatchObject();  // total sum of domain patches for allocation m_jl_boundary_list_patch_divided

    // create domain for each multigrid level
    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        auto slice_surface_list = new size_t*[number_of_patches];
        auto size_surface_list = new PatchObject();
        for (size_t patch = 0; patch < number_of_patches; patch++) {
            slice_surface_list[patch] = m_jl_surface_list_patch_divided[patch]->get_slice(level);
            size_surface_list->add_value(patch, m_jl_surface_list_patch_divided[patch]->get_slice_size(level));
        }
        size_t *slice_obstacle_list = m_jl_obstacle_list.get_slice(level);
        // TODO(issue 178) send only obstacles which are in the computational domain
        auto domain = new Domain(slice_obstacle_list,
                                 m_jl_obstacle_list.get_slice_size(level),
                                 slice_surface_list,
                                 *size_surface_list,
                                 level);
        // save domain object
        m_MG_domain_object_list[level] = domain;
        // count domain cells
        sum_domain_inner_cells += m_MG_domain_object_list[level]->get_size_inner_list();
        sum_domain_cells += m_MG_domain_object_list[level]->get_size_domain_list();
        PatchObject *domain_boundary_size = m_MG_domain_object_list[level]->get_size_boundary_list();
        *sum_domain_boundary += *domain_boundary_size;
        for (size_t patch = 0; patch < number_of_patches; patch++) {
            delete[] slice_surface_list[patch];
        }
        delete[] slice_surface_list;
        delete[] slice_obstacle_list;
        delete size_surface_list;
    }

    m_jl_domain_inner_list.set_size(sum_domain_inner_cells);
    m_jl_domain_list.set_size(sum_domain_cells);

    for (size_t patch = 0 ; patch < number_of_patches; patch++) {
        m_jl_domain_boundary_list_patch_divided[patch]->set_size((*sum_domain_boundary)[patch]);
    }

    for (size_t level = 0; level < m_multigrid_levels + 1; level++) {
        Domain *domain = m_MG_domain_object_list[level];
        m_jl_domain_inner_list.add_data(level, domain->get_size_inner_list(), domain->get_inner_list());

        m_jl_domain_list.add_data(level, domain->get_size_domain_list(), domain->get_domain_list());

        size_t **all_boundaries = domain->get_boundary_list();
        PatchObject *all_boundary_sizes = domain->get_size_boundary_list();
        for (size_t patch = 0; patch < number_of_patches; patch++) {
#ifndef BENCHMARKING
            m_logger->debug("add domain boundary data to SJL for patch '{}' level {} size {}",
                            PatchObject::get_patch_name(patch), level,
                            (*all_boundary_sizes)[level]);
#endif
            m_jl_domain_boundary_list_patch_divided[patch]->add_data(level, (*all_boundary_sizes)[patch], all_boundaries[patch]);
        }
    }
    delete sum_domain_boundary;
}

// ================================= Surface dominant restriction ==================================
// *************************************************************************************************
/// \brief  Calculates coarse level of multigrid of surfaces (dominant restriction)
/// \param level Multigrid level
// *************************************************************************************************
size_t Multigrid::surface_dominant_restriction(size_t level, PatchObject *sum_patches) {
    //TODO(issue 5): warning if two surfaces are merging which have different BC
    // - add merging
// SURFACES
    size_t sum = 0;
    // add index to m_surface_list if any of l-1 indices building the l index was a surface
    // dominant restriction
    DomainData *domain_data = DomainData::getInstance();

    Surface **surface_list_fine = *(m_MG_surface_object_list + (level - 1));
    auto surface_list_coarse = new Surface *[m_number_of_surface_objects];
    *(m_MG_surface_object_list + level) = surface_list_coarse;
    // loop through surfaces
    for (size_t id = 0; id < m_number_of_surface_objects; id++) {
        Surface *surface_fine = surface_list_fine[id];

        Patch patch = surface_fine->get_patch();

        auto start_coarse = new Coordinate<size_t>();
        auto *end_coarse = new Coordinate<size_t>();
        for (size_t axis = 0; axis < number_of_axes; axis++) {
            size_t start_fine = surface_fine->get_start_index(CoordinateAxis(axis));
            (*start_coarse)[axis] = static_cast<size_t>((start_fine + 1) / 2);

            size_t end_fine = surface_fine->get_end_index(CoordinateAxis(axis));
            (*end_coarse)[axis] = static_cast<size_t>((end_fine + 1) / 2);
        }

#ifndef BENCHMARKING
        for (size_t axis = 0; axis < number_of_axes; axis++) {
            size_t start_fine = surface_fine->get_start_index(CoordinateAxis(axis));
            size_t end_fine = surface_fine->get_end_index(CoordinateAxis(axis));
            if (end_fine - start_fine + 1 < domain_data->get_nx(level - 1)
                && (*end_coarse)[axis] - (*start_coarse)[axis] + 1 >= domain_data->get_number_of_inner_cells(CoordinateAxis(axis), level)) {
                m_logger->warn("Be cautious! Surface '{}' fills up boundary patch '{}' at level {}", surface_fine->get_name(), level);
            }
        }
        m_logger->debug("multigrid surface start fine {}|{}|{}",
                        surface_fine->get_start_index(X),
                        surface_fine->get_start_index(Y),
                        surface_fine->get_start_index(Z));
        m_logger->debug("multigrid surface end fine {}|{}|{}",
                        surface_fine->get_end_index(X),
                        surface_fine->get_end_index(Y),
                        surface_fine->get_end_index(Z));
        m_logger->debug("multigrid surface start coarse {}|{}|{}",
                        (*start_coarse)[X],
                        (*start_coarse)[Y],
                        (*start_coarse)[Z]);
        m_logger->debug("multigrid surface end coarse {}|{}|{}",
                        (*end_coarse)[X],
                        (*end_coarse)[Y],
                        (*end_coarse)[Z]);
#endif

        auto surface_coarse = new Surface(*start_coarse, *end_coarse, level, surface_fine->get_name(), patch);
        sum += surface_coarse->get_size_surface_list();
        *(surface_list_coarse + id) = surface_coarse;
        sum_patches->add_value(patch, surface_coarse->get_size_surface_list());
    }  // end surface id loop
    return sum;
}

// ================================= Obstacle dominant restriction =================================
// *************************************************************************************************
/// \brief  Calculates coarse level of multigrid of obstacles (dominant restriction)
///         !!!Be careful with obstacle on coarse level!!!
///         Obstacles on coarser level may overlap due to the dominant restriction, therefore
///         the obstacles are merged at these level
/// \param level Multigrid level
/// \param sum_patches stores the patch size of their respective patch for all obstacles of level
/// \param tmp_store_obstacle returns the merged obstacle list
// *************************************************************************************************
size_t Multigrid::obstacle_dominant_restriction(size_t level, PatchObject *sum_patches, size_t **tmp_store_obstacle) {
    // OBSTACLES
    // add index to m_obstacle_lists if any of l-1 indices building the l index was an obstacle
    // (dominant restriction)
    if (m_number_of_obstacle_objects == 0) {
        tmp_store_obstacle[level] = new size_t[0];
        return 0;
    }
    DomainData *domain_data = DomainData::getInstance();
    Obstacle **obstacle_list_fine = m_MG_obstacle_object_list[level - 1];
    auto obstacle_list_coarse = new Obstacle *[m_number_of_obstacle_objects];
    m_MG_obstacle_object_list[level] = obstacle_list_coarse;
    for (size_t id = 0; id < m_number_of_obstacle_objects; id++) {
        Obstacle *obstacle_fine = obstacle_list_fine[id];

        Coordinate<size_t> *start_fine = obstacle_fine->get_start_coordinates();
        Coordinate<size_t> *end_fine = obstacle_fine->get_end_coordinates();

        auto start_coarse = new Coordinate<size_t>();
        auto end_coarse = new Coordinate<size_t>();
        for (size_t axis = 0; axis < number_of_axes; axis++) {
            (*start_coarse)[axis] = static_cast<size_t>(((*start_fine)[axis] + 1) / 2);
            (*end_coarse)[axis] = static_cast<size_t>(((*end_fine)[axis] + 1) / 2);
        }

#ifndef BENCHMARKING
        if ((*end_fine)[CoordinateAxis::X] - (*start_fine)[CoordinateAxis::X] + 1 < domain_data->get_nx(level - 1)
            && (*end_coarse)[CoordinateAxis::X] - (*start_coarse)[CoordinateAxis::X] + 1 >= domain_data->get_nx(level)) {
            m_logger->warn("Be cautious! Obstacle '{}' fills up inner cells in x-direction at level {}", obstacle_fine->get_name(), level);
        }
        if ((*end_fine)[CoordinateAxis::Y] - (*start_fine)[CoordinateAxis::Y] + 1 < domain_data->get_ny(level - 1)
            && (*end_coarse)[CoordinateAxis::Y] - (*start_coarse)[CoordinateAxis::Y] + 1 >= domain_data->get_ny(level)) {
            m_logger->warn("Be cautious! Obstacle '{}' fills up inner cells in y-direction at level {}", obstacle_fine->get_name(), level);
        }
        if ((*end_fine)[CoordinateAxis::Z] - (*start_fine)[CoordinateAxis::Z] + 1 < domain_data->get_nz(level - 1)
            && (*end_coarse)[CoordinateAxis::Z] - (*start_coarse)[CoordinateAxis::Z] + 1 >= domain_data->get_nz(level)) {
            m_logger->warn("Be cautious! Obstacle '{}' fills up inner cells in z-direction at level {}", obstacle_fine->get_name(), level);
        }

        for (size_t c = 0; c < id; c++) {
            if (obstacle_list_coarse[c]->has_overlap((*start_coarse)[CoordinateAxis::X], (*end_coarse)[CoordinateAxis::X], (*start_coarse)[CoordinateAxis::Y], (*end_coarse)[CoordinateAxis::Y], (*start_coarse)[CoordinateAxis::Z], (*end_coarse)[CoordinateAxis::Z])) {
                m_logger->debug("overlapping of obstacle {} with obstacle {} on level {}",
                                obstacle_list_coarse[c]->get_name(), obstacle_fine->get_name(), level);
            }
        }
#endif

        auto obstacle_coarse = new Obstacle((*start_coarse), (*end_coarse),
                                            level, obstacle_fine->get_name());
        delete start_coarse;
        delete end_coarse;
#ifndef BENCHMARKING
        Coordinate<size_t> *strides = obstacle_coarse->get_strides();
        if ((*strides)[CoordinateAxis::X] <= 1) {
            m_logger->warn("Obstacle '{}' is too small with size 1 in x-direction at level {}. "
                           "Consider less multigrid level, a higher resolution at the finest grid "
                           "or expanding the obstacle. Otherwise only the right boundary condition "
                           "will be applied.", obstacle_fine->get_name(), level);
        }
        if ((*strides)[CoordinateAxis::Y] <= 1) {
            m_logger->warn("Obstacle '{}' is too small with size 1 in y-direction at level {}. "
                           "Consider less multigrid level, a higher resolution at the finest grid "
                           "or expanding the obstacle. Otherwise only the top boundary condition "
                           "will be applied.", obstacle_fine->get_name(), level);
        }
        if ((*strides)[CoordinateAxis::Z] <= 1) {
            m_logger->warn("Obstacle '{}' is too small with size 1 in z-direction at level {}. "
                           "Consider less multigrid level, a higher resolution at the finest grid "
                           "or expanding the obstacle. Otherwise only the back boundary condition "
                           "will be applied.", obstacle_fine->get_name(), level);
        }
#endif
        obstacle_list_coarse[id] = obstacle_coarse;
        *sum_patches += *obstacle_coarse->get_size_boundary_list();
    }

    // merge obstacles
    size_t *list = obstacle_list_coarse[0]->get_obstacle_list();
    size_t size = obstacle_list_coarse[0]->get_size_obstacle_list();
    std::vector<size_t> data;
    data.assign(list, list + size);

    for (size_t o = 1; o < m_number_of_obstacle_objects; o++) {
        Obstacle *obstacle = obstacle_list_coarse[o];
        data = Algorithm::merge_sort_with_duplicates(list, size, obstacle->get_obstacle_list(),
                                                     obstacle->get_size_obstacle_list());
        list = data.data();
        size = data.size();
    }

    //TODO(cvm) copy command
    auto obstacle_list_tmp = new size_t[size];
    std::copy(&data[0], &data[size], obstacle_list_tmp);
    tmp_store_obstacle[level] = obstacle_list_tmp;
    return size;
}

// ================================= Send boundary lists to GPU ====================================
// *************************************************************************************************
/// \brief  create boundary joined list and send them to GPU
// *************************************************************************************************
void Multigrid::send_domain_lists_to_GPU() {
    m_jl_domain_inner_list.copyin();
    m_jl_domain_list.copyin();
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_jl_domain_boundary_list_patch_divided[patch]->copyin();
    }
}

// ================================= Send surface list to GPU ======================================
// *************************************************************************************************
/// \brief  create surface joined list and send it to GPU
// *************************************************************************************************
void Multigrid::send_surface_lists_to_GPU() {
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_jl_surface_list_patch_divided[patch]->copyin();
    }
}

// ================================= Send obstacle lists to GPU ====================================
// *************************************************************************************************
/// \brief  send obstacle joined list to GPU
// *************************************************************************************************
void Multigrid::send_obstacle_lists_to_GPU() {
    m_jl_obstacle_list.copyin();
    for (size_t patch = 0; patch < number_of_patches; patch++) {
        m_jl_obstacle_boundary_list_patch_divided[patch]->copyin();
    }
}

// ================================= Apply boundary condition ======================================
// *************************************************************************************************
/// \brief  Apply boundary condition for obstacles, surfaces and domain
/// \param  field Field
/// \param  sync synchronous kernel launching (true, default: false)
// *************************************************************************************************
void Multigrid::apply_boundary_condition(Field &field, bool sync) {
    // order (shouldn't) matter. Only if there are still cells which have multiple
    // BC types (obstacle, boundary, surface), but these should be excluded
    // first obstacles BC, then domain BC and at last surface BC.
    if (m_number_of_obstacle_objects > 0) {
        for (size_t id = 0; id < m_number_of_obstacle_objects; ++id) {
            (static_cast<BoundaryDataController *> (m_bdc_obstacle[id]))->apply_boundary_condition_obstacle(field, m_jl_obstacle_boundary_list_patch_divided, id, sync);
        }
    }

    m_bdc_boundary->apply_boundary_condition(field, m_jl_domain_boundary_list_patch_divided, sync);

    if (m_number_of_surface_objects > 0) {
        for (size_t id = 0; id < m_number_of_surface_objects; ++id) {
            (static_cast<BoundaryDataController *> (m_bdc_surface[id]))->apply_boundary_condition_surface(field, m_jl_surface_list_patch_divided, id, sync);
        }
    }
}

//======================================== Update lists ====================================
// ***************************************************************************************
/// \brief  Updates lists of indices
// ***************************************************************************************
void Multigrid::update_lists() {
    //TODO(issue 178) rework update function especially with obstacles
    // develop concept where obstacles and surfaces are distinguished by whether
    // they are in the computational domain or not.
    // m_MG_domain_object_list[0]->update_lists()
    remove_domain_lists_from_GPU();
    create_multigrid_domain_lists();
    send_domain_lists_to_GPU();
}

void Multigrid::remove_domain_lists_from_GPU() {
    //TODO(issue 178) remove domain lists which are going to be replaced
    // probably the easiest way would be to implement a exit from gpu function
    // the joined lists
}

bool Multigrid::is_obstacle_cell(const size_t level,
                                 const Coordinate<size_t> &coords) {
    Obstacle **obstacle_list = m_MG_obstacle_object_list[level];
    for (size_t id = 0; id < m_number_of_obstacle_objects; id++) {
        Obstacle *obstacle = obstacle_list[id];
        if (obstacle->is_obstacle_cell(coords)) {
            return true;
        }
    }
    return false;
}
