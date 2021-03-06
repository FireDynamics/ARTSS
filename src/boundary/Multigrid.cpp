/// \file       Multigrid.h
/// \brief      Creates all lists needed for multigrid
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Multigrid.h"
#include <algorithm>
#include <string>
#include <vector>

Multigrid::Multigrid(BoundaryDataController *bdc_boundary) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_bdc_boundary = bdc_boundary;
    m_number_of_surface_objects = 0;
    m_number_of_obstacle_objects = 0;

    init();
    add_MG_lists();
    send_lists_to_GPU();
#ifndef BENCHMARKING
    print();
    control();
#endif
    m_data_boundary_patches_joined = new size_t *[number_of_patches];
    m_data_boundary_patches_joined[Patch::FRONT] = m_data_MG_boundary_front_level_joined;
    m_data_boundary_patches_joined[Patch::BACK] = m_data_MG_boundary_back_level_joined;
    m_data_boundary_patches_joined[Patch::BOTTOM] = m_data_MG_boundary_bottom_level_joined;
    m_data_boundary_patches_joined[Patch::TOP] = m_data_MG_boundary_top_level_joined;
    m_data_boundary_patches_joined[Patch::LEFT] = m_data_MG_boundary_left_level_joined;
    m_data_boundary_patches_joined[Patch::RIGHT] = m_data_MG_boundary_right_level_joined;
}

Multigrid::Multigrid(
        size_t number_of_surfaces, Surface **surface_list,
        size_t number_of_obstacles, Obstacle **obstacle_list,
        BoundaryDataController *bdc_boundary,
        BoundaryDataController **bdc_obstacles) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_bdc_boundary = bdc_boundary;
    m_number_of_surface_objects = number_of_surfaces;
    m_number_of_obstacle_objects = number_of_obstacles;

    m_levels = Domain::getInstance()->get_levels();  // multigrid level, 0 otherwise

    if (m_number_of_surface_objects > 0) {
        // list of surfaces for each level
        m_MG_surface_object_list = new Surface **[m_levels + 1];
        *(m_MG_surface_object_list) = surface_list;  // level 0

        // surface indices divided by level
        m_MG_surface_index_list = new size_t *[m_levels + 1];

        // start index of each surface in level joined list
        m_size_MG_surface_index_list_level = new size_t[m_levels * m_number_of_surface_objects + 1];
        // start index of first level in joined list = 0
        *(m_size_MG_surface_index_list_level) = 0;
    }

    if (m_number_of_obstacle_objects > 0) {
        // list of obstacles for each level
        m_MG_obstacle_object_list = new Obstacle **[m_levels + 1];
        *(m_MG_obstacle_object_list) = obstacle_list;  // level 0
        m_size_MG_obstacle_index_list_level = new size_t[m_levels + 1];

        // obstacle indices divided by level
        m_MG_obstacle_index_list = new size_t *[m_levels + 1];

        // start index of each obstacle in level joined list
        // (slice z = Front/Back, slice y = Bottom/Top, slice x = Left/Right)
        m_size_MG_obstacle_front_level = new size_t[(m_levels + 1) * m_number_of_obstacle_objects + 1];
        m_size_MG_obstacle_back_level = new size_t[(m_levels + 1) * m_number_of_obstacle_objects + 1];
        m_size_MG_obstacle_top_level = new size_t[(m_levels + 1) * m_number_of_obstacle_objects + 1];
        m_size_MG_obstacle_bottom_level = new size_t[(m_levels + 1) * m_number_of_obstacle_objects + 1];
        m_size_MG_obstacle_left_level = new size_t[(m_levels + 1) * m_number_of_obstacle_objects + 1];
        m_size_MG_obstacle_right_level = new size_t[(m_levels + 1) * m_number_of_obstacle_objects + 1];

        // size of obstacle level 0 / start index of level 1
        m_size_MG_obstacle_front_level[0] = 0;
        m_size_MG_obstacle_back_level[0] = 0;
        m_size_MG_obstacle_top_level[0] = 0;
        m_size_MG_obstacle_bottom_level[0] = 0;
        m_size_MG_obstacle_left_level[0] = 0;
        m_size_MG_obstacle_right_level[0] = 0;
        *(m_size_MG_obstacle_index_list_level) = 0;

        for (size_t i = 0; i < m_number_of_obstacle_objects; i++) {
            m_size_MG_obstacle_front_level[i + 1] =
                    obstacle_list[i]->get_size_obstacle_front() + m_size_MG_obstacle_front_level[i];
            m_size_MG_obstacle_back_level[i + 1] =
                    obstacle_list[i]->get_size_obstacle_back() + m_size_MG_obstacle_back_level[i];
            m_size_MG_obstacle_bottom_level[i + 1] =
                    obstacle_list[i]->get_size_obstacle_bottom() + m_size_MG_obstacle_bottom_level[i];
            m_size_MG_obstacle_top_level[i + 1] =
                    obstacle_list[i]->get_size_obstacle_top() + m_size_MG_obstacle_top_level[i];
            m_size_MG_obstacle_left_level[i + 1] =
                    obstacle_list[i]->get_size_obstacle_left() + m_size_MG_obstacle_left_level[i];
            m_size_MG_obstacle_right_level[i + 1] =
                    obstacle_list[i]->get_size_obstacle_right() + m_size_MG_obstacle_right_level[i];
            *(m_size_MG_obstacle_index_list_level) += obstacle_list[i]->get_size_obstacle_list();
        }
        m_bdc_obstacle = bdc_obstacles;
    }

    init();
    add_MG_lists();
    send_lists_to_GPU();
#ifndef BENCHMARKING
    print();
    control();
#endif
    m_data_boundary_patches_joined = new size_t *[number_of_patches];
    m_data_boundary_patches_joined[Patch::FRONT] = m_data_MG_boundary_front_level_joined;
    m_data_boundary_patches_joined[Patch::BACK] = m_data_MG_boundary_back_level_joined;
    m_data_boundary_patches_joined[Patch::BOTTOM] = m_data_MG_boundary_bottom_level_joined;
    m_data_boundary_patches_joined[Patch::TOP] = m_data_MG_boundary_top_level_joined;
    m_data_boundary_patches_joined[Patch::LEFT] = m_data_MG_boundary_left_level_joined;
    m_data_boundary_patches_joined[Patch::RIGHT] = m_data_MG_boundary_right_level_joined;

    if (m_number_of_obstacle_objects > 0) {
        m_data_obstacles_patches_joined = new size_t *[number_of_patches];
        m_data_obstacles_patches_joined[Patch::FRONT] = m_data_MG_obstacle_front_level_joined;
        m_data_obstacles_patches_joined[Patch::BACK] = m_data_MG_obstacle_back_level_joined;
        m_data_obstacles_patches_joined[Patch::BOTTOM] = m_data_MG_obstacle_bottom_level_joined;
        m_data_obstacles_patches_joined[Patch::TOP] = m_data_MG_obstacle_top_level_joined;
        m_data_obstacles_patches_joined[Patch::LEFT] = m_data_MG_obstacle_left_level_joined;
        m_data_obstacles_patches_joined[Patch::RIGHT] = m_data_MG_obstacle_right_level_joined;
    }
}

//======================================== Init ====================================================
// *************************************************************************************************
/// \brief  Initialize member variables (arrays)
// *************************************************************************************************
void Multigrid::init() {
    m_levels = Domain::getInstance()->get_levels();  // multigrid level, 0 otherwise

    // list of domain boundary for each level
    m_MG_boundary_object_list = new Boundary *[m_levels + 1];

    // start index of each level (difference equals size of respective element)
    m_size_MG_inner_index_list_level = new size_t[m_levels + 2];
    m_size_MG_boundary_index_list_level = new size_t[m_levels + 2];
    m_size_MG_boundary_slice_z_level = new size_t[m_levels + 2];
    m_size_MG_boundary_slice_y_level = new size_t[m_levels + 2];
    m_size_MG_boundary_slice_x_level = new size_t[m_levels + 2];

    // start index of first element = 0
    *(m_size_MG_inner_index_list_level) = 0;
    *(m_size_MG_boundary_index_list_level) = 0;

    *(m_size_MG_boundary_slice_z_level) = 0;
    *(m_size_MG_boundary_slice_y_level) = 0;
    *(m_size_MG_boundary_slice_x_level) = 0;
}

Multigrid::~Multigrid() {
    delete[] m_data_boundary_patches_joined;
    for (size_t level = 0; level < m_levels + 1; level++) {
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
        delete (*(m_MG_boundary_object_list + level));
    }
    delete[] m_MG_boundary_object_list;

    size_t size_inner_list = get_length_of_inner_index_list_joined();
    size_t size_boundary_list = get_length_of_boundary_index_list_joined();
#pragma acc exit data delete(m_data_MG_inner_list_level_joined[:size_inner_list])
#pragma acc exit data delete(m_data_MG_boundary_list_level_joined[:size_boundary_list])
    delete[] m_data_MG_inner_list_level_joined;
    delete[] m_data_MG_boundary_list_level_joined;

    size_t size_boundary_slice_z = get_length_of_boundary_slice_z_joined();
    size_t size_boundary_slice_y = get_length_of_boundary_slice_y_joined();
    size_t size_boundary_slice_x = get_length_of_boundary_slice_x_joined();
#pragma acc exit data delete(m_data_MG_boundary_front_level_joined[:size_boundary_slice_z])
#pragma acc exit data delete(m_data_MG_boundary_back_level_joined[:size_boundary_slice_z])
#pragma acc exit data delete(m_data_MG_boundary_top_level_joined[:size_boundary_slice_y])
#pragma acc exit data delete(m_data_MG_boundary_bottom_level_joined[:size_boundary_slice_y])
#pragma acc exit data delete(m_data_MG_boundary_left_level_joined[:size_boundary_slice_x])
#pragma acc exit data delete(m_data_MG_boundary_right_level_joined[:size_boundary_slice_x])
    delete[] m_data_MG_boundary_front_level_joined;
    delete[] m_data_MG_boundary_back_level_joined;
    delete[] m_data_MG_boundary_top_level_joined;
    delete[] m_data_MG_boundary_bottom_level_joined;
    delete[] m_data_MG_boundary_left_level_joined;
    delete[] m_data_MG_boundary_right_level_joined;

    delete[] m_size_MG_inner_index_list_level;
    delete[] m_size_MG_boundary_index_list_level;

    delete[] m_size_MG_boundary_slice_z_level;
    delete[] m_size_MG_boundary_slice_y_level;
    delete[] m_size_MG_boundary_slice_x_level;

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
        delete[] m_MG_obstacle_index_list;

        delete[] m_data_MG_obstacle_front_level_joined;
        delete[] m_data_MG_obstacle_back_level_joined;
        delete[] m_data_MG_obstacle_top_level_joined;
        delete[] m_data_MG_obstacle_bottom_level_joined;
        delete[] m_data_MG_obstacle_left_level_joined;
        delete[] m_data_MG_obstacle_right_level_joined;

        delete[] m_size_MG_obstacle_front_level;
        delete[] m_size_MG_obstacle_back_level;
        delete[] m_size_MG_obstacle_top_level;
        delete[] m_size_MG_obstacle_bottom_level;
        delete[] m_size_MG_obstacle_left_level;
        delete[] m_size_MG_obstacle_right_level;
    }
}

//======================================== Control =================================================
// *************************************************************************************************
/// \brief  Units test emergency solution
// *************************************************************************************************
void Multigrid::control() {
    std::string message;
    auto domain = Domain::getInstance();
    for (size_t level = 0; level < m_levels + 1; level++) {
        size_t nx = domain->get_nx(level);
        size_t ny = domain->get_ny(level);
        size_t nz = domain->get_nz(level);

        if (!message.empty()) {
            message += "For Level " + std::to_string(level) + "\n";
            message += "size control\n";
        }
        size_t bLen = (static_cast<Boundary *>(*(m_MG_boundary_object_list + level)))->get_size_inner_list();
        size_t cLen = get_last_index_of_inner_index_list(level) - get_first_index_of_inner_index_list(level) + 1;
        if (cLen != bLen) {
            size_t control = (domain->get_nx(level) - 2) * (domain->get_ny(level) - 2) * (domain->get_nz(level) - 2);
            message += "length calculated by first and last index of inner_list does not equals size of innerList of Boundary object "
                       + std::to_string(cLen) + "|" + std::to_string(bLen) + " control: " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundary_object_list + level)))->get_size_boundary_list();
        cLen = get_last_index_boundary_index_list(level) - get_first_index_of_boundary_index_list(level) + 1;
        if (cLen != bLen) {
            size_t control = (domain->get_nx(level) * domain->get_ny(level) * 2)
                             + (domain->get_nx(level) * (domain->get_nz(level) - 2)) * 2
                             + ((domain->get_ny(level) - 2) * (domain->get_nz(level) - 2)) * 2;
            message += "length calculated by first and last index of boundary_list does not equals size of boundaryList of Boundary object "
                       + std::to_string(cLen) + "|" + std::to_string(bLen) + " control: " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundary_object_list + level)))->get_size_boundary_left();
        cLen = get_last_index_of_boundary_slice_x(level) - get_first_index_of_boundary_slice_x(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_ny(level) * domain->get_nz(level);
            message += "length calculated by first and last index of boundary_slice_x does not equals size of boundary_slice_x of bLeft "
                       + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundary_object_list + level)))->get_size_boundary_right();
        cLen = get_last_index_of_boundary_slice_x(level) - get_first_index_of_boundary_slice_x(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_ny(level) * domain->get_nz(level);
            message += "length calculated by first and last index of boundary_slice_x does not equals size of boundary_slice_x of bRight "
                       + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundary_object_list + level)))->get_size_boundary_bottom();
        cLen = get_last_index_of_boundary_slice_y(level) - get_first_index_of_boundary_slice_y(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_nx(level) * domain->get_nz(level);
            message += "length calculated by first and last index of boundary_slice_y does not equals size of boundary_slice_y of bBottom "
                       + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundary_object_list + level)))->get_size_boundary_top();
        cLen = get_last_index_of_boundary_slice_y(level) - get_first_index_of_boundary_slice_y(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_nx(level) * domain->get_nz(level);
            message += "length calculated by first and last index of boundary_slice_y does not equals size of boundary_slice_y of bTop "
                       + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundary_object_list + level)))->get_size_boundary_front();
        cLen = get_last_index_of_boundary_slice_z(level) - get_first_index_of_boundary_slice_z(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_nx(level) * domain->get_ny(level);
            message += "length calculated by first and last index of boundary_slice_z does not equals size of boundary_slice_z of bFront "
                       + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundary_object_list + level)))->get_size_boundary_back();
        cLen = get_last_index_of_boundary_slice_z(level) - get_first_index_of_boundary_slice_z(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_nx(level) * domain->get_ny(level);
            message += "length calculated by first and last index of boundary_slice_z does not equals size of boundary_slice_z of bBack "
                       + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }

        size_t csize_inner = (domain->get_nx(level) - 2) * (domain->get_ny(level) - 2) * (domain->get_nz(level) - 2)
                             - get_size_obstacle_index_list(level);
        size_t bsize_inner = get_size_inner_list(level);
        if (csize_inner != bsize_inner) {
            message += "get_size_inner_list(level) does not equal (nx-2)*(ny-2)*(nz-2) "
                       + std::to_string(bsize_inner) + "|" + std::to_string(csize_inner) + "\n";
        }
        size_t cindex_inner_start = get_inner_list_level_joined_start(level);
        size_t cindex_inner_end = get_inner_list_level_joined_end(level);
        if (cindex_inner_end - cindex_inner_start + 1 != bsize_inner) {
            message += "get_size_inner_list(level) does not equal the difference between start and end "
                       + std::to_string(cindex_inner_start) + "|" + std::to_string(cindex_inner_end) + "\n";
        }
        size_t bsize_boundary = nx * ny * nz - (nx - 2) * (ny - 2) * (nz - 2);
        size_t csize_boundary = get_size_boundary_ist(level);
        if (csize_boundary != bsize_boundary) {
            message += "get_size_boundary_list(level) does not equal size-(nx-2)*(ny-2)*(nz-2) "
                       + std::to_string(bsize_boundary) + "|" + std::to_string(csize_boundary) + "\n";
        }

        if (!message.empty()) {
            message = "For Level " + std::to_string(level) + "\nsize control\n" + message;
        }
    }
    for (size_t level = 0; level < m_levels + 1; level++) {
        Boundary *b = *(m_MG_boundary_object_list + level);
        b->control(get_size_obstacle_index_list(level));
    }
    {
        size_t bsize_inner = 0;
        for (size_t level = 0; level < m_levels + 1; level++) {
            bsize_inner += get_size_inner_list(level);
        }
        size_t csize_inner = get_size_inner_list_level_joined();
        if (bsize_inner != csize_inner) {
            message += "get_size_inner_list_level_joined does not equal the sum of each inner list "
                       + std::to_string(bsize_inner) + "|" + std::to_string(csize_inner) + "\n";
        }
    }

    if (get_size_obstacle_list() != 0) {
        std::vector<size_t> v(m_data_MG_obstacle_list_zero_joined, m_data_MG_obstacle_list_zero_joined + get_size_obstacle_list());
        std::sort(v.begin(), v.end());
        std::vector<size_t>::iterator it = std::unique(v.begin(), v.begin() + get_size_obstacle_list());
        v.resize(std::distance(v.begin(), it));
        if (v.size() != get_size_obstacle_list()) {
            message += "obstacles are not allowed to overlap in level 0! difference: "
                       + std::to_string(get_size_obstacle_list() - v.size()) + "\n";

            for (size_t i = 0; i < m_number_of_obstacle_objects; i++) {
                Obstacle *o1 = m_MG_obstacle_object_list[0][i];
                size_t size1 = o1->get_size_obstacle_list();
                for (size_t j = i + 1; j < m_number_of_obstacle_objects; j++) {
                    Obstacle *o2 = m_MG_obstacle_object_list[0][j];
                    size_t size2 = o2->get_size_obstacle_list();
                    std::vector<size_t> tmp = Utility::mergeSortedListsToUniqueList(o1->get_obstacle_list(), size1, o2->get_obstacle_list(), size2);
                    if (tmp.size() != size1 + size2) {
                        message += "Obstacles " + std::to_string(i) + " and " + std::to_string(j) + "are overlapping\n";
                    }
                }
            }
        }
    }
    if (!message.empty()) {
        message = "################ MULTIGRID CONTROL ################\n" + message
                  + "---------------- MULTIGRID CONTROL END ----------------";
#ifndef BENCHMARKING
        m_logger->warn(message);
#endif
    }
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
    m_logger->debug("Levels: {}", m_levels);

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} inner_list starts at index {} and ends with index {}", level,
                        get_first_index_of_inner_index_list(level), get_last_index_of_inner_index_list(level));
        m_logger->debug("and the corresponding indices at this position: {}|{}",
                        *(m_data_MG_inner_list_level_joined + get_first_index_of_inner_index_list(level)),
                        *(m_data_MG_inner_list_level_joined + get_last_index_of_inner_index_list(level)));
    }
    m_logger->debug("Total length of inner_list: {}", get_length_of_inner_index_list_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} boundary_list starts at index {} and ends with index {}", level,
                        get_first_index_of_boundary_index_list(level), get_last_index_boundary_index_list(level));
        m_logger->debug("and the corresponding indices at this position: {} | {}",
                        *(m_data_MG_boundary_list_level_joined + get_first_index_of_boundary_index_list(level)),
                        *(m_data_MG_boundary_list_level_joined + get_last_index_boundary_index_list(level)));
    }
    m_logger->debug("Total length of boundary_list: {}", get_length_of_boundary_index_list_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} boundary_slice_z starts at index {} and ends with index {}", level,
                        get_first_index_of_boundary_slice_z(level), get_last_index_of_boundary_slice_z(level));
        m_logger->debug(" and the corresponding indices at this position for FRONT: {} | {}",
                        *(m_data_MG_boundary_front_level_joined + get_first_index_of_boundary_slice_z(level)),
                        *(m_data_MG_boundary_front_level_joined + get_last_index_of_boundary_slice_z(level)));
        m_logger->debug(" and the corresponding indices at this position for BACK : {} | {}",
                        *(m_data_MG_boundary_back_level_joined + get_first_index_of_boundary_slice_z(level)),
                        *(m_data_MG_boundary_back_level_joined + get_last_index_of_boundary_slice_z(level)));
    }
    m_logger->debug("Total length of boundary_slice_z: {}", get_length_of_boundary_slice_z_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} boundary_slice_y starts at index {} and ends with index {}", level,
                        get_first_index_of_boundary_slice_y(level), get_last_index_of_boundary_slice_y(level));
        m_logger->debug(" and the corresponding indices at this position for BOTTOM: {} | {}",
                        *(m_data_MG_boundary_bottom_level_joined + get_first_index_of_boundary_slice_y(level)),
                        *(m_data_MG_boundary_bottom_level_joined + get_last_index_of_boundary_slice_y(level)));
        m_logger->debug(" and the corresponding indices at this position for TOP   : {} | {}",
                        *(m_data_MG_boundary_top_level_joined + get_first_index_of_boundary_slice_y(level)),
                        *(m_data_MG_boundary_top_level_joined + get_last_index_of_boundary_slice_y(level)));
    }
    m_logger->debug("Total length of boundary_slice_y: {}", get_length_of_boundary_slice_y_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} boundary_slice_x starts at index {} and ends with index {}", level,
                        get_first_index_of_boundary_slice_x(level), get_last_index_of_boundary_slice_x(level));
        m_logger->debug(" and the corresponding indices at this position for LEFT : {} | {}",
                        *(m_data_MG_boundary_left_level_joined + get_first_index_of_boundary_slice_x(level)),
                        *(m_data_MG_boundary_left_level_joined + get_last_index_of_boundary_slice_x(level)));
        m_logger->debug(" and the corresponding indices at this position for RIGHT: {} | {}",
                        *(m_data_MG_boundary_right_level_joined + get_first_index_of_boundary_slice_x(level)),
                        *(m_data_MG_boundary_right_level_joined + get_last_index_of_boundary_slice_x(level)));
    }
    m_logger->debug("Total length of boundary_slice_x: {}", get_length_of_boundary_slice_x_joined());

    for (size_t level = 0; level < m_levels; level++) {
        m_logger->debug("For Level {} obstacle_list has {} elements",
                        level, get_size_obstacle_index_list(level));
    }

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oFront starts at index {} and ends with index {} with length {}",
                        level,
                        get_first_index_of_obstacle_front(level, 0),
                        get_last_index_of_obstacle_front(level, m_number_of_obstacle_objects - 1),
                        get_length_of_obstacle_front(level));
    }
    m_logger->debug("Total length of oFront: ", get_length_of_obstacle_front_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oBack starts at index {} and ends with index {} with length {}",
                        level,
                        get_first_index_of_obstacle_back(level, 0),
                        get_last_index_of_obstacle_back(level, m_number_of_obstacle_objects - 1),
                        get_length_of_obstacle_back(level));
    }
    m_logger->debug("Total length of oBack: {}", get_length_of_obstacle_back_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oBottom starts at index {} and ends with index {} with length {}",
                        level,
                        get_first_index_of_obstacle_bottom(level, 0),
                        get_last_index_of_obstacle_bottom(level, m_number_of_obstacle_objects - 1),
                        get_length_of_obstacle_bottom(level));
    }
    m_logger->debug("Total length of oBottom: {}", get_length_of_obstacle_bottom_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oTop starts at index {} and ends with index {} with length {}",
                        level,
                        get_first_index_of_obstacle_top(level, 0),
                        get_last_index_of_obstacle_top(level, m_number_of_obstacle_objects - 1),
                        get_length_of_obstacle_top(level));
    }
    m_logger->debug("Total length of oTop: {}", get_length_of_obstacle_top_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oLeft starts at index {} and ends with index {} with length {}",
                        level,
                        get_first_index_of_obstacle_left(level, 0),
                        get_last_index_of_obstacle_left(level, m_number_of_obstacle_objects - 1),
                        get_length_of_obstacle_left(level));
    }
    m_logger->debug("Total length of oLeft: {}", get_length_of_obstacle_left_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oRight starts at index {} and ends with index {} with length {}",
                        level,
                        get_first_index_of_obstacle_right(level, 0),
                        get_last_index_of_obstacle_right(level, m_number_of_obstacle_objects - 1),
                        get_length_of_obstacle_right(level));
    }
    m_logger->debug("Total length of oRight: {}", get_length_of_obstacle_right_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
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
    // create Boundary object of level 0
    Boundary *b;
    if (m_number_of_obstacle_objects > 0) {
        calc_obstacles(*(m_MG_obstacle_object_list));
        b = new Boundary(*(m_MG_obstacle_object_list), m_number_of_obstacle_objects,
                         get_size_obstacle_index_list(0));
    } else {
        b = new Boundary();
    }
    // set size of respective lists
    m_size_MG_inner_index_list_level[1] = b->get_size_inner_list();
    m_size_MG_boundary_index_list_level[1] = b->get_size_boundary_list();
    m_size_MG_boundary_slice_z_level[1] = b->get_size_boundary_front();
    m_size_MG_boundary_slice_y_level[1] = b->get_size_boundary_top();
    m_size_MG_boundary_slice_x_level[1] = b->get_size_boundary_left();
    // save boundary object in multigrid list
    *(m_MG_boundary_object_list) = b;

    if (m_number_of_surface_objects > 0) {
        calc_surfaces(*(m_MG_surface_object_list));
    }

    // create boundary object, surfaces and obstacles for each multigrid level
    for (size_t level = 1; level < m_levels + 1; level++) {
        surface_dominant_restriction(level);
        Obstacle **obstacle_list = obstacle_dominant_restriction(level);

        Boundary *boundary;
        if (m_number_of_obstacle_objects > 0) {
            boundary = new Boundary(obstacle_list, m_number_of_obstacle_objects,
                                    get_size_obstacle_index_list(level), level);
        } else {
            boundary = new Boundary(level);
        }
        *(m_MG_boundary_object_list + level) = boundary;

        m_size_MG_inner_index_list_level[level + 1] = m_size_MG_inner_index_list_level[level]
                                                      + boundary->get_size_inner_list();
        m_size_MG_boundary_index_list_level[level + 1] = m_size_MG_boundary_index_list_level[level]
                                                         + boundary->get_size_boundary_list();
        m_size_MG_boundary_slice_z_level[level + 1] = m_size_MG_boundary_slice_z_level[level]
                                                      + boundary->get_size_boundary_front();
        m_size_MG_boundary_slice_y_level[level + 1] = m_size_MG_boundary_slice_y_level[level]
                                                      + boundary->get_size_boundary_top();
        m_size_MG_boundary_slice_x_level[level + 1] = m_size_MG_boundary_slice_x_level[level]
                                                      + boundary->get_size_boundary_left();
    }
}

// ================================= Calc obstacles ================================================
// *************************************************************************************************
/// \brief  create obstacles of level 0, no check if obstacles are overlapping
/// \param obstacle_object_list List of obstacle objects
// *************************************************************************************************
void Multigrid::calc_obstacles(Obstacle **obstacle_object_list) {
    if (m_number_of_obstacle_objects > 0) {
        size_t level = 0;
        size_t o_size = get_size_obstacle_index_list(level);
        size_t *obstacle_index_list = new size_t[o_size];
        size_t counter = 0;
        for (size_t o = 0; o < m_number_of_obstacle_objects; o++) {
            Obstacle *obstacle_tmp = obstacle_object_list[o];
            size_t len_all = obstacle_tmp->get_size_obstacle_list();
            for (size_t i = 0; i < len_all; i++) {
                *(obstacle_index_list + counter) = obstacle_tmp->get_obstacle_list()[i];
                counter++;
            }
        }
        *(m_MG_obstacle_index_list) = obstacle_index_list;
    }
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
        Domain *domain = Domain::getInstance();
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

            Surface *surface_coarse = new Surface(surfaceID, startIndex_coarse,
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
Obstacle **Multigrid::obstacle_dominant_restriction(size_t level) {
// OBSTACLES
    // add index to m_obstacle_lists if any of l-1 indices building the l index was an obstacle
    // (dominant restriction)
    if (m_number_of_obstacle_objects == 0) {
        return nullptr;
    }
    // TODO define lists with obstacle size
    Domain *domain = Domain::getInstance();
    Obstacle **obstacle_list_fine = *(m_MG_obstacle_object_list + (level - 1));
    Obstacle **obstacle_list_coarse = new Obstacle *[m_number_of_obstacle_objects];
    *(m_MG_obstacle_object_list + level) = obstacle_list_coarse;
    *(m_size_MG_obstacle_index_list_level + level) = 0;
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
        if (i2_fine - i1_fine + 1 < domain->get_nx(level - 1) - 2
            && i2_coarse - i1_coarse + 1 >= domain->get_nx(level) - 2) {
            m_logger->warn("Be cautious! Obstacle fills up inner cells in x-direction at level {}", level);
        }
        if (j2_fine - j1_fine + 1 < domain->get_ny(level - 1) - 2
            && j2_coarse - j1_coarse + 1 >= domain->get_ny(level) - 2) {
            m_logger->warn("Be cautious! Obstacle fills up inner cells in y-direction at level {}", level);
        }
        if (k2_fine - k1_fine + 1 < domain->get_nz(level - 1) - 2
            && k2_coarse - k1_coarse + 1 >= domain->get_nz(level) - 2) {
            m_logger->warn("Be cautious! Obstacle fills up inner cells in z-direction at level {}", level);
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
        *(obstacle_list_coarse + id) = obstacle_coarse;

        size_t index = level * m_number_of_obstacle_objects + id + 1;
        m_size_MG_obstacle_front_level[index] = obstacle_coarse->get_size_obstacle_front()
                                                + m_size_MG_obstacle_front_level[index - 1];
        m_size_MG_obstacle_back_level[index] = obstacle_coarse->get_size_obstacle_back()
                                               + m_size_MG_obstacle_back_level[index - 1];
        m_size_MG_obstacle_bottom_level[index] = obstacle_coarse->get_size_obstacle_bottom()
                                                 + m_size_MG_obstacle_bottom_level[index - 1];
        m_size_MG_obstacle_top_level[index] = obstacle_coarse->get_size_obstacle_top()
                                              + m_size_MG_obstacle_top_level[index - 1];
        m_size_MG_obstacle_left_level[index] = obstacle_coarse->get_size_obstacle_left()
                                               + m_size_MG_obstacle_left_level[index - 1];
        m_size_MG_obstacle_right_level[index] = obstacle_coarse->get_size_obstacle_right()
                                                + m_size_MG_obstacle_right_level[index - 1];
        *(m_size_MG_obstacle_index_list_level + level) += obstacle_coarse->get_size_obstacle_list();
    } //end obstacle id loop

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
    *(m_size_MG_obstacle_index_list_level + level) = size;
    *(m_MG_obstacle_index_list + level) = obstacle_list_tmp;
    return obstacle_list_coarse;
}

// ================================= Send lists to GPU =============================================
// *************************************************************************************************
/// \brief  create joined list and send them to GPU
// *************************************************************************************************
void Multigrid::send_lists_to_GPU() {
    send_surface_lists_to_GPU();
    send_boundary_lists_to_GPU();
    send_obstacle_lists_to_GPU();
}

// ================================= Send boundary lists to GPU ====================================
// *************************************************************************************************
/// \brief  create boundary joined list and send them to GPU
// *************************************************************************************************
void Multigrid::send_boundary_lists_to_GPU() {
    size_t size_inner_list = get_length_of_inner_index_list_joined();
    size_t size_boundary_list = get_length_of_boundary_index_list_joined();
    size_t size_boundary_slice_z = get_length_of_boundary_slice_z_joined();
    size_t size_boundary_slice_y = get_length_of_boundary_slice_y_joined();
    size_t size_boundary_slice_x = get_length_of_boundary_slice_x_joined();

    m_data_MG_inner_list_level_joined = new size_t[size_inner_list];
    m_data_MG_boundary_list_level_joined = new size_t[size_boundary_list];
    m_data_MG_boundary_front_level_joined = new size_t[size_boundary_slice_z];
    m_data_MG_boundary_back_level_joined = new size_t[size_boundary_slice_z];
    m_data_MG_boundary_top_level_joined = new size_t[size_boundary_slice_y];
    m_data_MG_boundary_bottom_level_joined = new size_t[size_boundary_slice_y];
    m_data_MG_boundary_left_level_joined = new size_t[size_boundary_slice_x];
    m_data_MG_boundary_right_level_joined = new size_t[size_boundary_slice_x];

    size_t counter_inner_list = 0;
    size_t counter_boundary_list = 0;

    size_t counter_boundary_slice_z = 0;
    size_t counter_boundary_slice_y = 0;
    size_t counter_boundary_slice_x = 0;

    for (size_t level = 0; level < m_levels + 1; level++) {
        Boundary *boundary = *(m_MG_boundary_object_list + level);
        for (size_t i = 0; i < boundary->get_size_inner_list(); i++) {
            *(m_data_MG_inner_list_level_joined + counter_inner_list) = boundary->get_inner_list()[i];
            counter_inner_list++;
        }
        for (size_t i = 0; i < boundary->get_size_boundary_list(); i++) {
            *(m_data_MG_boundary_list_level_joined + counter_boundary_list) = boundary->get_boundary_list()[i];
            counter_boundary_list++;
        }
        for (size_t i = 0; i < boundary->get_size_boundary_front(); i++) {
            *(m_data_MG_boundary_front_level_joined + counter_boundary_slice_z) = boundary->get_boundary_front()[i];
            *(m_data_MG_boundary_back_level_joined + counter_boundary_slice_z) = boundary->get_boundary_back()[i];
            counter_boundary_slice_z++;
        }
        for (size_t i = 0; i < boundary->get_size_boundary_top(); i++) {
            *(m_data_MG_boundary_top_level_joined + counter_boundary_slice_y) = boundary->get_boundary_top()[i];
            *(m_data_MG_boundary_bottom_level_joined + counter_boundary_slice_y) = boundary->get_boundary_bottom()[i];
            counter_boundary_slice_y++;
        }
        for (size_t i = 0; i < boundary->get_size_boundary_left(); i++) {
            *(m_data_MG_boundary_left_level_joined + counter_boundary_slice_x) = boundary->get_boundary_left()[i];
            *(m_data_MG_boundary_right_level_joined + counter_boundary_slice_x) = boundary->get_boundary_right()[i];
            counter_boundary_slice_x++;
        }
    }
#pragma acc enter data copyin(m_data_MG_inner_list_level_joined[:size_inner_list])
#pragma acc enter data copyin(m_data_MG_boundary_list_level_joined[:size_boundary_list])
#pragma acc enter data copyin(m_data_MG_boundary_front_level_joined[:size_boundary_slice_z])
#pragma acc enter data copyin(m_data_MG_boundary_back_level_joined[:size_boundary_slice_z])
#pragma acc enter data copyin(m_data_MG_boundary_top_level_joined[:size_boundary_slice_y])
#pragma acc enter data copyin(m_data_MG_boundary_bottom_level_joined[:size_boundary_slice_y])
#pragma acc enter data copyin(m_data_MG_boundary_left_level_joined[:size_boundary_slice_x])
#pragma acc enter data copyin(m_data_MG_boundary_right_level_joined[:size_boundary_slice_x])
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
        for (size_t level = 0; level < m_levels + 1; level++) {
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
#pragma acc enter data copyin(m_data_MG_surface_list_level_joined[:size_surface_list])
    }
}

// ================================= Send obstacle lists to GPU ====================================
// *************************************************************************************************
/// \brief  create obstacle joined list and send them to GPU
// *************************************************************************************************
void Multigrid::send_obstacle_lists_to_GPU() {
    if (m_number_of_obstacle_objects > 0) {
        size_t size_oFront = get_length_of_obstacle_front_joined();
        size_t size_oBack = get_length_of_obstacle_back_joined();
        size_t size_oBottom = get_length_of_obstacle_bottom_joined();
        size_t size_oTop = get_length_of_obstacle_top_joined();
        size_t size_oLeft = get_length_of_obstacle_left_joined();
        size_t size_oRight = get_length_of_obstacle_right_joined();

        m_data_MG_obstacle_front_level_joined = new size_t[size_oFront];
        m_data_MG_obstacle_back_level_joined = new size_t[size_oBack];
        m_data_MG_obstacle_bottom_level_joined = new size_t[size_oBottom];
        m_data_MG_obstacle_top_level_joined = new size_t[size_oTop];
        m_data_MG_obstacle_left_level_joined = new size_t[size_oLeft];
        m_data_MG_obstacle_right_level_joined = new size_t[size_oRight];

        size_t counter_oFront = 0;
        size_t counter_oBack = 0;
        size_t counter_oBottom = 0;
        size_t counter_oTop = 0;
        size_t counter_oLeft = 0;
        size_t counter_oRight = 0;

        for (size_t level = 0; level < m_levels + 1; level++) {
            Obstacle **obstacle_list = *(m_MG_obstacle_object_list + level);
            for (size_t o = 0; o < m_number_of_obstacle_objects; o++) {
                Obstacle *obstacle = *(obstacle_list + o);
                for (size_t i = 0; i < obstacle->get_size_obstacle_front(); i++) {
                    m_data_MG_obstacle_front_level_joined[counter_oFront] = obstacle->get_obstacle_front()[i];
                    counter_oFront++;
                }
                for (size_t i = 0; i < obstacle->get_size_obstacle_back(); i++) {
                    *(m_data_MG_obstacle_back_level_joined + counter_oBack) = obstacle->get_obstacle_back()[i];
                    counter_oBack++;
                }
                for (size_t i = 0; i < obstacle->get_size_obstacle_bottom(); i++) {
                    *(m_data_MG_obstacle_bottom_level_joined + counter_oBottom) = obstacle->get_obstacle_bottom()[i];
                    counter_oBottom++;
                }
                for (size_t i = 0; i < obstacle->get_size_obstacle_top(); i++) {
                    *(m_data_MG_obstacle_top_level_joined + counter_oTop) = obstacle->get_obstacle_top()[i];
                    counter_oTop++;
                }
                for (size_t i = 0; i < obstacle->get_size_obstacle_left(); i++) {
                    *(m_data_MG_obstacle_left_level_joined + counter_oLeft) = obstacle->get_obstacle_left()[i];
                    counter_oLeft++;
                }
                for (size_t i = 0; i < obstacle->get_size_obstacle_right(); i++) {
                    *(m_data_MG_obstacle_right_level_joined + counter_oRight) = obstacle->get_obstacle_right()[i];
                    counter_oRight++;
                }
            }
        }

        m_data_MG_obstacle_list_zero_joined = m_MG_obstacle_index_list[0];  // TODO(issue 33) wrong because only one obstacle is used?
        size_t size_obstacle_list = get_size_obstacle_list();
#pragma acc enter data copyin(m_data_MG_obstacle_front_level_joined[:size_oFront])
#pragma acc enter data copyin(m_data_MG_obstacle_back_level_joined[:size_oBack])
#pragma acc enter data copyin(m_data_MG_obstacle_top_level_joined[:size_oTop])
#pragma acc enter data copyin(m_data_MG_obstacle_bottom_level_joined[:size_oBottom])
#pragma acc enter data copyin(m_data_MG_obstacle_left_level_joined[:size_oLeft])
#pragma acc enter data copyin(m_data_MG_obstacle_right_level_joined[:size_oRight])
#pragma acc enter data copyin(m_data_MG_obstacle_list_zero_joined[:size_obstacle_list])
    }
}

// ================================= Apply boundary condition ======================================
// *************************************************************************************************
/// \brief  Apply boundary condition for obstacles, surfaces and domain
/// \param  d Field
/// \param  level Multigrid level
/// \param  f Field type
/// \param  sync synchronous kernel launching (true, default: false)
// *************************************************************************************************
void Multigrid::apply_boundary_condition(real *d, size_t level, FieldType f, bool sync) {
    if (m_number_of_surface_objects > 0) {
        Surface **surface_list = *(m_MG_surface_object_list + level);
        for (size_t id = 0; id < m_number_of_surface_objects; ++id) {
            (static_cast<Surface *>(*(surface_list + id)))->applyBoundaryConditions(d, f, level, sync);
        }
    }
    if (m_number_of_obstacle_objects > 0) {
        for (size_t id = 0; id < m_number_of_obstacle_objects - 1; ++id) {
            size_t opatch_start[] = {get_first_index_of_obstacle_front(level, id),
                                     get_first_index_of_obstacle_back(level, id),
                                     get_first_index_of_obstacle_bottom(level, id),
                                     get_first_index_of_obstacle_top(level, id),
                                     get_first_index_of_obstacle_left(level, id),
                                     get_first_index_of_obstacle_right(level, id)};
            size_t opatch_end[] = {get_first_index_of_obstacle_front(level, id + 1),
                                   get_first_index_of_obstacle_back(level, id + 1),
                                   get_first_index_of_obstacle_bottom(level, id + 1),
                                   get_first_index_of_obstacle_top(level, id + 1),
                                   get_first_index_of_obstacle_left(level, id + 1),
                                   get_first_index_of_obstacle_right(level, id + 1)};
            (static_cast<BoundaryDataController *> (*(m_bdc_obstacle + id)))->apply_boundary_condition_obstacle(d, m_data_obstacles_patches_joined, opatch_start, opatch_end, f, level, id, sync);
        }
        size_t opatch_start[] = {get_first_index_of_obstacle_front(level, m_number_of_obstacle_objects - 1),
                                 get_first_index_of_obstacle_back(level, m_number_of_obstacle_objects - 1),
                                 get_first_index_of_obstacle_bottom(level, m_number_of_obstacle_objects - 1),
                                 get_first_index_of_obstacle_top(level, m_number_of_obstacle_objects - 1),
                                 get_first_index_of_obstacle_left(level, m_number_of_obstacle_objects - 1),
                                 get_first_index_of_obstacle_right(level, m_number_of_obstacle_objects - 1)};
        size_t opatch_end[] = {get_first_index_of_obstacle_front(level + 1, 0),
                               get_first_index_of_obstacle_back(level + 1, 0),
                               get_first_index_of_obstacle_bottom(level + 1, 0),
                               get_first_index_of_obstacle_top(level + 1, 0),
                               get_first_index_of_obstacle_left(level + 1, 0),
                               get_first_index_of_obstacle_right(level + 1, 0)};
        (static_cast<BoundaryDataController *>(*(m_bdc_obstacle + m_number_of_obstacle_objects - 1)))->apply_boundary_condition_obstacle(d, m_data_obstacles_patches_joined, opatch_start, opatch_end, f, level, m_number_of_obstacle_objects - 1, sync);
    }

    size_t patch_start[] = {get_first_index_of_boundary_slice_z(level), get_first_index_of_boundary_slice_z(level),
                            get_first_index_of_boundary_slice_y(level), get_first_index_of_boundary_slice_y(level),
                            get_first_index_of_boundary_slice_x(level), get_first_index_of_boundary_slice_x(level)};
    size_t patch_end[] = {get_first_index_of_boundary_slice_z(level + 1), get_first_index_of_boundary_slice_z(level + 1),
                          get_first_index_of_boundary_slice_y(level + 1), get_first_index_of_boundary_slice_y(level + 1),
                          get_first_index_of_boundary_slice_x(level + 1), get_first_index_of_boundary_slice_x(level + 1)};
    m_bdc_boundary->apply_boundary_condition(d, m_data_boundary_patches_joined, patch_start, patch_end, f, level, sync);
}

//======================================== Update lists ====================================
// ***************************************************************************************
/// \brief  Updates lists of indices
// ***************************************************************************************
void Multigrid::update_lists() {
    remove_boundary_lists_from_GPU();

    *(m_size_MG_inner_index_list_level) = 0;
    *(m_size_MG_boundary_index_list_level) = 0;

    *(m_size_MG_boundary_slice_z_level) = 0;
    *(m_size_MG_boundary_slice_y_level) = 0;
    *(m_size_MG_boundary_slice_x_level) = 0;

    if (m_number_of_obstacle_objects > 0) {
        for (size_t level = 0; level < m_levels + 1; level++) {
            Boundary *boundary = *(m_MG_boundary_object_list + level);
            boundary->update_lists(*(m_MG_obstacle_object_list + level),
                                   m_number_of_obstacle_objects,
                                   get_size_obstacle_index_list(level));
            m_size_MG_inner_index_list_level[level + 1] = m_size_MG_inner_index_list_level[level]
                                                          + boundary->get_size_inner_list();
            m_size_MG_boundary_index_list_level[level + 1] = m_size_MG_boundary_index_list_level[level]
                                                             + boundary->get_size_boundary_list();
            m_size_MG_boundary_slice_z_level[level + 1] = m_size_MG_boundary_slice_z_level[level]
                                                          + boundary->get_size_boundary_front();
            m_size_MG_boundary_slice_y_level[level + 1] = m_size_MG_boundary_slice_y_level[level]
                                                          + boundary->get_size_boundary_top();
            m_size_MG_boundary_slice_x_level[level + 1] = m_size_MG_boundary_slice_x_level[level]
                                                          + boundary->get_size_boundary_left();
        }
    } else {
        for (size_t level = 0; level < m_levels + 1; level++) {
            Boundary *boundary = *(m_MG_boundary_object_list + level);
            boundary->update_lists();
            m_size_MG_inner_index_list_level[level + 1] = m_size_MG_inner_index_list_level[level]
                                                          + boundary->get_size_inner_list();
            m_size_MG_boundary_index_list_level[level + 1] = m_size_MG_boundary_index_list_level[level]
                                                             + boundary->get_size_boundary_list();
            m_size_MG_boundary_slice_z_level[level + 1] = m_size_MG_boundary_slice_z_level[level]
                                                          + boundary->get_size_boundary_front();
            m_size_MG_boundary_slice_y_level[level + 1] = m_size_MG_boundary_slice_y_level[level]
                                                          + boundary->get_size_boundary_top();
            m_size_MG_boundary_slice_x_level[level + 1] = m_size_MG_boundary_slice_x_level[level]
                                                          + boundary->get_size_boundary_left();
        }
    }
    send_boundary_lists_to_GPU();
    m_data_boundary_patches_joined[Patch::FRONT] = m_data_MG_boundary_front_level_joined;
    m_data_boundary_patches_joined[Patch::BACK] = m_data_MG_boundary_back_level_joined;
    m_data_boundary_patches_joined[Patch::BOTTOM] = m_data_MG_boundary_bottom_level_joined;
    m_data_boundary_patches_joined[Patch::TOP] = m_data_MG_boundary_top_level_joined;
    m_data_boundary_patches_joined[Patch::LEFT] = m_data_MG_boundary_left_level_joined;
    m_data_boundary_patches_joined[Patch::RIGHT] = m_data_MG_boundary_right_level_joined;
}

void Multigrid::remove_boundary_lists_from_GPU() {
    size_t size_inner_list = get_length_of_inner_index_list_joined();
    size_t size_boundary_list = get_length_of_boundary_index_list_joined();
#pragma acc exit data delete(m_data_MG_inner_list_level_joined[:size_inner_list])
#pragma acc exit data delete(m_data_MG_boundary_list_level_joined[:size_boundary_list])
    delete[] m_data_MG_inner_list_level_joined;
    delete[] m_data_MG_boundary_list_level_joined;

    size_t size_boundary_slice_z = get_length_of_boundary_slice_z_joined();
    size_t size_boundary_slice_y = get_length_of_boundary_slice_y_joined();
    size_t size_boundary_slice_x = get_length_of_boundary_slice_x_joined();
#pragma acc exit data delete(m_data_MG_boundary_front_level_joined[:size_boundary_slice_z])
#pragma acc exit data delete(m_data_MG_boundary_back_level_joined[:size_boundary_slice_z])
#pragma acc exit data delete(m_data_MG_boundary_top_level_joined[:size_boundary_slice_y])
#pragma acc exit data delete(m_data_MG_boundary_bottom_level_joined[:size_boundary_slice_y])
#pragma acc exit data delete(m_data_MG_boundary_left_level_joined[:size_boundary_slice_x])
#pragma acc exit data delete(m_data_MG_boundary_right_level_joined[:size_boundary_slice_x])
    delete[] m_data_MG_boundary_front_level_joined;
    delete[] m_data_MG_boundary_back_level_joined;
    delete[] m_data_MG_boundary_top_level_joined;
    delete[] m_data_MG_boundary_bottom_level_joined;
    delete[] m_data_MG_boundary_left_level_joined;
    delete[] m_data_MG_boundary_right_level_joined;
}

//======================================== Private getter ==========================================

// *************************************************************************************************
/// \brief  get size of obstacle index list for level
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_size_obstacle_index_list(size_t level) const {
    size_t size_obstacle_list = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_obstacle_list = *(m_size_MG_obstacle_index_list_level + level);
    }
    return size_obstacle_list;
}

// *************************************************************************************************
/// \brief  get length of the joined inner cell list
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_inner_index_list_joined() const {
    return get_first_index_of_inner_index_list(m_levels + 1) + 1;
}

// *************************************************************************************************
/// \brief  get the index of joined inner cell list, where the first cell of level l is
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_first_index_of_inner_index_list(size_t level) const {
    return *(m_size_MG_inner_index_list_level + level);
}

// *************************************************************************************************
/// \brief  get the index of joined inner cell list, where the last cell of level l is
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_last_index_of_inner_index_list(size_t level) const {
    return *(m_size_MG_inner_index_list_level + level + 1) - 1;
}

// boundary_list

// *************************************************************************************************
/// \brief  get the index of joined inner cell list, where the first cell of level l is
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_boundary_index_list_joined() const {
    return get_first_index_of_boundary_index_list(m_levels + 1) + 1;
}

// *************************************************************************************************
/// \brief  get the index of joined boundary cell list, where the first cell of level l is
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_first_index_of_boundary_index_list(size_t level) const {
    return *(m_size_MG_boundary_index_list_level + level);
}

// *************************************************************************************************
/// \brief  get the index of joined boundary cell list, where the last cell of level l is
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_last_index_boundary_index_list(size_t level) const {
    return *(m_size_MG_boundary_index_list_level + level + 1) - 1;
}

// surface_list
size_t Multigrid::get_length_of_surface_index_list_joined() const {
    size_t len = 0;
    if (m_number_of_surface_objects > 0) {
        len = get_first_index_of_surface_index_list(m_levels + 1) + 1;
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

// bSlices

// *************************************************************************************************
/// \brief  get length of left/right joined boundary list
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_boundary_slice_x_joined() const {
    return get_first_index_of_boundary_slice_x(m_levels + 1) + 1;
}

// *************************************************************************************************
/// \brief  get length of top/bottom joined boundary list
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_boundary_slice_y_joined() const {
    return get_first_index_of_boundary_slice_y(m_levels + 1) + 1;
}

// *************************************************************************************************
/// \brief  get length of front/back joined boundary list
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_boundary_slice_z_joined() const {
    return get_first_index_of_boundary_slice_z(m_levels + 1) + 1;
}

// *************************************************************************************************
/// \brief  get the index of joined boundary cell list of left/right patch of the first cell
/// of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_first_index_of_boundary_slice_x(size_t level) const {
    return *(m_size_MG_boundary_slice_x_level + level);
}

// *************************************************************************************************
/// \brief  get the index of joined boundary cell list of top/bottom patch of the first cell
/// of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_first_index_of_boundary_slice_y(size_t level) const {
    return *(m_size_MG_boundary_slice_y_level + level);
}

// *************************************************************************************************
/// \brief  get the index of joined boundary cell list of front/back patch of the first cell
/// of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_first_index_of_boundary_slice_z(size_t level) const {
    return *(m_size_MG_boundary_slice_z_level + level);
}

// *************************************************************************************************
/// \brief  get the index of joined boundary cell list of left/right patch of the last cell
/// of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_last_index_of_boundary_slice_x(size_t level) const {
    return *(m_size_MG_boundary_slice_x_level + level + 1) - 1;
}

// *************************************************************************************************
/// \brief  get the index of joined boundary cell list of top/bottom patch of the last cell
/// of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_last_index_of_boundary_slice_y(size_t level) const {
    return *(m_size_MG_boundary_slice_y_level + level + 1) - 1;
}

// *************************************************************************************************
/// \brief  get the index of joined boundary cell list of front/back patch of the last cell
/// of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_last_index_of_boundary_slice_z(size_t level) const {
    return *(m_size_MG_boundary_slice_z_level + level + 1) - 1;
}

// oSlices

// *************************************************************************************************
/// \brief  get length of m_data_MG_obstacle_front_level_joined
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_front_joined() const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_t index = (m_levels + 1) * m_number_of_obstacle_objects;
        if (m_size_MG_obstacle_front_level[index] > 0) {
            len = m_size_MG_obstacle_front_level[index] + 1;
        }
    }
    return len;
}

// *************************************************************************************************
/// \brief  get length of m_data_MG_obstacle_back_level_joined
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_back_joined() const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_t index = (m_levels + 1) * m_number_of_obstacle_objects;
        if (m_size_MG_obstacle_back_level[index] > 0) {
            len = m_size_MG_obstacle_back_level[index] + 1;
        }
    }
    return len;
}

// *************************************************************************************************
/// \brief  get length of m_data_MG_obstacle_bottom_level_joined
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_bottom_joined() const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_t index = (m_levels + 1) * m_number_of_obstacle_objects;
        if (m_size_MG_obstacle_bottom_level[index] > 0) {
            len = m_size_MG_obstacle_bottom_level[index] + 1;
        }
    }
    return len;
}

// *************************************************************************************************
/// \brief  get length of m_data_MG_obstacle_top_level_joined
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_top_joined() const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_t index = (m_levels + 1) * m_number_of_obstacle_objects;
        if (m_size_MG_obstacle_top_level[index] > 0) {
            len = m_size_MG_obstacle_top_level[index] + 1;
        }
    }
    return len;
}

// *************************************************************************************************
/// \brief  get length of m_data_MG_obstacle_left_level_joined
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_left_joined() const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_t index = (m_levels + 1) * m_number_of_obstacle_objects;
        if (m_size_MG_obstacle_left_level[index] > 0) {
            len = m_size_MG_obstacle_left_level[index] + 1;
        }
    }
    return len;
}

// *************************************************************************************************
/// \brief  get length of m_data_MG_obstacle_right_level_joined
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_right_joined() const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_t index = (m_levels + 1) * m_number_of_obstacle_objects;
        if (m_size_MG_obstacle_right_level[index] > 0) {
            len = m_size_MG_obstacle_right_level[index] + 1;
        }
    }
    return len;
}

// *************************************************************************************************
/// \brief  get amount of Front obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_front(size_t level) const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_t last = get_last_index_of_obstacle_front(level, m_number_of_obstacle_objects - 1);
        if (last > 0) {
            len = last - get_first_index_of_obstacle_front(level, 0) + 1;
        }
    }
    return len;
}

// *************************************************************************************************
/// \brief get amount of Back obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_back(size_t level) const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_t last = get_last_index_of_obstacle_back(level, m_number_of_obstacle_objects - 1);
        if (last > 0) {
            len = last - get_first_index_of_obstacle_back(level, 0) + 1;
        }
    }
    return len;
}

// *************************************************************************************************
/// \brief  get amount of Bottom obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_bottom(size_t level) const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_t last = get_last_index_of_obstacle_bottom(level, m_number_of_obstacle_objects - 1);
        if (last > 0) {
            len = last - get_first_index_of_obstacle_bottom(level, 0) + 1;
        }
    }
    return len;
}

// *************************************************************************************************
/// \brief  get amount of Top obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_top(size_t level) const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_t last = get_last_index_of_obstacle_top(level, m_number_of_obstacle_objects - 1);
        if (last > 0) {
            len = last - get_first_index_of_obstacle_top(level, 0) + 1;
        }
    }
    return len;
}

// *************************************************************************************************
/// \brief  get amount of Left obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_left(size_t level) const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        size_t last = get_last_index_of_obstacle_left(level, m_number_of_obstacle_objects - 1);
        if (last > 0) {
            len = last - get_first_index_of_obstacle_left(level, 0) + 1;
        }
    }
    return len;
}

// *************************************************************************************************
/// \brief  get amount of Front obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_length_of_obstacle_right(size_t level) const {
    size_t len = 0;
    if (m_number_of_obstacle_objects > 0) {
        len = get_first_index_of_obstacle_right(level + 1, 0)
              - get_first_index_of_obstacle_right(level, 0);
    }
    return len;
}

// *************************************************************************************************
/// \brief  get the index of joined obstacle cell list of Front patch of the first cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_first_index_of_obstacle_front(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = *(m_size_MG_obstacle_front_level + level * m_number_of_obstacle_objects + id);
    }
    return index;
}

// *************************************************************************************************
/// \brief  get the index of joined obstacle cell list of Back patch of the first cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_first_index_of_obstacle_back(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = *(m_size_MG_obstacle_back_level + level * m_number_of_obstacle_objects + id);
    }
    return index;
}

// *************************************************************************************************
/// \brief  get the index of joined obstacle cell list of Bottom patch of the first cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_first_index_of_obstacle_bottom(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = *(m_size_MG_obstacle_bottom_level + level * m_number_of_obstacle_objects + id);
    }
    return index;
}

// *************************************************************************************************
/// \brief  get the index of joined obstacle cell list of Top patch of the first cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_first_index_of_obstacle_top(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = *(m_size_MG_obstacle_top_level + level * m_number_of_obstacle_objects + id);
    }
    return index;
}

// *************************************************************************************************
/// \brief  get the index of joined obstacle cell list of Left patch of the first cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_first_index_of_obstacle_left(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = *(m_size_MG_obstacle_left_level + level * m_number_of_obstacle_objects + id);
    }
    return index;
}

// *************************************************************************************************
/// \brief get the index of joined obstacle cell list of right patch of the first cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_first_index_of_obstacle_right(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = *(m_size_MG_obstacle_right_level + level * m_number_of_obstacle_objects + id);
    }
    return index;
}

// *************************************************************************************************
/// \brief get the index of joined obstacle cell list of front patch of the last cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_last_index_of_obstacle_front(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = m_size_MG_obstacle_front_level[level * m_number_of_obstacle_objects + id + 1];
    }
    return index;
}

// *************************************************************************************************
/// \brief get the index of joined obstacle cell list of back patch of the last cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_last_index_of_obstacle_back(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = m_size_MG_obstacle_back_level[level * m_number_of_obstacle_objects + id + 1];
    }
    return index;
}

// *************************************************************************************************
/// \brief get the index of joined obstacle cell list of bottom patch of the last cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// *************************************************************************************************
size_t Multigrid::get_last_index_of_obstacle_bottom(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = m_size_MG_obstacle_bottom_level[level * m_number_of_obstacle_objects + id + 1] - 1;
    }
    return index;
}

// *************************************************************************************************
/// \brief get the index of joined obstacle cell list of top patch of the last cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return int
// *************************************************************************************************
size_t Multigrid::get_last_index_of_obstacle_top(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = m_size_MG_obstacle_top_level[level * m_number_of_obstacle_objects + id + 1] - 1;
    }
    return index;
}

// *************************************************************************************************
/// \brief get the index of joined obstacle cell list of left patch of the last cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return int
// *************************************************************************************************
size_t Multigrid::get_last_index_of_obstacle_left(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = m_size_MG_obstacle_left_level[level * m_number_of_obstacle_objects + id + 1];
    }
    return index;
}

// *************************************************************************************************
/// \brief get the index of joined obstacle cell list of right patch of the last cell of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return int
// *************************************************************************************************
size_t Multigrid::get_last_index_of_obstacle_right(size_t level, size_t id) const {
    size_t index = 0;
    if (m_number_of_obstacle_objects > 0) {
        index = m_size_MG_obstacle_right_level[level * m_number_of_obstacle_objects + id + 1];
    }
    return index;
}

// ---------------- public getter
size_t Multigrid::get_size_inner_list(size_t level) const {
    return get_last_index_of_inner_index_list(level) - get_first_index_of_inner_index_list(level) + 1;
}

size_t Multigrid::get_size_boundary_ist(size_t level) const {
    return get_last_index_boundary_index_list(level) - get_first_index_of_boundary_index_list(level) + 1;
}

size_t Multigrid::get_size_obstacle_list() const {
    return get_size_obstacle_index_list(0);
}

size_t *Multigrid::get_obstacle_list() const {
    if (m_number_of_obstacle_objects > 0) {
        return m_data_MG_obstacle_list_zero_joined;
    } else {
        return nullptr;
    }
}

size_t Multigrid::get_inner_list_level_joined_start(size_t level) const {
    return *(m_size_MG_inner_index_list_level + level);
}

size_t Multigrid::get_inner_list_level_joined_end(size_t level) const {
    return get_inner_list_level_joined_start(level + 1) - 1;
}

size_t Multigrid::get_boundary_list_level_joined_start(size_t level) const {
    return *(m_size_MG_boundary_index_list_level + level);
}

size_t Multigrid::get_boundary_list_level_joined_end(size_t level) const {
    return get_boundary_list_level_joined_start(level + 1) - 1;
}

size_t Multigrid::get_obstacle_stride_x(size_t id, size_t level) const {
    return (static_cast<Obstacle *>(m_MG_obstacle_object_list[level][id]))->get_stride_x();
}

size_t Multigrid::get_obstacle_stride_y(size_t id, size_t level) const {
    return (static_cast<Obstacle *>(m_MG_obstacle_object_list[level][id]))->get_stride_y();
}

size_t Multigrid::get_obstacle_stride_z(size_t id, size_t level) const {
    return (static_cast<Obstacle *>(m_MG_obstacle_object_list[level][id]))->get_stride_z();
}
