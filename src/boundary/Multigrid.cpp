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
    m_numberOfSurfaces = 0;
    m_numberOfObstacles = 0;

    init();
    addMGLists();
    sendListsToGPU();
#ifndef BENCHMARKING
    print();
    control();
#endif
    m_data_boundary_patches_joined = new size_t *[number_of_patches];
    m_data_boundary_patches_joined[Patch::FRONT] = m_data_MG_bFront_level_joined;
    m_data_boundary_patches_joined[Patch::BACK] = m_data_MG_bBack_level_joined;
    m_data_boundary_patches_joined[Patch::BOTTOM] = m_data_MG_bBottom_level_joined;
    m_data_boundary_patches_joined[Patch::TOP] = m_data_MG_bTop_level_joined;
    m_data_boundary_patches_joined[Patch::LEFT] = m_data_MG_bLeft_level_joined;
    m_data_boundary_patches_joined[Patch::RIGHT] = m_data_MG_bRight_level_joined;
}

Multigrid::Multigrid(size_t number_of_surfaces, Surface **surface_list, size_t number_of_obstacles, Obstacle **obstacle_list, BoundaryDataController *bdc_boundary, BoundaryDataController **bdc_obstacles) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_bdc_boundary = bdc_boundary;
    m_numberOfSurfaces = number_of_surfaces;
    m_numberOfObstacles = number_of_obstacles;

    m_levels = Domain::getInstance()->get_levels();  // multigrid level, 0 otherwise

    if (m_numberOfSurfaces > 0) {
        // list of surfaces for each level
        m_MG_surfaceList = new Surface **[m_levels + 1];
        *(m_MG_surfaceList) = surface_list;  // level 0

        // surface indices divided by level
        m_MG_sList = new size_t *[m_levels + 1];

        // start index of each surface in level joined list
        m_size_MG_sList_level = new size_t[m_levels * m_numberOfSurfaces + 1];
        // start index of first level in joined list = 0
        *(m_size_MG_sList_level) = 0;
    }

    if (m_numberOfObstacles > 0) {
        // list of obstacles for each level
        m_MG_obstacleList = new Obstacle **[m_levels + 1];
        *(m_MG_obstacleList) = obstacle_list;  // level 0
        m_size_MG_oList_level = new size_t[m_levels + 1];

        // obstacle indices divided by level
        m_MG_oList = new size_t *[m_levels + 1];
        // m_MG_oFront = new size_t *[m_levels + 1];
        // m_MG_oBack = new size_t *[m_levels + 1];
        // m_MG_oBottom = new size_t *[m_levels + 1];
        // m_MG_oTop = new size_t *[m_levels + 1];
        // m_MG_oLeft = new size_t *[m_levels + 1];
        // m_MG_oRight = new size_t *[m_levels + 1];

        // start index of each obstacle in level joined list (SliceZ = Front/Back, SliceY = Bottom/Top, SliceX = Left/Right)
        m_size_MG_oFront_level = new size_t[(m_levels + 1) * m_numberOfObstacles + 1];
        m_size_MG_oBack_level = new size_t[(m_levels + 1) * m_numberOfObstacles + 1];
        m_size_MG_oTop_level = new size_t[(m_levels + 1) * m_numberOfObstacles + 1];
        m_size_MG_oBottom_level = new size_t[(m_levels + 1) * m_numberOfObstacles + 1];
        m_size_MG_oLeft_level = new size_t[(m_levels + 1) * m_numberOfObstacles + 1];
        m_size_MG_oRight_level = new size_t[(m_levels + 1) * m_numberOfObstacles + 1];

        // size of obstacle level 0 / start index of level 1
        m_size_MG_oFront_level[0] = 0;
        m_size_MG_oBack_level[0] = 0;
        m_size_MG_oTop_level[0] = 0;
        m_size_MG_oBottom_level[0] = 0;
        m_size_MG_oLeft_level[0] = 0;
        m_size_MG_oRight_level[0] = 0;
        *(m_size_MG_oList_level) = 0;
        // TODO notwendig?
        for (size_t i = 0; i < m_numberOfObstacles; i++) {
            m_size_MG_oFront_level[i + 1] = obstacle_list[i]->getSize_obstacleFront() + m_size_MG_oFront_level[i];
            m_size_MG_oBack_level[i + 1] = obstacle_list[i]->getSize_obstacleBack() + m_size_MG_oBack_level[i];
            m_size_MG_oBottom_level[i + 1] = obstacle_list[i]->getSize_obstacleBottom() + m_size_MG_oBottom_level[i];
            m_size_MG_oTop_level[i + 1] = obstacle_list[i]->getSize_obstacleTop() + m_size_MG_oTop_level[i];
            m_size_MG_oLeft_level[i + 1] = obstacle_list[i]->getSize_obstacleLeft() + m_size_MG_oLeft_level[i];
            m_size_MG_oRight_level[i + 1] = obstacle_list[i]->getSize_obstacleRight() + m_size_MG_oRight_level[i];
            *(m_size_MG_oList_level) += obstacle_list[i]->getSize_obstacleList();
        }
        m_bdc_obstacle = bdc_obstacles;
    }

    init();
    addMGLists();
    sendListsToGPU();
#ifndef BENCHMARKING
    // print();
    control();
#endif
    m_data_boundary_patches_joined = new size_t *[number_of_patches];
    m_data_boundary_patches_joined[Patch::FRONT] = m_data_MG_bFront_level_joined;
    m_data_boundary_patches_joined[Patch::BACK] = m_data_MG_bBack_level_joined;
    m_data_boundary_patches_joined[Patch::BOTTOM] = m_data_MG_bBottom_level_joined;
    m_data_boundary_patches_joined[Patch::TOP] = m_data_MG_bTop_level_joined;
    m_data_boundary_patches_joined[Patch::LEFT] = m_data_MG_bLeft_level_joined;
    m_data_boundary_patches_joined[Patch::RIGHT] = m_data_MG_bRight_level_joined;

    if (m_numberOfObstacles > 0) {
        m_data_obstacles_patches_joined = new size_t *[number_of_patches];
        m_data_obstacles_patches_joined[Patch::FRONT] = m_data_MG_oFront_level_joined;
        m_data_obstacles_patches_joined[Patch::BACK] = m_data_MG_oBack_level_joined;
        m_data_obstacles_patches_joined[Patch::BOTTOM] = m_data_MG_oBottom_level_joined;
        m_data_obstacles_patches_joined[Patch::TOP] = m_data_MG_oTop_level_joined;
        m_data_obstacles_patches_joined[Patch::LEFT] = m_data_MG_oLeft_level_joined;
        m_data_obstacles_patches_joined[Patch::RIGHT] = m_data_MG_oRight_level_joined;
    }
}

//======================================== Init ====================================
// ***************************************************************************************
/// \brief  Initialize member variables (arrays)
// ***************************************************************************************
void Multigrid::init() {
    m_levels = Domain::getInstance()->get_levels();  // multigrid level, 0 otherwise

    // list of domain boundary for each level
    m_MG_boundaryList = new Boundary *[m_levels + 1];

    // start index of each level (difference equals size of respective element)
    m_size_MG_iList_level = new size_t[m_levels + 2];
    m_size_MG_bList_level = new size_t[m_levels + 2];
    m_size_MG_bSliceZ_level = new size_t[m_levels + 2];
    m_size_MG_bSliceY_level = new size_t[m_levels + 2];
    m_size_MG_bSliceX_level = new size_t[m_levels + 2];

    // start index of first element = 0
    *(m_size_MG_iList_level) = 0;
    *(m_size_MG_bList_level) = 0;

    *(m_size_MG_bSliceZ_level) = 0;
    *(m_size_MG_bSliceY_level) = 0;
    *(m_size_MG_bSliceX_level) = 0;
}

Multigrid::~Multigrid() {
    delete[] m_data_boundary_patches_joined;
    for (size_t level = 0; level < m_levels + 1; level++) {
        if (m_numberOfSurfaces > 0) {
            Surface **surface_level = *(m_MG_surfaceList + level);
            for (size_t surface = 0; surface < m_numberOfSurfaces; surface++) {
                delete (*(surface_level + surface));
            }
            delete[] surface_level;
            delete (*(m_MG_sList + level));
        }
        if (m_numberOfObstacles > 0) {
            Obstacle **obstacle_level = *(m_MG_obstacleList);
            for (size_t obstacle = 0; obstacle < m_numberOfObstacles; obstacle++) {
                delete (*(obstacle_level + obstacle));
            }
            delete[] obstacle_level;

        //    delete (*(m_MG_oList + level));
        //    delete (*(m_MG_oFront + level));
        //    delete (*(m_MG_oBack + level));
        //    delete (*(m_MG_oBottom + level));
        //    delete (*(m_MG_oTop + level));
        //    delete (*(m_MG_oLeft + level));
        //    delete (*(m_MG_oRight + level));
        }
        delete (*(m_MG_boundaryList + level));
    }
    delete[] m_MG_boundaryList;

    size_t size_iList = getLen_iList_joined();
    size_t size_bList = getLen_bList_joined();
#pragma acc exit data delete(m_data_MG_iList_level_joined[:size_iList])
#pragma acc exit data delete(m_data_MG_bList_level_joined[:size_bList])
    delete[] m_data_MG_inner_list_level_joined;
    delete[] m_data_MG_bList_level_joined;

    size_t size_bSliceZ = getLen_bSliceZ_joined();
    size_t size_bSliceY = getLen_bSliceY_joined();
    size_t size_bSliceX = getLen_bSliceX_joined();
#pragma acc exit data delete(m_data_MG_bFront_level_joined[:size_bSliceZ])
#pragma acc exit data delete(m_data_MG_bBack_level_joined[:size_bSliceZ])
#pragma acc exit data delete(m_data_MG_bTop_level_joined[:size_bSliceY])
#pragma acc exit data delete(m_data_MG_bBottom_level_joined[:size_bSliceY])
#pragma acc exit data delete(m_data_MG_bLeft_level_joined[:size_bSliceX])
#pragma acc exit data delete(m_data_MG_bRight_level_joined[:size_bSliceX])
    delete[] m_data_MG_bFront_level_joined;
    delete[] m_data_MG_bBack_level_joined;
    delete[] m_data_MG_bTop_level_joined;
    delete[] m_data_MG_bBottom_level_joined;
    delete[] m_data_MG_bLeft_level_joined;
    delete[] m_data_MG_bRight_level_joined;

    delete[] m_size_MG_iList_level;
    delete[] m_size_MG_bList_level;

    delete[] m_size_MG_bSliceZ_level;
    delete[] m_size_MG_bSliceY_level;
    delete[] m_size_MG_bSliceX_level;

    if (m_numberOfSurfaces > 0) {
        size_t size_sList = getLen_sList_joined();
#pragma acc exit data delete(m_data_MG_sList_level_joined[:size_sList])
        delete[] m_MG_surfaceList;
        delete[] m_MG_sList;
        delete[] m_size_MG_sList_level;
        delete[] m_data_MG_sList_level_joined;
    }
    if (m_numberOfObstacles > 0) {
        delete[] m_MG_obstacleList;
        delete[] m_MG_oList;
     //   delete[] m_MG_oFront;
     //   delete[] m_MG_oBack;
     //   delete[] m_MG_oBottom;
     //   delete[] m_MG_oTop;
     //   delete[] m_MG_oLeft;
     //   delete[] m_MG_oRight;

//        size_t size_oFront = getLen_oFront_joined();
//        size_t size_oBack = getLen_oBack_joined();
//        size_t size_oBottom = getLen_oBottom_joined();
//        size_t size_oTop = getLen_oTop_joined();
//        size_t size_oLeft = getLen_oLeft_joined();
//        size_t size_oRight = getLen_oRight_joined();
// #pragma acc exit data delete(m_data_MG_oFront_level_joined[:size_oFront])
// #pragma acc exit data delete(m_data_MG_oBack_level_joined[:size_oBack])
// #pragma acc exit data delete(m_data_MG_oTop_level_joined[:size_oTop])
// #pragma acc exit data delete(m_data_MG_oBottom_level_joined[:size_oBottom])
// #pragma acc exit data delete(m_data_MG_oLeft_level_joined[:size_oLeft])
// #pragma acc exit data delete(m_data_MG_oRight_level_joined[:size_oRight])
        delete[] m_data_MG_oFront_level_joined;
        delete[] m_data_MG_oBack_level_joined;
        delete[] m_data_MG_oTop_level_joined;
        delete[] m_data_MG_oBottom_level_joined;
        delete[] m_data_MG_oLeft_level_joined;
        delete[] m_data_MG_oRight_level_joined;

        delete[] m_size_MG_oFront_level;
        delete[] m_size_MG_oBack_level;
        delete[] m_size_MG_oTop_level;
        delete[] m_size_MG_oBottom_level;
        delete[] m_size_MG_oLeft_level;
        delete[] m_size_MG_oRight_level;
    }
}

//======================================== Control ====================================
// ***************************************************************************************
/// \brief  Units test emergency solution
// ***************************************************************************************
void Multigrid::control() {
    // TODO control
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
        size_t bLen = (static_cast<Boundary *>(*(m_MG_boundaryList + level)))->get_size_inner_list();
        size_t cLen = getLastIndex_iList(level) - getFirstIndex_iList(level) + 1;
        if (cLen != bLen) {
            size_t control = (domain->get_nx(level) - 2) * (domain->get_ny(level) - 2) * (domain->get_nz(level) - 2);
            message += "length calculated by first and last index of iList does not equals size of innerList of Boundary object " + std::to_string(cLen) + "|" + std::to_string(bLen) + " control: " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundaryList + level)))->get_size_boundary_list();
        cLen = getLastIndex_bList(level) - getFirstIndex_bList(level) + 1;
        if (cLen != bLen) {
            size_t control = (domain->get_nx(level) * domain->get_ny(level) * 2) + (domain->get_nx(level) * (domain->get_nz(level) - 2)) * 2 + ((domain->get_ny(level) - 2) * (domain->get_nz(level) - 2)) * 2;
            message += "length calculated by first and last index of bList does not equals size of boundaryList of Boundary object " + std::to_string(cLen) + "|" + std::to_string(bLen) + " control: " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundaryList + level)))->get_size_boundary_left();
        cLen = getLastIndex_bSliceX(level) - getFirstIndex_bSliceX(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_ny(level) * domain->get_nz(level);
            message += "length calculated by first and last index of bSliceX does not equals size of bSliceX of bLeft " + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundaryList + level)))->get_size_boundary_right();
        cLen = getLastIndex_bSliceX(level) - getFirstIndex_bSliceX(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_ny(level) * domain->get_nz(level);
            message += "length calculated by first and last index of bSliceX does not equals size of bSliceX of bRight " + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundaryList + level)))->get_size_boundary_bottom();
        cLen = getLastIndex_bSliceY(level) - getFirstIndex_bSliceY(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_nx(level) * domain->get_nz(level);
            message += "length calculated by first and last index of bSliceY does not equals size of bSliceY of bBottom " + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundaryList + level)))->get_size_boundary_top();
        cLen = getLastIndex_bSliceY(level) - getFirstIndex_bSliceY(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_nx(level) * domain->get_nz(level);
            message += "length calculated by first and last index of bSliceY does not equals size of bSliceY of bTop " + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundaryList + level)))->get_size_boundary_front();
        cLen = getLastIndex_bSliceZ(level) - getFirstIndex_bSliceZ(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_nx(level) * domain->get_ny(level);
            message += "length calculated by first and last index of bSliceZ does not equals size of bSliceZ of bFront " + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }
        bLen = (static_cast<Boundary *>(*(m_MG_boundaryList + level)))->get_size_boundary_back();
        cLen = getLastIndex_bSliceZ(level) - getFirstIndex_bSliceZ(level) + 1;
        if (cLen != bLen) {
            size_t control = domain->get_nx(level) * domain->get_ny(level);
            message += "length calculated by first and last index of bSliceZ does not equals size of bSliceZ of bBack " + std::to_string(cLen) + "|" + std::to_string(bLen) + " control " + std::to_string(control) + "\n";
        }

        size_t csize_inner = (domain->get_nx(level) - 2) * (domain->get_ny(level) - 2) * (domain->get_nz(level) - 2) - getSize_oList(level);
        size_t bsize_inner = get_size_inner_list(level);
        if (csize_inner != bsize_inner) {
            message += "get_size_inner_list(level) does not equal (nx-2)*(ny-2)*(nz-2) " + std::to_string(bsize_inner) + "|" + std::to_string(csize_inner) + "\n";
        }
        size_t cindex_inner_start = getInnerList_level_joined_start(level);
        size_t cindex_inner_end = getInnerList_level_joined_end(level);
        if (cindex_inner_end - cindex_inner_start + 1 != bsize_inner) {
            message += "get_size_inner_list(level) does not equal the difference between start and end " + std::to_string(cindex_inner_start) + "|" + std::to_string(cindex_inner_end) + "\n";
        }
        size_t bsize_boundary = nx * ny * nz - (nx - 2) * (ny - 2) * (nz - 2);
        size_t csize_boundary = get_size_boundary_ist(level);
        if (csize_boundary != bsize_boundary) {
            message += "get_size_boundary_list(level) does not equal size-(nx-2)*(ny-2)*(nz-2) " + std::to_string(bsize_boundary) + "|" + std::to_string(csize_boundary) + "\n";
        }

        if (!message.empty()) {
            message = "For Level " + std::to_string(level) + "\nsize control\n" + message;
        }
    }
    for (size_t level = 0; level < m_levels + 1; level++) {
        Boundary *b = *(m_MG_boundaryList + level);
        b->control(getSize_oList(level));
    }
    {
        size_t bsize_inner = 0;
        for (size_t level = 0; level < m_levels + 1; level++) {
            bsize_inner += get_size_inner_list(level);
        }
        size_t csize_inner = getSize_innerList_level_joined();
        if (bsize_inner != csize_inner) {
            message += "get_size_inner_list_level_joined does not equal the sum of each inner list " + std::to_string(bsize_inner) + "|" + std::to_string(csize_inner) + "\n";
        }
    }

    {
        std::vector<size_t> v(m_data_MG_oList_zero_joined, m_data_MG_oList_zero_joined + get_size_obstacle_list());
        std::sort(v.begin(), v.end());
        std::vector<size_t>::iterator it = std::unique(v.begin(), v.begin() + get_size_obstacle_list());
        v.resize(std::distance(v.begin(), it));
        if (v.size() != get_size_obstacle_list()) {
            message += "obstacles are not allowed to overlap in level 0! difference: " + std::to_string(get_size_obstacle_list() - v.size()) + "\n";

            for (size_t i = 0; i < m_numberOfObstacles; i++) {
                Obstacle *o1 = m_MG_obstacleList[0][i];
                size_t  size1 = o1->getSize_obstacleList();
                for (size_t j = i+1; j < m_numberOfObstacles; j++) {
                    Obstacle *o2 = m_MG_obstacleList[0][j];
                    size_t size2 = o2->getSize_obstacleList();
                    std::vector<size_t> tmp = Utility::mergeSortedListsToUniqueList(o1->getObstacleList(), size1, o2->getObstacleList(), size2);
                    if (tmp.size() != size1 + size2) {
                        message += "Obstacles " + std::to_string(i) + " and " + std::to_string(j) + "are overlapping\n";
                    }
                }
            }
        }
    }
    if (!message.empty()) {
        message = "################ MULTIGRID CONTROL ################\n" + message + "---------------- MULTIGRID CONTROL END ----------------";
#ifndef BENCHMARKING
        m_logger->warn(message);
#endif
    }
}

//======================================== Print ====================================
// ***************************************************************************************
/// \brief  Print multigrid infos
// ***************************************************************************************
void Multigrid::print() {
#ifndef BENCHMARKING
    m_logger->debug("################ MULTIGRID ################");
    m_logger->debug("Number of Obstacles: {}, Number of Surfaces: {}",
            m_numberOfObstacles, m_numberOfSurfaces);
    m_logger->debug("Levels: {}", m_levels);

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} iList starts at index {} and ends with index {}",
                level, getFirstIndex_iList(level), getLastIndex_iList(level));
        m_logger->debug("and the corresponding indices at this position: {}|{}",
                *(m_data_MG_inner_list_level_joined + getFirstIndex_iList(level)),
                *(m_data_MG_inner_list_level_joined + getLastIndex_iList(level)));
    }
    m_logger->debug("Total length of iList: {}", getLen_iList_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} bList starts at index {} and ends with index {}",
                level, getFirstIndex_bList(level), getLastIndex_bList(level));
        m_logger->debug("and the corresponding indices at this position: {} | {}",
                *(m_data_MG_bList_level_joined + getFirstIndex_bList(level)),
                *(m_data_MG_bList_level_joined + getLastIndex_bList(level)));
    }
    m_logger->debug("Total length of bList: {}", getLen_bList_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} bSliceZ starts at index {} and ends with index {}",
                level, getFirstIndex_bSliceZ(level), getLastIndex_bSliceZ(level));
        m_logger->debug(" and the corresponding indices at this position for FRONT: {} | {}",
                *(m_data_MG_bFront_level_joined + getFirstIndex_bSliceZ(level)),
                *(m_data_MG_bFront_level_joined + getLastIndex_bSliceZ(level)));
        m_logger->debug(" and the corresponding indices at this position for BACK : {} | {}",
                *(m_data_MG_bBack_level_joined + getFirstIndex_bSliceZ(level)),
                *(m_data_MG_bBack_level_joined + getLastIndex_bSliceZ(level)));
    }
    m_logger->debug("Total length of bSliceZ: {}", getLen_bSliceZ_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} bSliceY starts at index {} and ends with index {}",
            level, getFirstIndex_bSliceY(level), getLastIndex_bSliceY(level));
        m_logger->debug(" and the corresponding indices at this position for BOTTOM: {} | {}",
                *(m_data_MG_bBottom_level_joined + getFirstIndex_bSliceY(level)),
                *(m_data_MG_bBottom_level_joined + getLastIndex_bSliceY(level)));
        m_logger->debug(" and the corresponding indices at this position for TOP   : {} | {}",
                *(m_data_MG_bTop_level_joined + getFirstIndex_bSliceY(level)),
                *(m_data_MG_bTop_level_joined + getLastIndex_bSliceY(level)));
    }
    m_logger->debug("Total length of bSliceY: {}", getLen_bSliceY_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} bSliceX starts at index {} and ends with index {}",
                level, getFirstIndex_bSliceX(level), getLastIndex_bSliceX(level));
        m_logger->debug(" and the corresponding indices at this position for LEFT : {} | {}",
                *(m_data_MG_bLeft_level_joined + getFirstIndex_bSliceX(level)),
                *(m_data_MG_bLeft_level_joined + getLastIndex_bSliceX(level)));
        m_logger->debug(" and the corresponding indices at this position for RIGHT: {} | {}",
                *(m_data_MG_bRight_level_joined + getFirstIndex_bSliceX(level)),
                *(m_data_MG_bRight_level_joined + getLastIndex_bSliceX(level)));
    }
    m_logger->debug("Total length of bSliceX: {}", getLen_bSliceX_joined());

    for (size_t level = 0; level < m_levels; level++) {
        m_logger->debug("For Level {} obstacleList has {} elements",
                level, getSize_oList(level));
    }

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oFront starts at index {} and ends with index {} with length {}",
                level ,
                getFirstIndex_oFront(level, 0),
                getLastIndex_oFront(level, m_numberOfObstacles - 1) ,
                getLen_oFront(level));
    }
    m_logger->debug("Total length of oFront: ", getLen_oFront_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oBack starts at index {} and ends with index {} with length {}",
                level,
                getFirstIndex_oBack(level, 0),
                getLastIndex_oBack(level, m_numberOfObstacles - 1),
                getLen_oBack(level));
    }
    m_logger->debug("Total length of oBack: {}", getLen_oBack_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oBottom starts at index {} and ends with index {} with length {}",
                level,
                getFirstIndex_oBottom(level, 0),
                getLastIndex_oBottom(level, m_numberOfObstacles - 1),
                getLen_oBottom(level));
    }
    m_logger->debug("Total length of oBottom: {}", getLen_oBottom_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oTop starts at index {} and ends with index {} with length {}",
                level,
                getFirstIndex_oTop(level, 0),
                getLastIndex_oTop(level, m_numberOfObstacles - 1),
                getLen_oTop(level));
    }
    m_logger->debug("Total length of oTop: {}", getLen_oTop_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oLeft starts at index {} and ends with index {} with length {}",
                level,
                getFirstIndex_oLeft(level, 0),
                getLastIndex_oLeft(level, m_numberOfObstacles - 1),
                getLen_oLeft(level));
    }
    m_logger->debug("Total length of oLeft: {}", getLen_oLeft_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} oRight starts at index {} and ends with index {} with length {}",
                level,
                getFirstIndex_oRight(level, 0),
                getLastIndex_oRight(level, m_numberOfObstacles - 1),
                getLen_oRight(level));
    }
    m_logger->debug("Total length of oRight: {}", getLen_oRight_joined());

    for (size_t level = 0; level < m_levels + 1; level++) {
        m_logger->debug("For Level {} sList starts at index {} and ends with index {}",
                level,
                getFirstIndex_sList(level),
                getLastIndex_sList(level));
    }
    m_logger->debug("Total length of sList: {}", getLen_sList_joined());
    m_logger->debug("---------------- MULTIGRID END ----------------");
#endif
}

// ================================= Add MG lists ==========================================
// ***************************************************************************************
/// \brief  adds lists (outer, inner, surfaces, obstacles) in case of grid restriction (dominant)
// ***************************************************************************************
void Multigrid::addMGLists() {
    // create Boundary object of level 0
    Boundary *b;
    if (m_numberOfObstacles > 0) {
        calcObstacles(*(m_MG_obstacleList));
        b = new Boundary(*(m_MG_obstacleList), m_numberOfObstacles, getSize_oList(0));
    } else {
        b = new Boundary();
    }
    // set size of respective lists
    m_size_MG_iList_level[1] = b->get_size_inner_list();
    m_size_MG_bList_level[1] = b->get_size_boundary_list();
    m_size_MG_bSliceZ_level[1] = b->get_size_boundary_front();
    m_size_MG_bSliceY_level[1] = b->get_size_boundary_top();
    m_size_MG_bSliceX_level[1] = b->get_size_boundary_left();
    // save boundary object in multigrid list
    *(m_MG_boundaryList) = b;

    if (m_numberOfSurfaces > 0) {
        calcSurfaces(*(m_MG_surfaceList));
    }

    // create boundary object, surfaces and obstacles for each multigrid level
    for (size_t level = 1; level < m_levels + 1; level++) {
        surfaceDominantRestriction(level);
        Obstacle **obstacleList = obstacleDominantRestriction(level);

        Boundary *boundary;
        if (m_numberOfObstacles > 0) {
            boundary = new Boundary(obstacleList, m_numberOfObstacles, getSize_oList(level), level);
        } else {
            boundary = new Boundary(level);
        }
        *(m_MG_boundaryList + level) = boundary;

        m_size_MG_iList_level[level + 1] = m_size_MG_iList_level[level] + boundary->get_size_inner_list();
        m_size_MG_bList_level[level + 1] = m_size_MG_bList_level[level] + boundary->get_size_boundary_list();
        m_size_MG_bSliceZ_level[level + 1] = m_size_MG_bSliceZ_level[level] + boundary->get_size_boundary_front();
        m_size_MG_bSliceY_level[level + 1] = m_size_MG_bSliceY_level[level] + boundary->get_size_boundary_top();
        m_size_MG_bSliceX_level[level + 1] = m_size_MG_bSliceX_level[level] + boundary->get_size_boundary_left();
    }
}

// ================================= Calc obstacles ==========================================
// ***************************************************************************************
/// \brief  create obstacles of level 0, no check if obstacles are overlapping
/// \param obstacleList List of obstacle objects
// ***************************************************************************************
void Multigrid::calcObstacles(Obstacle **obstacleList) {
    if (m_numberOfObstacles > 0) {
        size_t level = 0;
        size_t o_size = getSize_oList(level);
        size_t *oList = new size_t[o_size];
        size_t counter = 0;
        for (size_t o = 0; o < m_numberOfObstacles; o++) {
            Obstacle *obstacle_tmp = obstacleList[o];
            size_t len_all = obstacle_tmp->getSize_obstacleList();
            for (size_t i = 0; i < len_all; i++) {
                *(oList + counter) = obstacle_tmp->getObstacleList()[i];
                counter++;
            }
        }
        *(m_MG_oList) = oList;
    }
}

// ================================= Calc surfaces ==========================================
// ***************************************************************************************
/// \brief  Create surfaces of level 0
/// \param surfaceList List of surface objects
// ***************************************************************************************
void Multigrid::calcSurfaces(Surface **surfaceList) {
    if (m_numberOfSurfaces > 0) {
        size_t *sList = new size_t[getFirstIndex_sList(1)];
        size_t counter_s = 0;
        for (size_t s = 0; s < m_numberOfSurfaces; s++) {
            Surface *surface_tmp = surfaceList[s];
            for (size_t i = 0; i < surface_tmp->getSize_surfaceList(); i++) {
                *(sList + counter_s) = surface_tmp->getSurfaceList()[i];
                counter_s++;
            }
            *(m_size_MG_sList_level + s + 1) = counter_s;
        }
        *(m_MG_sList) = sList;
    }
}

// ================================= Surface dominant restriction ==========================================
// ***************************************************************************************
/// \brief  Calculates coarse level of multigrid of surfaces (dominant restriction)
/// \param level Multigrid level
// ***************************************************************************************
void Multigrid::surfaceDominantRestriction(size_t level) {
// SURFACES
    if (m_numberOfSurfaces > 0) {
// add index to m_surface_list if any of l-1 indices building the l index was a surface (dominant restriction)
        Domain *domain = Domain::getInstance();
        size_t Nx = domain->get_Nx(level);
        size_t Ny = domain->get_Ny(level);

        Surface **surfaceList_fine = *(m_MG_surfaceList + (level - 1));
        Surface **surfaceList_coarse = new Surface *[m_numberOfSurfaces];
        *(m_MG_surfaceList + level) = surfaceList_coarse;
        // loop through surfaces
        for (size_t surfaceID = 0; surfaceID < m_numberOfSurfaces; ++surfaceID) {
            Surface *surface_fine = surfaceList_fine[surfaceID];

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

            Surface *surface_coarse = new Surface(surfaceID, startIndex_coarse, stride_x_coarse, stride_y_coarse, stride_z_coarse, level);
            *(surfaceList_coarse + surfaceID) = surface_coarse;
            size_t index = level * m_numberOfSurfaces + surfaceID + 1;
            *(m_size_MG_sList_level + index) = *(m_size_MG_sList_level + index - 1) + surface_coarse->getSize_surfaceList();
#ifndef BENCHMARKING
            m_logger->debug("control multigrid surface index: {} {}",
                    *(m_size_MG_sList_level + index),
                    surface_coarse->getSize_surfaceList());
#endif
        }  // end surface id loop
    }
}

// ================================= Obstacle dominant restriction ==========================================
// ***************************************************************************************
/// \brief  Calculates coarse level of multigrid of obstacles (dominant restriction)
/// \param level Multigrid level
/// \return Obstacle** list of created obstacles
// ***************************************************************************************
Obstacle** Multigrid::obstacleDominantRestriction(size_t level) {
// OBSTACLES
    // add index to m_oLists if any of l-1 indices building the l index was an obstacle (dominant restriction)
    if (m_numberOfObstacles == 0) {
        return nullptr;
    }
    // TODO define lists with obstacle size
    Domain *domain = Domain::getInstance();
    Obstacle **obstacleList_fine = *(m_MG_obstacleList + (level - 1));
    Obstacle **obstacleList_coarse = new Obstacle *[m_numberOfObstacles];
    *(m_MG_obstacleList + level) = obstacleList_coarse;
    *(m_size_MG_oList_level + level) = 0;
    for (size_t id = 0; id < m_numberOfObstacles; id++) {
        Obstacle *obstacle_fine = obstacleList_fine[id];
        size_t i1_fine = obstacle_fine->getCoordinates_i1();
        size_t j1_fine = obstacle_fine->getCoordinates_j1();
        size_t k1_fine = obstacle_fine->getCoordinates_k1();
        size_t i2_fine = obstacle_fine->getCoordinates_i2();
        size_t j2_fine = obstacle_fine->getCoordinates_j2();
        size_t k2_fine = obstacle_fine->getCoordinates_k2();

        size_t i1_coarse = (i1_fine + 1) / 2;
        size_t j1_coarse = (j1_fine + 1) / 2;
        size_t k1_coarse = (k1_fine + 1) / 2;
        size_t i2_coarse = (i2_fine + 1) / 2;
        size_t j2_coarse = (j2_fine + 1) / 2;
        size_t k2_coarse = (k2_fine + 1) / 2;

        // TODO exit?
#ifndef BENCHMARKING
        if (i2_fine - i1_fine + 1< domain->get_nx(level - 1) - 2 && i2_coarse - i1_coarse + 1 >= domain->get_nx(level) - 2) {
            m_logger->warn("Be cautious! Obstacle fills up inner cells in x-direction at level {}", level);
        }
        if (j2_fine - j1_fine +1< domain->get_ny(level - 1) - 2 && j2_coarse - j1_coarse + 1 >= domain->get_ny(level) - 2) {
            m_logger->warn("Be cautious! Obstacle fills up inner cells in y-direction at level {}", level);
        }
        if (k2_fine - k1_fine +1< domain->get_nz(level - 1) - 2 && k2_coarse - k1_coarse + 1 >= domain->get_nz(level) - 2) {
            m_logger->warn("Be cautious! Obstacle fills up inner cells in z-direction at level {}", level);
        }
#endif
        // for (size_t c = 0; c < id; c++) {
        //     controlObstacleOverlap(obstacleList_coarse[c], &i1_coarse, &i2_coarse, &j1_coarse, &j2_coarse, &k1_coarse, &k2_coarse);
        // }

        Obstacle *obstacle_coarse = new Obstacle(i1_coarse, j1_coarse, k1_coarse, i2_coarse, j2_coarse, k2_coarse, level, obstacle_fine->get_name());
        *(obstacleList_coarse + id) = obstacle_coarse;

        size_t index = level * m_numberOfObstacles + id + 1;
        m_size_MG_oFront_level[index] = obstacle_coarse->getSize_obstacleFront() + m_size_MG_oFront_level[index - 1];
        m_size_MG_oBack_level[index] = obstacle_coarse->getSize_obstacleBack() + m_size_MG_oBack_level[index - 1];
        m_size_MG_oBottom_level[index] = obstacle_coarse->getSize_obstacleBottom() + m_size_MG_oBottom_level[index - 1];
        m_size_MG_oTop_level[index] = obstacle_coarse->getSize_obstacleTop() + m_size_MG_oTop_level[index - 1];
        m_size_MG_oLeft_level[index] = obstacle_coarse->getSize_obstacleLeft() + m_size_MG_oLeft_level[index - 1];
        m_size_MG_oRight_level[index] = obstacle_coarse->getSize_obstacleRight() + m_size_MG_oRight_level[index - 1];
        *(m_size_MG_oList_level + level) += obstacle_coarse->getSize_obstacleList();
    } //end obstacle id loop

    size_t *list = obstacleList_coarse[0]->getObstacleList();
    size_t size = obstacleList_coarse[0]->getSize_obstacleList();
    std::vector<size_t> data;
    data.assign(list, list+size);

    for (size_t o = 1; o < m_numberOfObstacles; o++) {
        Obstacle *obstacle = obstacleList_coarse[o];
        data = Utility::mergeSortedListsToUniqueList(list, size, obstacle->getObstacleList(), obstacle->getSize_obstacleList());
        list = data.data();
        size = data.size();
    }

    size_t *oList_tmp = new size_t[size];
    std::copy(&data[0], &data[size], oList_tmp);
    *(m_size_MG_oList_level + level) = size;
    *(m_MG_oList + level) = oList_tmp;
    return obstacleList_coarse;
}

// ================================= Send lists to GPU ==========================================
// ***************************************************************************************
/// \brief  control and correct when obstacles are overlapping caused by the dominant restriction
// ***************************************************************************************
void Multigrid::controlObstacleOverlap(Obstacle *o, size_t *i1, size_t *i2, size_t *j1, size_t *j2, size_t *k1, size_t *k2) {
    size_t o_i1 = o->getCoordinates_i1();
    size_t o_i2 = o->getCoordinates_i2();
    size_t o_j1 = o->getCoordinates_j1();
    size_t o_j2 = o->getCoordinates_j2();
    size_t o_k1 = o->getCoordinates_k1();
    size_t o_k2 = o->getCoordinates_k2();

    if (o_i1 <= *i1 && *i1 <= o_i2) {
        *i1 = o_i2 + 1;
    }
    if (o_i1 <= *i2 && *i2 <= o_i2) {
        *i2 = o_i1 - 1;
    }
    if (o_j1 <= *j1 && *j1 <= o_j2) {
        *j1 = o_j2 + 1;
    }
    if (o_j1 <= *j2 && *j2 <= o_j2) {
        *j2 = o_j1 - 1;
    }
    if (o_k1 <= *k1 && *k1 <= o_k2) {
        *k1 = o_k2 + 1;
    }
    if (o_k1 <= *k2 && *k2 <= o_k2) {
        *k2 = o_k1 - 1;
    }
}

// ================================= Send lists to GPU ==========================================
// ***************************************************************************************
/// \brief  create joined list and send them to GPU
// ***************************************************************************************
void Multigrid::sendListsToGPU() {
    sendSurfaceListsToGPU();
    sendBoundaryListsToGPU();
    sendObstacleListsToGPU();
}

// ================================= Send boundary lists to GPU ====================================
// ***************************************************************************************
/// \brief  create boundary joined list and send them to GPU
// ***************************************************************************************
void Multigrid::sendBoundaryListsToGPU() {
    size_t size_iList = getLen_iList_joined();
    size_t size_bList = getLen_bList_joined();
    size_t size_bSliceZ = getLen_bSliceZ_joined();
    size_t size_bSliceY = getLen_bSliceY_joined();
    size_t size_bSliceX = getLen_bSliceX_joined();

    m_data_MG_inner_list_level_joined = new size_t[size_iList];
    m_data_MG_bList_level_joined = new size_t[size_bList];
    m_data_MG_bFront_level_joined = new size_t[size_bSliceZ];
    m_data_MG_bBack_level_joined = new size_t[size_bSliceZ];
    m_data_MG_bTop_level_joined = new size_t[size_bSliceY];
    m_data_MG_bBottom_level_joined = new size_t[size_bSliceY];
    m_data_MG_bLeft_level_joined = new size_t[size_bSliceX];
    m_data_MG_bRight_level_joined = new size_t[size_bSliceX];

    size_t counter_iList = 0;
    size_t counter_bList = 0;

    size_t counter_bSliceZ = 0;
    size_t counter_bSliceY = 0;
    size_t counter_bSliceX = 0;

    for (size_t level = 0; level < m_levels + 1; level++) {
        Boundary *boundary = *(m_MG_boundaryList + level);
        for (size_t i = 0; i < boundary->get_size_inner_list(); i++) {
            *(m_data_MG_inner_list_level_joined + counter_iList) = boundary->get_inner_list()[i];
            counter_iList++;
        }
        for (size_t i = 0; i < boundary->get_size_boundary_list(); i++) {
            *(m_data_MG_bList_level_joined + counter_bList) = boundary->get_boundary_list()[i];
            counter_bList++;
        }
        for (size_t i = 0; i < boundary->get_size_boundary_front(); i++) {
            *(m_data_MG_bFront_level_joined + counter_bSliceZ) = boundary->get_boundary_front()[i];
            *(m_data_MG_bBack_level_joined + counter_bSliceZ) = boundary->get_boundary_back()[i];
            counter_bSliceZ++;
        }
        for (size_t i = 0; i < boundary->get_size_boundary_top(); i++) {
            *(m_data_MG_bTop_level_joined + counter_bSliceY) = boundary->get_boundary_top()[i];
            *(m_data_MG_bBottom_level_joined + counter_bSliceY) = boundary->get_boundary_bottom()[i];
            counter_bSliceY++;
        }
        for (size_t i = 0; i < boundary->get_size_boundary_left(); i++) {
            *(m_data_MG_bLeft_level_joined + counter_bSliceX) = boundary->get_boundary_left()[i];
            *(m_data_MG_bRight_level_joined + counter_bSliceX) = boundary->get_boundary_right()[i];
            counter_bSliceX++;
        }
    }
#pragma acc enter data copyin(m_data_MG_iList_level_joined[:size_iList])
#pragma acc enter data copyin(m_data_MG_bList_level_joined[:size_bList])
#pragma acc enter data copyin(m_data_MG_bFront_level_joined[:size_bSliceZ])
#pragma acc enter data copyin(m_data_MG_bBack_level_joined[:size_bSliceZ])
#pragma acc enter data copyin(m_data_MG_bTop_level_joined[:size_bSliceY])
#pragma acc enter data copyin(m_data_MG_bBottom_level_joined[:size_bSliceY])
#pragma acc enter data copyin(m_data_MG_bLeft_level_joined[:size_bSliceX])
#pragma acc enter data copyin(m_data_MG_bRight_level_joined[:size_bSliceX])
}

// ================================= Send surface list to GPU ==========================================
// ***************************************************************************************
/// \brief  create surface joined list and send it to GPU
// ***************************************************************************************
void Multigrid::sendSurfaceListsToGPU() {
    size_t counter_sList = 0;

    size_t size_sList = 0;
    if (m_numberOfSurfaces > 0) {
        size_sList = getLen_sList_joined();
        m_data_MG_sList_level_joined = new size_t[size_sList];
#ifndef BENCHMARKING
        m_logger->debug("control sendMGListsToGPU size surface {}", size_sList);
#endif
        for (size_t level = 0; level < m_levels + 1; level++) {
            Surface **surfaceList = *(m_MG_surfaceList + level);
            for (size_t s = 0; s < m_numberOfSurfaces; s++) {
                Surface *surface = *(surfaceList + s);
                for (size_t i = 0; i < surface->getSize_surfaceList(); i++) {
                    *(m_data_MG_sList_level_joined + counter_sList) = surface->getSurfaceList()[i];
                    counter_sList++;
                }
            }
        }
#ifndef BENCHMARKING
        m_logger->debug("control sendMGListsToGPU surface {} | {}",
                counter_sList, size_sList);
#endif
#pragma acc enter data copyin(m_data_MG_sList_level_joined[:size_sList])
    }
}

// ================================= Send obstacle lists to GPU ==========================================
// ***************************************************************************************
/// \brief  create obstacle joined list and send them to GPU
// ***************************************************************************************
void Multigrid::sendObstacleListsToGPU() {
    if (m_numberOfObstacles > 0) {
        size_t size_oFront = getLen_oFront_joined();
        size_t size_oBack = getLen_oBack_joined();
        size_t size_oBottom = getLen_oBottom_joined();
        size_t size_oTop = getLen_oTop_joined();
        size_t size_oLeft = getLen_oLeft_joined();
        size_t size_oRight = getLen_oRight_joined();

        m_data_MG_oFront_level_joined = new size_t[size_oFront];
        m_data_MG_oBack_level_joined = new size_t[size_oBack];
        m_data_MG_oBottom_level_joined = new size_t[size_oBottom];
        m_data_MG_oTop_level_joined = new size_t[size_oTop];
        m_data_MG_oLeft_level_joined = new size_t[size_oLeft];
        m_data_MG_oRight_level_joined = new size_t[size_oRight];

        size_t counter_oFront = 0;
        size_t counter_oBack = 0;
        size_t counter_oBottom = 0;
        size_t counter_oTop = 0;
        size_t counter_oLeft = 0;
        size_t counter_oRight = 0;

        for (size_t level = 0; level < m_levels + 1; level++) {
            Obstacle **obstacleList = *(m_MG_obstacleList + level);
            for (size_t o = 0; o < m_numberOfObstacles; o++) {
                Obstacle *obstacle = *(obstacleList + o);
                for (size_t i = 0; i < obstacle->getSize_obstacleFront(); i++) {
                    m_data_MG_oFront_level_joined[counter_oFront] = obstacle->getObstacleFront()[i];
                    counter_oFront++;
                }
                for (size_t i = 0; i < obstacle->getSize_obstacleBack(); i++) {
                    *(m_data_MG_oBack_level_joined + counter_oBack) = obstacle->getObstacleBack()[i];
                    counter_oBack++;
                }
                for (size_t i = 0; i < obstacle->getSize_obstacleBottom(); i++) {
                    *(m_data_MG_oBottom_level_joined + counter_oBottom) = obstacle->getObstacleBottom()[i];
                    counter_oBottom++;
                }
                for (size_t i = 0; i < obstacle->getSize_obstacleTop(); i++) {
                    *(m_data_MG_oTop_level_joined + counter_oTop) = obstacle->getObstacleTop()[i];
                    counter_oTop++;
                }
                for (size_t i = 0; i < obstacle->getSize_obstacleLeft(); i++) {
                    *(m_data_MG_oLeft_level_joined + counter_oLeft) = obstacle->getObstacleLeft()[i];
                    counter_oLeft++;
                }
                for (size_t i = 0; i < obstacle->getSize_obstacleRight(); i++) {
                    *(m_data_MG_oRight_level_joined + counter_oRight) = obstacle->getObstacleRight()[i];
                    counter_oRight++;
                }
            }
        }

        m_data_MG_oList_zero_joined = m_MG_oList[0];// TODO wrong because only one obstacle is used ?
        size_t size_oList = get_size_obstacle_list();
#pragma acc enter data copyin(m_data_MG_oFront_level_joined[:size_oFront])
#pragma acc enter data copyin(m_data_MG_oBack_level_joined[:size_oBack])
#pragma acc enter data copyin(m_data_MG_oTop_level_joined[:size_oTop])
#pragma acc enter data copyin(m_data_MG_oBottom_level_joined[:size_oBottom])
#pragma acc enter data copyin(m_data_MG_oLeft_level_joined[:size_oLeft])
#pragma acc enter data copyin(m_data_MG_oRight_level_joined[:size_oRight])
#pragma acc enter data copyin(m_data_MG_oList_zero_joined[:size_oList])
    }
}

// ================================= Apply boundary condition ==========================================
// ***************************************************************************************
/// \brief  Apply boundary condition for obstacles, surfaces and domain
/// \param  d Field
/// \param  level Multigrid level
/// \param  f Field type
/// \param  sync synchronous kernel launching (true, default: false)
// ***************************************************************************************
void Multigrid::applyBoundaryCondition(real *d, size_t level, FieldType f, bool sync) {
    if (m_numberOfSurfaces > 0) {
        Surface **surfaceList = *(m_MG_surfaceList + level);
        for (size_t id = 0; id < m_numberOfSurfaces; ++id) {
            ((Surface *) *(surfaceList + id))->applyBoundaryConditions(d, f, level, sync);
        }
    }
    if (m_numberOfObstacles > 0) {
        for (size_t id = 0; id < m_numberOfObstacles - 1; ++id) {
            size_t opatch_start[] = {getFirstIndex_oFront(level, id), getFirstIndex_oBack(level, id), getFirstIndex_oBottom(level, id), getFirstIndex_oTop(level, id), getFirstIndex_oLeft(level, id), getFirstIndex_oRight(level, id)};
            size_t opatch_end[] = {getFirstIndex_oFront(level, id+1), getFirstIndex_oBack(level, id+1), getFirstIndex_oBottom(level, id+1), getFirstIndex_oTop(level, id+1), getFirstIndex_oLeft(level, id+1), getFirstIndex_oRight(level, id+1)};
            (static_cast<BoundaryDataController *> (*(m_bdc_obstacle + id)))->apply_boundary_condition_obstacle(d, m_data_obstacles_patches_joined, opatch_start, opatch_end, f, level, id, sync);
        }
        size_t opatch_start[] = {getFirstIndex_oFront(level, m_numberOfObstacles-1), getFirstIndex_oBack(level, m_numberOfObstacles-1), getFirstIndex_oBottom(level, m_numberOfObstacles-1), getFirstIndex_oTop(level, m_numberOfObstacles-1), getFirstIndex_oLeft(level, m_numberOfObstacles-1), getFirstIndex_oRight(level, m_numberOfObstacles-1)};
        size_t opatch_end[] = {getFirstIndex_oFront(level+1, 0), getFirstIndex_oBack(level+1, 0), getFirstIndex_oBottom(level+1, 0), getFirstIndex_oTop(level+1, 0), getFirstIndex_oLeft(level+1, 0), getFirstIndex_oRight(level+1, 0)};
        ((BoundaryDataController *) *(m_bdc_obstacle + m_numberOfObstacles - 1))->apply_boundary_condition_obstacle(d, m_data_obstacles_patches_joined, opatch_start, opatch_end, f, level, m_numberOfObstacles - 1, sync);
    }

    size_t patch_start[] = {getFirstIndex_bSliceZ(level), getFirstIndex_bSliceZ(level), getFirstIndex_bSliceY(level), getFirstIndex_bSliceY(level), getFirstIndex_bSliceX(level), getFirstIndex_bSliceX(level)};
    size_t patch_end[] = {getFirstIndex_bSliceZ(level + 1), getFirstIndex_bSliceZ(level + 1), getFirstIndex_bSliceY(level + 1), getFirstIndex_bSliceY(level + 1), getFirstIndex_bSliceX(level + 1), getFirstIndex_bSliceX(level + 1)};
    m_bdc_boundary->apply_boundary_condition(d, m_data_boundary_patches_joined, patch_start, patch_end, f, level, sync);
}

//======================================== Update lists ====================================
// ***************************************************************************************
/// \brief  Updates lists of indices
// ***************************************************************************************
void Multigrid::updateLists() {
    removeBoundaryListsFromGPU();

    *(m_size_MG_iList_level) = 0;
    *(m_size_MG_bList_level) = 0;

    *(m_size_MG_bSliceZ_level) = 0;
    *(m_size_MG_bSliceY_level) = 0;
    *(m_size_MG_bSliceX_level) = 0;

    if (m_numberOfObstacles > 0) {
        for (size_t level = 0; level < m_levels + 1; level++) {
            Boundary *boundary = *(m_MG_boundaryList + level);
            boundary->update_lists(*(m_MG_obstacleList + level), m_numberOfObstacles, getSize_oList(level));
            m_size_MG_iList_level[level + 1] = m_size_MG_iList_level[level] + boundary->get_size_inner_list();
            m_size_MG_bList_level[level + 1] = m_size_MG_bList_level[level] + boundary->get_size_boundary_list();
            m_size_MG_bSliceZ_level[level + 1] = m_size_MG_bSliceZ_level[level] + boundary->get_size_boundary_front();
            m_size_MG_bSliceY_level[level + 1] = m_size_MG_bSliceY_level[level] + boundary->get_size_boundary_top();
            m_size_MG_bSliceX_level[level + 1] = m_size_MG_bSliceX_level[level] + boundary->get_size_boundary_left();
        }
    } else {
        for (size_t level = 0; level < m_levels + 1; level++) {
            Boundary *boundary = *(m_MG_boundaryList + level);
            boundary->update_lists();
            m_size_MG_iList_level[level + 1] = m_size_MG_iList_level[level] + boundary->get_size_inner_list();
            m_size_MG_bList_level[level + 1] = m_size_MG_bList_level[level] + boundary->get_size_boundary_list();
            m_size_MG_bSliceZ_level[level + 1] = m_size_MG_bSliceZ_level[level] + boundary->get_size_boundary_front();
            m_size_MG_bSliceY_level[level + 1] = m_size_MG_bSliceY_level[level] + boundary->get_size_boundary_top();
            m_size_MG_bSliceX_level[level + 1] = m_size_MG_bSliceX_level[level] + boundary->get_size_boundary_left();
        }
    }
    sendBoundaryListsToGPU();
    m_data_boundary_patches_joined[Patch::FRONT] = m_data_MG_bFront_level_joined;
    m_data_boundary_patches_joined[Patch::BACK] = m_data_MG_bBack_level_joined;
    m_data_boundary_patches_joined[Patch::BOTTOM] = m_data_MG_bBottom_level_joined;
    m_data_boundary_patches_joined[Patch::TOP] = m_data_MG_bTop_level_joined;
    m_data_boundary_patches_joined[Patch::LEFT] = m_data_MG_bLeft_level_joined;
    m_data_boundary_patches_joined[Patch::RIGHT] = m_data_MG_bRight_level_joined;
}

void Multigrid::removeBoundaryListsFromGPU() {
    size_t size_iList = getLen_iList_joined();
    size_t size_bList = getLen_bList_joined();
#pragma acc exit data delete(m_data_MG_iList_level_joined[:size_iList])
#pragma acc exit data delete(m_data_MG_bList_level_joined[:size_bList])
    delete[] m_data_MG_inner_list_level_joined;
    delete[] m_data_MG_bList_level_joined;

    size_t size_bSliceZ = getLen_bSliceZ_joined();
    size_t size_bSliceY = getLen_bSliceY_joined();
    size_t size_bSliceX = getLen_bSliceX_joined();
#pragma acc exit data delete(m_data_MG_bFront_level_joined[:size_bSliceZ])
#pragma acc exit data delete(m_data_MG_bBack_level_joined[:size_bSliceZ])
#pragma acc exit data delete(m_data_MG_bTop_level_joined[:size_bSliceY])
#pragma acc exit data delete(m_data_MG_bBottom_level_joined[:size_bSliceY])
#pragma acc exit data delete(m_data_MG_bLeft_level_joined[:size_bSliceX])
#pragma acc exit data delete(m_data_MG_bRight_level_joined[:size_bSliceX])
    delete[] m_data_MG_bFront_level_joined;
    delete[] m_data_MG_bBack_level_joined;
    delete[] m_data_MG_bTop_level_joined;
    delete[] m_data_MG_bBottom_level_joined;
    delete[] m_data_MG_bLeft_level_joined;
    delete[] m_data_MG_bRight_level_joined;
}

//======================================== Private getter ====================================

// ***************************************************************************************
/// \brief  get size of obstacle list for level
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getSize_oList(size_t level) {
    size_t size_oList = 0;
    if (m_numberOfObstacles > 0) {
        size_oList = *(m_size_MG_oList_level + level);
    }
    return size_oList;
}

// ***************************************************************************************
/// \brief  get length of the joined inner cell List
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_iList_joined() {
    return getFirstIndex_iList(m_levels + 1) + 1;
}

// ***************************************************************************************
/// \brief  get the index of joined inner cell list, where the first cell of of level l is
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getFirstIndex_iList(size_t level) {
    return *(m_size_MG_iList_level + level);
}

// ***************************************************************************************
/// \brief  get the index of joined inner cell list, where the last cell of level l is
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLastIndex_iList(size_t level) {
    return *(m_size_MG_iList_level + level + 1) - 1;
}

// bList

// ***************************************************************************************
/// \brief  get the index of joined inner cell list, where the first cell of of level l is
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_bList_joined() {
    return getFirstIndex_bList(m_levels + 1) + 1;
}

// ***************************************************************************************
/// \brief  get the index of joined boundary cell list, where the first cell of of level l is
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getFirstIndex_bList(size_t level) {
    return *(m_size_MG_bList_level + level);
}

// ***************************************************************************************
/// \brief  get the index of joined boundary cell list, where the last cell of of level l is
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLastIndex_bList(size_t level) {
    return *(m_size_MG_bList_level + level + 1) - 1;
}

// sList
size_t Multigrid::getLen_sList_joined() {
    size_t len = 0;
    if (m_numberOfSurfaces > 0) {
        len = getFirstIndex_sList(m_levels + 1) + 1;
    }
    return len;
}

size_t Multigrid::getFirstIndex_sList(size_t level) {
    size_t index = 0;
    if (m_numberOfSurfaces > 0) {
        index = *(m_size_MG_sList_level + level * m_numberOfSurfaces);
    }
    return index;
}

size_t Multigrid::getLastIndex_sList(size_t level) {
    size_t index = 0;
    if (m_numberOfSurfaces > 0) {
        index = *(m_size_MG_sList_level + (level + 1) * m_numberOfSurfaces) - 1;
    }
    return index;
}

// bSlices

// ***************************************************************************************
/// \brief  get length of left/right joined boundary list
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_bSliceX_joined() {
    return getFirstIndex_bSliceX(m_levels + 1) + 1;
}

// ***************************************************************************************
/// \brief  get length of top/bottom joined boundary list
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_bSliceY_joined() {
    return getFirstIndex_bSliceY(m_levels + 1) + 1;
}

// ***************************************************************************************
/// \brief  get length of front/back joined boundary list
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_bSliceZ_joined() {
    return getFirstIndex_bSliceZ(m_levels + 1) + 1;
}

// ***************************************************************************************
/// \brief  get the index of joined boundary cell list of left/right patch of the first cell of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getFirstIndex_bSliceX(size_t level) {
    return *(m_size_MG_bSliceX_level + level);
}

// ***************************************************************************************
/// \brief  get the index of joined boundary cell list of top/bottom patch of the first cell of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getFirstIndex_bSliceY(size_t level) {
    return *(m_size_MG_bSliceY_level + level);
}

// ***************************************************************************************
/// \brief  get the index of joined boundary cell list of front/back patch of the first cell of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getFirstIndex_bSliceZ(size_t level) {
    return *(m_size_MG_bSliceZ_level + level);
}

// ***************************************************************************************
/// \brief  get the index of joined boundary cell list of left/right patch of the last cell of of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLastIndex_bSliceX(size_t level) {
    return *(m_size_MG_bSliceX_level + level + 1) - 1;
}

// ***************************************************************************************
/// \brief  get the index of joined boundary cell list of top/bottom patch of the last cell of of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLastIndex_bSliceY(size_t level) {
    return *(m_size_MG_bSliceY_level + level + 1) - 1;
}

// ***************************************************************************************
/// \brief  get the index of joined boundary cell list of front/back patch of the last cell of of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLastIndex_bSliceZ(size_t level) {
    return *(m_size_MG_bSliceZ_level + level + 1) - 1;
}

// oSlices

// ***************************************************************************************
/// \brief  get length of m_data_MG_oFront_level_joined
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oFront_joined() {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        size_t index = (m_levels + 1) * m_numberOfObstacles;
        if (m_size_MG_oFront_level[index] > 0) {
            len = m_size_MG_oFront_level[index] + 1;
        }
    }
    return len;
}

// ***************************************************************************************
/// \brief  get length of m_data_MG_oBack_level_joined
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oBack_joined() {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        size_t index = (m_levels + 1) * m_numberOfObstacles;
        if (m_size_MG_oBack_level[index] > 0) {
            len = m_size_MG_oBack_level[index] + 1;
        }
    }
    return len;
}

// ***************************************************************************************
/// \brief  get length of m_data_MG_oBottom_level_joined
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oBottom_joined() {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        size_t index = (m_levels + 1) * m_numberOfObstacles;
        if (m_size_MG_oBottom_level[index] > 0) {
            len = m_size_MG_oBottom_level[index] + 1;
        }
    }
    return len;
}

// ***************************************************************************************
/// \brief  get length of m_data_MG_oTop_level_joined
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oTop_joined() {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        size_t index = (m_levels + 1) * m_numberOfObstacles;
        if (m_size_MG_oTop_level[index] > 0) {
            len = m_size_MG_oTop_level[index] + 1;
        }
    }
    return len;
}

// ***************************************************************************************
/// \brief  get length of m_data_MG_oLeft_level_joined
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oLeft_joined() {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        size_t index = (m_levels + 1) * m_numberOfObstacles;
        if (m_size_MG_oLeft_level[index] > 0) {
            len = m_size_MG_oLeft_level[index] + 1;
        }
    }
    return len;
}

// ***************************************************************************************
/// \brief  get length of m_data_MG_oRight_level_joined
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oRight_joined() {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        size_t index = (m_levels + 1) * m_numberOfObstacles;
        if (m_size_MG_oRight_level[index] > 0) {
            len = m_size_MG_oRight_level[index] + 1;
        }
    }
    return len;
}

// ***************************************************************************************
/// \brief  get amount of Front obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oFront(size_t level) {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        size_t last = getLastIndex_oFront(level, m_numberOfObstacles - 1);
        if (last > 0) {
            len = last - getFirstIndex_oFront(level, 0) + 1;
        }
    }
    return len;
}

// ***************************************************************************************
/// \brief get amount of Back obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oBack(size_t level) {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        size_t last = getLastIndex_oBack(level, m_numberOfObstacles - 1);
        if (last > 0) {
            len = last - getFirstIndex_oBack(level, 0) + 1;
        }
    }
    return len;
}

// ***************************************************************************************
/// \brief  get amount of Bottom obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oBottom(size_t level) {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        size_t last = getLastIndex_oBottom(level, m_numberOfObstacles - 1);
        if (last > 0) {
            len = last - getFirstIndex_oBottom(level, 0) + 1;
        }
    }
    return len;
}

// ***************************************************************************************
/// \brief  get amount of Top obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oTop(size_t level) {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        size_t last = getLastIndex_oTop(level, m_numberOfObstacles - 1);
        if (last > 0) {
            len = last - getFirstIndex_oTop(level, 0) + 1;
        }
    }
    return len;
}

// ***************************************************************************************
/// \brief  get amount of Left obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oLeft(size_t level) {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        size_t last = getLastIndex_oLeft(level, m_numberOfObstacles - 1);
        if (last > 0) {
            len = last - getFirstIndex_oLeft(level, 0) + 1;
        }
    }
    return len;
}

// ***************************************************************************************
/// \brief  get amount of Front obstacle cells of level l
/// \param level Multigrid level
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getLen_oRight(size_t level) {
    size_t len = 0;
    if (m_numberOfObstacles > 0) {
        len = getFirstIndex_oRight(level + 1, 0) - getFirstIndex_oRight(level, 0);
    }
    return len;
}

// ***************************************************************************************
/// \brief  get the index of joined obstacle cell list of Front patch of the first cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getFirstIndex_oFront(size_t level, size_t id) {
    size_t index = 0;
    if (m_numberOfObstacles > 0) {
        index = *(m_size_MG_oFront_level + level * m_numberOfObstacles + id);
    }
    return index;
}

// ***************************************************************************************
/// \brief  get the index of joined obstacle cell list of Back patch of the first cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getFirstIndex_oBack(size_t level, size_t id) {
    size_t index = 0;
    if (m_numberOfObstacles > 0) {
        index = *(m_size_MG_oBack_level + level * m_numberOfObstacles + id);
    }
    return index;
}

// ***************************************************************************************
/// \brief  get the index of joined obstacle cell list of Bottom patch of the first cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getFirstIndex_oBottom(size_t level, size_t id) {
    size_t index = 0;
    if (m_numberOfObstacles > 0) {
        index = *(m_size_MG_oBottom_level + level * m_numberOfObstacles + id);
    }
    return index;
}

// ***************************************************************************************
/// \brief  get the index of joined obstacle cell list of Top patch of the first cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getFirstIndex_oTop(size_t level, size_t id) {
    size_t index = 0;
    if (m_numberOfObstacles > 0) {
        index = *(m_size_MG_oTop_level + level * m_numberOfObstacles + id);
    }
    return index;
}

// ***************************************************************************************
/// \brief  get the index of joined obstacle cell list of Left patch of the first cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getFirstIndex_oLeft(size_t level, size_t id) {
    size_t index = 0;
    if (m_numberOfObstacles > 0) {
        index = *(m_size_MG_oLeft_level + level * m_numberOfObstacles + id);
    }
    return index;
}

// ***************************************************************************************
/// \brief get the index of joined obstacle cell list of right patch of the first cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// ***************************************************************************************
size_t Multigrid::getFirstIndex_oRight(size_t level, size_t id) {
    size_t index = 0;
    if (m_numberOfObstacles > 0) {
        index = *(m_size_MG_oRight_level + level * m_numberOfObstacles + id);
    }
    return index;
}

// ***************************************************************************************
/// \brief get the index of joined obstacle cell list of front patch of the last cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// ***************************************************************************************
int Multigrid::getLastIndex_oFront(size_t level, size_t id) {
    int index = 0;
    if (m_numberOfObstacles > 0) {
        index = m_size_MG_oFront_level[level * m_numberOfObstacles + id + 1];
    }
    return index;
}

// ***************************************************************************************
/// \brief get the index of joined obstacle cell list of back patch of the last cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// ***************************************************************************************
int Multigrid::getLastIndex_oBack(size_t level, size_t id) {
    int index = 0;
    if (m_numberOfObstacles > 0) {
        index = m_size_MG_oBack_level[level * m_numberOfObstacles + id + 1];
    }
    return index;
}

// ***************************************************************************************
/// \brief get the index of joined obstacle cell list of bottom patch of the last cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return size_t
// ***************************************************************************************
int Multigrid::getLastIndex_oBottom(size_t level, size_t id) {
    size_t index = 0;
    if (m_numberOfObstacles > 0) {
        index = m_size_MG_oBottom_level[level * m_numberOfObstacles + id + 1] - 1;
    }
    return index;
}

// ***************************************************************************************
/// \brief get the index of joined obstacle cell list of top patch of the last cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return int
// ***************************************************************************************
int Multigrid::getLastIndex_oTop(size_t level, size_t id) {
    int index = 0;
    if (m_numberOfObstacles > 0) {
        index = m_size_MG_oTop_level[level * m_numberOfObstacles + id + 1] - 1;
    }
    return index;
}

// ***************************************************************************************
/// \brief get the index of joined obstacle cell list of left patch of the last cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return int
// ***************************************************************************************
int Multigrid::getLastIndex_oLeft(size_t level, size_t id) {
    int index = 0;
    if (m_numberOfObstacles > 0) {
        index = m_size_MG_oLeft_level[level * m_numberOfObstacles + id + 1];
    }
    return index;
}

// ***************************************************************************************
/// \brief get the index of joined obstacle cell list of right patch of the last cell of of level l
/// \param level Multigrid level
/// \param id ID of obstacle
/// \return int
// ***************************************************************************************
int Multigrid::getLastIndex_oRight(size_t level, size_t id) {
    int index = 0;
    if (m_numberOfObstacles > 0) {
        index = m_size_MG_oRight_level[level * m_numberOfObstacles + id + 1];
    }
    return index;
}

// ---------------- public getter
size_t Multigrid::get_size_inner_list(size_t level) {
    return getLastIndex_iList(level) - getFirstIndex_iList(level) + 1;
}

size_t Multigrid::get_size_boundary_ist(size_t level) {
    return getLastIndex_bList(level) - getFirstIndex_bList(level) + 1;
}

size_t Multigrid::get_size_obstacle_list() {
    return getSize_oList(0);
}

size_t* Multigrid::get_obstacle_list() {
    if (m_numberOfObstacles > 0) {
        return m_data_MG_oList_zero_joined;
    } else {
        return nullptr;
    }
}

size_t Multigrid::getInnerList_level_joined_start(size_t level) {
    return *(m_size_MG_iList_level + level);
}

size_t Multigrid::getInnerList_level_joined_end(size_t level) {
    return getInnerList_level_joined_start(level + 1) - 1;
}

size_t Multigrid::getBoundaryList_level_joined_start(size_t level) {
    return *(m_size_MG_bList_level + level);
}

size_t Multigrid::getBoundaryList_level_joined_end(size_t level) {
    return getBoundaryList_level_joined_start(level + 1) - 1;
}

size_t Multigrid::getObstacleStrideX(size_t id, size_t level) {
    return (static_cast<Obstacle *>(m_MG_obstacleList[level][id]))->getStrideX();
}

size_t Multigrid::getObstacleStrideY(size_t id, size_t level) {
    return (static_cast<Obstacle *>(m_MG_obstacleList[level][id]))->getStrideY();
}

size_t Multigrid::getObstacleStrideZ(size_t id, size_t level) {
    return (static_cast<Obstacle *>(m_MG_obstacleList[level][id]))->getStrideZ();
}

