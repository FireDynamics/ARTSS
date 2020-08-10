/// \file 		BoundaryController.cpp
/// \brief 		Controll class for boundary
/// \date 		Oct 01, 2020
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "BoundaryController.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include <tuple>
#include "../utility/Utility.h"
#include <algorithm>

BoundaryController *BoundaryController::singleton = nullptr; // Singleton


BoundaryController::BoundaryController() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_bdc_boundary = new BoundaryDataController();
    readXML();
    if (m_numberOfObstacles + m_numberOfSurfaces > 0) {
        m_multigrid = new Multigrid(m_numberOfSurfaces, m_surfaceList, m_numberOfObstacles, m_obstacleList, m_bdc_boundary, m_bdc_obstacles);
    } else {
        m_multigrid = new Multigrid(m_bdc_boundary);
    }
#ifndef BENCHMARKING
    printBoundaries();
#endif
}

// ================================= Read XML =============================================
// ***************************************************************************************
/// \brief  Reads in all parameters of boundary, obstacles and surfaces
// ***************************************************************************************
void BoundaryController::readXML() {
    auto params = Parameters::getInstance();
    parseBoundaryParameter(params->get_first_child("boundaries"));
    parseObstacleParameter(params->get_first_child("obstacles"));
    parseSurfaceParameter( params->get_first_child("surfaces"));
}

// ================================= Parser =============================================
// ***************************************************************************************
/// \brief  parses boundaries of domain from XML file
/// \param	xmlParameter pointer to XMLElement to start with
// ***************************************************************************************
void BoundaryController::parseBoundaryParameter(tinyxml2::XMLElement *xmlParameter) {
// BOUNDARY
    auto curElem = xmlParameter->FirstChildElement();
    while (curElem) {
        m_bdc_boundary->addBoundaryData(curElem);
        curElem = curElem->NextSiblingElement();
    } // end while
}

// ================================= Parser =============================================
// ***************************************************************************************
/// \brief  parses surfaces from XML file
/// \param	xmlParameter pointer to XMLElement to start with
// ***************************************************************************************
void BoundaryController::parseSurfaceParameter(tinyxml2::XMLElement *xmlParameter) {
// SURFACES
//TODO
    m_hasSurfaces = (Parameters::getInstance()->get("surfaces/enabled") == "Yes");
    if (m_hasSurfaces) {
        std::vector<Surface *> surfaces;
        auto curElem = xmlParameter->FirstChildElement();
        while (curElem) {
            Surface *o = new Surface(curElem);
            surfaces.push_back(o);
            curElem = curElem->NextSiblingElement();
        } // end while
        m_numberOfSurfaces = surfaces.size();
        m_surfaceList = surfaces.data();
    }
}

// ================================= Parser =============================================
// ***************************************************************************************
/// \brief  parses obstacles from XML file
/// \param	xmlParameter pointer to XMLElement to start with
// ***************************************************************************************
void BoundaryController::parseObstacleParameter(tinyxml2::XMLElement *xmlParameter) {
// OBSTACLES
    m_hasObstacles = (Parameters::getInstance()->get("obstacles/enabled") == "Yes");
    if (m_hasObstacles) {
        std::vector<Obstacle *> obstacles;
        std::vector<BoundaryDataController *> bdc_obstacles;
        auto curElem_obstacle = xmlParameter->FirstChildElement();
        while (curElem_obstacle) {
            BoundaryDataController *bdc = new BoundaryDataController();
            auto curElem = curElem_obstacle->FirstChildElement();
            real ox1;
            real ox2;
            real oy1;
            real oy2;
            real oz1;
            real oz2;
            while (curElem) {
                std::string nodeName = curElem->Value();
                if (nodeName == "boundary") {
                    bdc->addBoundaryData(curElem);
                } else if (nodeName == "geometry") {
                    ox1 = curElem->DoubleAttribute("ox1");
                    ox2 = curElem->DoubleAttribute("ox2");
                    oy1 = curElem->DoubleAttribute("oy1");
                    oy2 = curElem->DoubleAttribute("oy2");
                    oz1 = curElem->DoubleAttribute("oz1");
                    oz2 = curElem->DoubleAttribute("oz2");
                } else {
#ifndef BENCHMARKING
                    m_logger->warn("Ignoring unknown node {}", nodeName);
#endif
                }
                curElem = curElem->NextSiblingElement();
            }
            Obstacle *o = new Obstacle(ox1, ox2, oy1, oy2, oz1, oz2);
            obstacles.push_back(o);
            bdc_obstacles.push_back(bdc);
            curElem_obstacle = curElem_obstacle->NextSiblingElement();
        } // end while
        m_numberOfObstacles = obstacles.size();
        m_obstacleList = new Obstacle *[m_numberOfObstacles];
        m_bdc_obstacles = new BoundaryDataController *[m_numberOfObstacles];
        for (size_t i = 0; i < m_numberOfObstacles; i++) {
            *(m_obstacleList + i) = obstacles[i];
            *(m_bdc_obstacles + i) = bdc_obstacles[i];
        }
    }
}

BoundaryController::~BoundaryController() {
    delete (m_multigrid);
    delete (m_bdc_boundary);
    for (size_t i = 0; i < m_numberOfObstacles; i++) {
        delete (m_bdc_obstacles[i]);
    }
    //for (size_t surface = 0; surface < m_numberOfSurfaces; surface++) {
    //    delete (*(m_surfaceList + surface));
    //}
    //delete[] m_surfaceList;
    //for (size_t obstacle = 0; obstacle < m_numberOfObstacles; obstacle++) {
    //    delete (*(m_obstacleList + obstacle));
    //}
    //delete[] m_obstacleList;
}


BoundaryController *BoundaryController::getInstance() {
    if (singleton == nullptr) {
        singleton = new BoundaryController();
    }
    return singleton;
}

// ================================= Printer =============================================
// ***************************************************************************************
/// \brief  prints boundaries (outer, inner, surfaces)
// ***************************************************************************************
void BoundaryController::printBoundaries() {
#ifdef BENCHMARKING
    return;
#else
    m_logger->info("-- Info summary");
    Domain::getInstance()->print();
    m_bdc_boundary->print();
    for (size_t i = 0; i < m_numberOfObstacles; i++) {
        m_obstacleList[i]->print();
        m_bdc_obstacles[i]->print();
    }
    for (size_t i = 0; i < m_numberOfSurfaces; i++) {
        m_surfaceList[i]->print();
    }
#endif
}

//======================================== Update lists ====================================
// ***************************************************************************************
/// \brief  Updates lists of indices
// ***************************************************************************************
void BoundaryController::updateLists() {
    m_multigrid->updateLists();
}

//======================================== Update lists ====================================
// ***************************************************************************************
/// \brief  Apply boundary for level 0
/// \param  d data field
/// \param  f Field
/// \param  sync synchronous kernel launching (true, default: false)
// ***************************************************************************************
void BoundaryController::applyBoundary(real *d, FieldType f, bool sync) {
    applyBoundary(d, 0, f, sync);
}

// ================================= Apply BCs in level l > 0 ===========================================
// ***************************************************************************************
/// \brief  applies zero-value boundary conditions to all boundaries (domain, surfaces, obstacles) for level l
/// \param 	d 			data field
/// \param 	level 		Multigrid level
/// \param 	f 			type of output pointer
/// \param  sync    synchronization (default: false)
// ***************************************************************************************
void BoundaryController::applyBoundary(real *d, size_t level, FieldType f, bool sync) {
    m_multigrid->applyBoundaryCondition(d, level, f, sync);
}

size_t BoundaryController::getSize_innerList_level_joined() {
    return m_multigrid->getSize_innerList_level_joined();
}

size_t BoundaryController::getSize_boundaryList_level_joined() {
    return m_multigrid->getSize_boundaryList_level_joined();
}

size_t* BoundaryController::get_obstacleList() {
    return m_multigrid->get_obstacleList();
}

size_t BoundaryController::getSize_boundaryList() {
    return m_multigrid->getSize_boundaryList();
}

size_t BoundaryController::getSize_innerList() {
    return m_multigrid->getSize_innerList();
}

size_t BoundaryController::getSize_obstacleList() {
    return m_multigrid->getSize_obstacleList();
}

size_t *BoundaryController::get_innerList_level_joined() {
    return m_multigrid->getInnerList_level_joined();
}

size_t BoundaryController::get_innerList_level_joined_start(size_t level) {
    return m_multigrid->getInnerList_level_joined_start(level);
}

size_t BoundaryController::get_innerList_level_joined_end(size_t level) {
    return m_multigrid->getInnerList_level_joined_end(level);
}

size_t *BoundaryController::get_boundaryList_level_joined() {
    return m_multigrid->getBoundaryList_level_joined();
}

size_t BoundaryController::get_boundaryList_level_joined_start(size_t level) {
    return m_multigrid->getBoundaryList_level_joined_start(level);
}

size_t BoundaryController::get_boundaryList_level_joined_end(size_t level) {
    return m_multigrid->getBoundaryList_level_joined_end(level);
}

size_t BoundaryController::getObstacleStrideX(size_t id, size_t level){
    return m_multigrid->getObstacleStrideX(id, level);
}
size_t BoundaryController::getObstacleStrideY(size_t id, size_t level){
    return m_multigrid->getObstacleStrideY(id, level);
}
size_t BoundaryController::getObstacleStrideZ(size_t id, size_t level){
    return m_multigrid->getObstacleStrideZ(id, level);
}

std::vector<FieldType> BoundaryController::get_used_fields() {
    return m_bdc_boundary->get_used_fields();
}
