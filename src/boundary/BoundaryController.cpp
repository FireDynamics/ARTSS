/// \file       BoundaryController.cpp
/// \brief      Controll class for boundary
/// \date       Oct 01, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "BoundaryController.h"
#include <string>

BoundaryController *BoundaryController::singleton = nullptr;  // Singleton


BoundaryController::BoundaryController() {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_bdc_boundary = new BoundaryDataController();
    read_XML();
    if (m_number_of_obstacles + m_number_of_surfaces > 0) {
        m_multigrid = new Multigrid(m_number_of_surfaces, m_surface_list,
                                    m_number_of_obstacles, m_obstacle_list,
                                    m_bdc_boundary, m_bdc_obstacles);
    } else {
        m_multigrid = new Multigrid(m_bdc_boundary);
    }
#ifndef BENCHMARKING
    print_boundaries();
#endif
}

// ================================= Read XML ======================================================
// *************************************************************************************************
/// \brief  Reads in all parameters of boundary, obstacles and surfaces
// *************************************************************************************************
void BoundaryController::read_XML() {
#ifndef BENCHMARKING
    m_logger->debug("start parsing XML");
#endif

    auto params = Parameters::getInstance();
    parse_boundary_parameter(params->get_first_child("boundaries"));
    parse_obstacle_parameter(params->get_first_child("obstacles"));
    detect_neighbouring_obstacles();
    parse_surface_parameter(params->get_first_child("surfaces"));
#ifndef BENCHMARKING
    m_logger->debug("finished parsing XML");
#endif
}

// ================================= Parser ========================================================
// *************************************************************************************************
/// \brief  parses boundaries of domain from XML file
/// \param  xmlParameter pointer to XMLElement to start with
// *************************************************************************************************
void BoundaryController::parse_boundary_parameter(tinyxml2::XMLElement *xmlParameter) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing boundary parameter");
#endif
// BOUNDARY
    auto curElem = xmlParameter->FirstChildElement();
    while (curElem) {
        m_bdc_boundary->add_boundary_data(curElem);
        curElem = curElem->NextSiblingElement();
    }
#ifndef BENCHMARKING
    m_logger->debug("finished parsing boundary parameter");
#endif
}

// ================================= Parser ========================================================
// *************************************************************************************************
/// \brief  parses surfaces from XML file
/// \param  xmlParameter pointer to XMLElement to start with
// *************************************************************************************************
void BoundaryController::parse_surface_parameter(tinyxml2::XMLElement *xmlParameter) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing surface parameter");
#endif
    // SURFACES
    // TODO(issue 5): surfaces
    m_has_surfaces = (Parameters::getInstance()->get("surfaces/enabled") == "Yes");
    if (m_has_surfaces) {
        std::vector<Surface *> surfaces;
        auto cur_elem = xmlParameter->FirstChildElement();
        while (cur_elem) {
            Surface *o = new Surface(cur_elem);
            surfaces.push_back(o);
            cur_elem = cur_elem->NextSiblingElement();
        }
        m_number_of_surfaces = surfaces.size();
        m_surface_list = surfaces.data();
    }
#ifndef BENCHMARKING
    m_logger->debug("finished parsing surface parameter");
#endif
}

// ================================= Parser ========================================================
// *************************************************************************************************
/// \brief  parses obstacles from XML file
/// \param  xmlParameter pointer to XMLElement to start with
// *************************************************************************************************
void BoundaryController::parse_obstacle_parameter(tinyxml2::XMLElement *xmlParameter) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing obstacle parameter");
#endif
// OBSTACLES
    m_has_obstacles = (Parameters::getInstance()->get("obstacles/enabled") == "Yes");
    if (m_has_obstacles) {
        std::vector<Obstacle *> obstacles;
        std::vector<BoundaryDataController *> bdc_obstacles;
        auto cur_elem_obstacle = xmlParameter->FirstChildElement();
        while (cur_elem_obstacle) {
            std::string name = cur_elem_obstacle->Attribute("name");
#ifndef BENCHMARKING
            m_logger->debug("read obstacle '{}'", name);
#endif
            BoundaryDataController *bdc = new BoundaryDataController();
            auto cur_elem = cur_elem_obstacle->FirstChildElement();
            real ox1;
            real ox2;
            real oy1;
            real oy2;
            real oz1;
            real oz2;
            while (cur_elem) {
                std::string nodeName = cur_elem->Value();
                if (nodeName == "boundary") {
                    bdc->add_boundary_data(cur_elem);
                } else if (nodeName == "geometry") {
                    ox1 = cur_elem->DoubleAttribute("ox1");
                    ox2 = cur_elem->DoubleAttribute("ox2");
                    oy1 = cur_elem->DoubleAttribute("oy1");
                    oy2 = cur_elem->DoubleAttribute("oy2");
                    oz1 = cur_elem->DoubleAttribute("oz1");
                    oz2 = cur_elem->DoubleAttribute("oz2");
                } else {
#ifndef BENCHMARKING
                    m_logger->warn("Ignoring unknown node {}", nodeName);
#endif
                }
                cur_elem = cur_elem->NextSiblingElement();
            }
            Obstacle *o = new Obstacle(ox1, ox2, oy1, oy2, oz1, oz2, name);
            obstacles.push_back(o);
            bdc_obstacles.push_back(bdc);
            cur_elem_obstacle = cur_elem_obstacle->NextSiblingElement();
        }
        m_number_of_obstacles = obstacles.size();
        m_obstacle_list = new Obstacle *[m_number_of_obstacles];
        m_bdc_obstacles = new BoundaryDataController *[m_number_of_obstacles];
        for (size_t i = 0; i < m_number_of_obstacles; i++) {
            *(m_obstacle_list + i) = obstacles[i];
            *(m_bdc_obstacles + i) = bdc_obstacles[i];
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("finished parsing obstacle parameter");
#endif
}

BoundaryController::~BoundaryController() {
    delete (m_multigrid);
    delete (m_bdc_boundary);
    for (size_t i = 0; i < m_number_of_obstacles; i++) {
        delete (m_bdc_obstacles[i]);
    }
}


BoundaryController *BoundaryController::getInstance() {
    if (singleton == nullptr) {
        singleton = new BoundaryController();
    }
    return singleton;
}

// ================================= Printer =======================================================
// *************************************************************************************************
/// \brief  prints boundaries (outer, inner, surfaces)
// *************************************************************************************************
void BoundaryController::print_boundaries() {
#ifndef BENCHMARKING
    m_logger->info("-- Info summary");
    Domain::getInstance()->print();
    m_bdc_boundary->print();
    for (size_t i = 0; i < m_number_of_obstacles; i++) {
        m_obstacle_list[i]->print();
        m_bdc_obstacles[i]->print();
    }
    for (size_t i = 0; i < m_number_of_surfaces; i++) {
        m_surface_list[i]->print();
    }
#endif
}

//======================================== Update lists ============================================
// *************************************************************************************************
/// \brief  Updates lists of indices
// *************************************************************************************************
void BoundaryController::update_lists() {
    m_multigrid->update_lists();
}

//======================================== Update lists ============================================
// *************************************************************************************************
/// \brief  Apply boundary for level 0
/// \param  d data field
/// \param  f Field
/// \param  sync synchronous kernel launching (true, default: false)
// *************************************************************************************************
void BoundaryController::apply_boundary(real *d, FieldType f, bool sync) {
    apply_boundary(d, 0, f, sync);
}

// ================================= Apply BCs in level l > 0 ======================================
// *************************************************************************************************
/// \brief  applies zero-value boundary conditions to all boundaries (domain, surfaces, obstacles)
/// for level l
/// \param  d           data field
/// \param  level       Multigrid level
/// \param  f           type of output pointer
/// \param  sync    synchronization (default: false)
// *************************************************************************************************
void BoundaryController::apply_boundary(real *d, size_t level, FieldType f, bool sync) {
    m_multigrid->apply_boundary_condition(d, level, f, sync);
}

// ================================= detect neighbouring obstacles =================================
// *************************************************************************************************
/// \brief  handling of neighbouring obstacles, removes cells with circular constraints which never
/// change their value
// *************************************************************************************************
void BoundaryController::detect_neighbouring_obstacles() {
#ifndef BENCHMARKING
    m_logger->debug("start detecting neighbouring obstacles");
#endif
    for (size_t o1 = 0; o1 < m_number_of_obstacles; o1++) {
        Obstacle *obstacle1 = m_obstacle_list[o1];
        for (size_t o2 = o1 + 1; o2 < m_number_of_obstacles; o2++) {
            Obstacle *obstacle2 = m_obstacle_list[o2];
            Obstacle::remove_circular_constraints(obstacle1, obstacle2);
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("finished detecting neighbouring obstacles");
#endif
}

size_t BoundaryController::get_size_inner_list_level_joined() const {
    return m_multigrid->get_size_inner_list_level_joined();
}

size_t BoundaryController::get_size_boundary_list_level_joined() const {
    return m_multigrid->get_size_boundary_list_level_joined();
}

size_t *BoundaryController::get_obstacle_list() const {
    return m_multigrid->get_obstacle_list();
}

size_t BoundaryController::get_size_boundary_list() const {
    return m_multigrid->get_size_boundary_list();
}

size_t BoundaryController::get_size_inner_list() const {
    return m_multigrid->get_size_inner_list();
}

size_t BoundaryController::get_size_obstacle_list() const {
    return m_multigrid->get_size_obstacle_list();
}

size_t *BoundaryController::get_inner_list_level_joined() const {
    return m_multigrid->get_inner_list_level_joined();
}

size_t BoundaryController::get_inner_list_level_joined_start(size_t level) const {
    return m_multigrid->get_inner_list_level_joined_start(level);
}

size_t BoundaryController::get_inner_list_level_joined_end(size_t level) const {
    return m_multigrid->get_inner_list_level_joined_end(level);
}

size_t *BoundaryController::get_boundary_list_level_joined() const {
    return m_multigrid->get_boundary_list_level_joined();
}

size_t BoundaryController::get_boundary_list_level_joined_start(size_t level) const {
    return m_multigrid->get_boundary_list_level_joined_start(level);
}

size_t BoundaryController::get_boundary_list_level_joined_end(size_t level) const {
    return m_multigrid->get_boundary_list_level_joined_end(level);
}

size_t BoundaryController::get_obstacle_stride_x(size_t id, size_t level) const {
    return m_multigrid->get_obstacle_stride_x(id, level);
}

size_t BoundaryController::get_obstacle_stride_y(size_t id, size_t level) const {
    return m_multigrid->get_obstacle_stride_y(id, level);
}

size_t BoundaryController::get_obstacle_stride_z(size_t id, size_t level) const {
    return m_multigrid->get_obstacle_stride_z(id, level);
}

std::vector<FieldType> BoundaryController::get_used_fields() const {
    return m_bdc_boundary->get_used_fields();
}
