/// \file       DomainController.cpp
/// \brief      Controller class for boundary
/// \date       Oct 01, 2020
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "DomainController.h"
#include <string>

DomainController *DomainController::singleton = nullptr;  // Singleton


DomainController::DomainController(Settings::Settings const &settings) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_bdc_boundary = new BoundaryDataController();
    read_XML(settings);
    size_t multigrid_level = DomainData::getInstance()->get_levels();
    m_multigrid = new Multigrid(m_number_of_surfaces, m_surface_list,
                                m_number_of_obstacles, m_obstacle_list,
                                m_bdc_boundary, m_bdc_obstacles, m_bdc_surfaces,
                                multigrid_level);
#ifndef BENCHMARKING
    print_boundaries();
#endif
}

// ================================= Read XML ======================================================
// *************************************************************************************************
/// \brief  Reads in all parameters of boundary, obstacles and surfaces
// *************************************************************************************************
void DomainController::read_XML(Settings::Settings const &settings) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing XML");
#endif
    m_has_obstacles = settings.get_bool("obstacles/enabled");
    m_has_surfaces = settings.get_bool("surfaces/enabled");

    parse_boundary_parameter(settings.get_boundaries());
    parse_obstacle_parameter(settings.get_obstacles());
    detect_neighbouring_obstacles();
    parse_surface_parameter(settings.get_surfaces());
#ifndef BENCHMARKING
    m_logger->debug("finished parsing XML");
#endif
}

// ================================= Parser ========================================================
// *************************************************************************************************
/// \brief  parses boundaries of domain from XML file
/// \param  xmlParameter pointer to XMLElement to start with
// *************************************************************************************************
void DomainController::parse_boundary_parameter(const std::vector<Settings::BoundarySetting> &boundaries) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing boundary parameter");
#endif
// BOUNDARY
    for (const auto &boundary : boundaries) {
        m_bdc_boundary->add_boundary_data(boundary);
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
void DomainController::parse_surface_parameter(const std::vector<Settings::SurfaceSetting> &surfaces) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing surface parameter");
#endif
    // SURFACES
    if (m_has_surfaces) {
        std::vector<Surface *> ret_surfaces;
        std::vector<BoundaryDataController *> bdc_surfaces;
        for (const auto &surface : surfaces) {

            std::string name = surface.get_name();
            real x1 = surface.get_sx1();
            real x2 = surface.get_sx2();
            real y1 = surface.get_sy1();
            real y2 = surface.get_sy2();
            real z1 = surface.get_sz1();
            real z2 = surface.get_sz2();
            auto s = new Surface(x1, x2, y1, y2, z1, z2, name);

            auto bdc = new BoundaryDataController();
            for (const auto &boundary : surface.get_boundaries()) {
                bdc->add_boundary_data(boundary);
            }

            //TODO(issue 5) store surface boundary conditions. surfaces only have one patch
            ret_surfaces.push_back(s);
            bdc_surfaces.push_back(bdc);
        }
        m_number_of_surfaces = ret_surfaces.size();
        m_surface_list = new Surface *[m_number_of_surfaces];
        m_bdc_surfaces = new BoundaryDataController *[m_number_of_surfaces];
        //TODO(cvm) vector surfaces will be destroyed after this method therefore I need to preserve the data
        //std::memcpy(m_surface_list, ret_surfaces.data(), ret_surfaces.size());
        for (size_t i = 0; i < m_number_of_surfaces; i++) {
            m_surface_list[i] = ret_surfaces[i];
            m_bdc_surfaces[i] = bdc_surfaces[i];
        }
    } else {
        m_surface_list = new Surface *[m_number_of_surfaces];
        m_bdc_surfaces = new BoundaryDataController *[m_number_of_surfaces];
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
void DomainController::parse_obstacle_parameter(const std::vector<Settings::ObstacleSetting> &obstacles) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing obstacle parameter");
#endif
// OBSTACLES
    if (m_has_obstacles) {
        std::vector<Obstacle *> ret_obstacles;
        std::vector<BoundaryDataController *> bdc_obstacles;

        for (const auto &obstacle : obstacles) {
            std::string name = obstacle.get_name();
#ifndef BENCHMARKING
            m_logger->debug("read obstacle '{}'", name);
#endif
            auto bdc = new BoundaryDataController();
            real ox1 = obstacle.get_ox1();
            real ox2 = obstacle.get_ox2();
            real oy1 = obstacle.get_oy1();
            real oy2 = obstacle.get_oy2();
            real oz1 = obstacle.get_oz1();
            real oz2 = obstacle.get_oz2();

            for (const auto &boundary : obstacle.get_boundaries()) {
                bdc->add_boundary_data(boundary);
            }

            auto o = new Obstacle(ox1, ox2, oy1, oy2, oz1, oz2, name);
            ret_obstacles.push_back(o);
            bdc_obstacles.push_back(bdc);
        }

        m_number_of_obstacles = ret_obstacles.size();
        m_obstacle_list = new Obstacle *[m_number_of_obstacles];
        m_bdc_obstacles = new BoundaryDataController *[m_number_of_obstacles];
        //TODO(cvm) vector surfaces will be destroyed after this method therefore I need to preserve the data
        for (size_t i = 0; i < m_number_of_obstacles; i++) {
            m_obstacle_list[i] = ret_obstacles[i];
            m_bdc_obstacles[i] = bdc_obstacles[i];
        }
    } else {
        m_obstacle_list = new Obstacle *[m_number_of_obstacles];
        m_bdc_obstacles = new BoundaryDataController *[m_number_of_obstacles];
    }
#ifndef BENCHMARKING
    m_logger->debug("finished parsing obstacle parameter");
#endif
}

DomainController::~DomainController() {
    delete m_multigrid;
    delete m_bdc_boundary;
    for (size_t i = 0; i < m_number_of_obstacles; i++) {
        delete m_bdc_obstacles[i];
    }
}


DomainController *DomainController::getInstance(Settings::Settings const &settings) {
    if (singleton == nullptr) {
        singleton = new DomainController(settings);
    }
    return singleton;
}

// ================================= Printer =======================================================
// *************************************************************************************************
/// \brief  prints boundaries (outer, inner, surfaces)
// *************************************************************************************************
void DomainController::print_boundaries() {
#ifndef BENCHMARKING
    m_logger->info("-- Info summary");
    DomainData::getInstance()->print();
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
void DomainController::update_lists() {
    m_multigrid->update_lists();
}

// =========================================== Apply BCs ===========================================
// *************************************************************************************************
/// \brief  applies boundary conditions to all boundaries (domain, surfaces, obstacles)
/// \param  field   field
/// \param  sync    synchronization (default: false)
// *************************************************************************************************
void DomainController::apply_boundary(Field &field, bool sync) {
    m_multigrid->apply_boundary_condition(field, sync);
}

// ================================= detect neighbouring obstacles =================================
// *************************************************************************************************
/// \brief  handling of neighbouring obstacles, removes cells with circular constraints which never
/// change their value
// *************************************************************************************************
void DomainController::detect_neighbouring_obstacles() {
#ifndef BENCHMARKING
    m_logger->debug("start detecting neighbouring obstacles");
#endif
    for (size_t o1 = 0; o1 < m_number_of_obstacles; o1++) {
        Obstacle *obstacle1 = m_obstacle_list[o1];
        for (size_t o2 = o1 + 1; o2 < m_number_of_obstacles; o2++) {
            Obstacle *obstacle2 = m_obstacle_list[o2];
#ifndef BENCHMARKING
            m_logger->debug("scan for neighbouring cells for '{}' and '{}'", obstacle1->get_name(), obstacle2->get_name());
#endif
            Obstacle::remove_circular_constraints(obstacle1, obstacle2);
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("finished detecting neighbouring obstacles");
#endif
}

std::vector<FieldType> DomainController::get_used_fields() const {
    return m_bdc_boundary->get_used_fields();
}
