/// \file       BoundaryController.cpp
/// \brief      Controller class for boundary
/// \date       Oct 01, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "BoundaryController.h"
#include <string>

BoundaryController *BoundaryController::singleton = nullptr;  // Singleton


BoundaryController::BoundaryController(Settings::Settings const &settings) :
        m_settings(settings) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    m_bdc_boundary = new BoundaryDataController();
    auto [surfaces, bdc_surfaces] = read_XML();
    size_t multigrid_level = DomainData::getInstance()->get_levels();
    m_multigrid = new Multigrid(surfaces, bdc_surfaces,
                                m_number_of_obstacles, m_obstacle_list,
                                m_bdc_boundary, m_bdc_obstacles,
                                multigrid_level);
#ifndef BENCHMARKING
    print_boundaries();
#endif
}

// ================================= Read XML ======================================================
// *************************************************************************************************
/// \brief  Reads in all parameters of boundary, obstacles and surfaces
// *************************************************************************************************
return_surface BoundaryController::read_XML() {
#ifndef BENCHMARKING
    m_logger->debug("start parsing XML");
#endif

    parse_boundary_parameter(m_settings.get_boundaries());
    parse_obstacle_parameter(m_settings.get_obstacles());
    detect_neighbouring_obstacles();
    return_surface ret_surfaces = parse_surface_parameter(m_settings.get_surfaces());
#ifndef BENCHMARKING
    m_logger->debug("finished parsing XML");
#endif
    return ret_surfaces;
}

// ================================= Parser ========================================================
// *************************************************************************************************
/// \brief  parses boundaries of domain from XML file
/// \param  xmlParameter pointer to XMLElement to start with
// *************************************************************************************************
void BoundaryController::parse_boundary_parameter(const std::vector<Settings::BoundarySetting> &boundaries) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing boundary parameter");
#endif
// BOUNDARY
    for (const auto &boundary: boundaries) {
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
return_surface BoundaryController::parse_surface_parameter(const std::vector<Settings::SurfaceSetting> &surface_setting) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing surface parameter");
#endif
    // SURFACES
    // TODO(issue 5): surfaces
    m_has_surfaces = m_settings.get_bool("surfaces/enabled");
    std::vector<Surface> surfaces;
    std::vector<BoundaryDataController> bdc_surfaces;
    if (m_has_surfaces) {
        surfaces.reserve(surfaces.size());
        bdc_surfaces.reserve(surfaces.size());
        for (const auto &surface: surface_setting) {
            std::string name = surface.get_name();
            Patch patch = PatchObject::match_patch(surface.get_patch());
            real sx1 = surface.get_sx1();
            real sx2 = surface.get_sx2();
            real sy1 = surface.get_sy1();
            real sy2 = surface.get_sy2();
            real sz1 = surface.get_sz1();
            real sz2 = surface.get_sz2();
            surfaces.emplace_back(sx1, sx2, sy1, sy2, sz1, sz2, name, patch);

            BoundaryDataController bdc;
            for (const auto &boundary: surface.get_boundaries()) {
                bdc.add_boundary_data(boundary);
            }
            bdc_surfaces.emplace_back(bdc);
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("finished parsing surface parameter");
#endif
    surfaces.shrink_to_fit();
    bdc_surfaces.shrink_to_fit();
    return {surfaces, bdc_surfaces};
}

// ================================= Parser ========================================================
// *************************************************************************************************
/// \brief  parses obstacles from XML file
/// \param  xmlParameter pointer to XMLElement to start with
// *************************************************************************************************
void BoundaryController::parse_obstacle_parameter(const std::vector<Settings::ObstacleSetting> &obstacles) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing obstacle parameter");
#endif
// OBSTACLES
    m_has_obstacles = m_settings.get_bool("obstacles/enabled");
    if (m_has_obstacles) {
        std::vector<Obstacle *> ret_obstacles;
        std::vector<BoundaryDataController *> bdc_obstacles;

        for (const Settings::ObstacleSetting &obstacle: obstacles) {
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

            for (const auto &bound: obstacle.get_boundaries()) {
                bdc->add_boundary_data(bound);
            }

            auto o = new Obstacle(ox1, ox2, oy1, oy2, oz1, oz2, name);
            ret_obstacles.emplace_back(o);
            bdc_obstacles.emplace_back(bdc);
        }

        m_number_of_obstacles = ret_obstacles.size();
        m_obstacle_list = new Obstacle *[m_number_of_obstacles];
        m_bdc_obstacles = new BoundaryDataController *[m_number_of_obstacles];
        for (size_t i = 0; i < m_number_of_obstacles; i++) {
            *(m_obstacle_list + i) = ret_obstacles[i];
            *(m_bdc_obstacles + i) = bdc_obstacles[i];
        }
    } else {
        m_obstacle_list = new Obstacle *[m_number_of_obstacles];
        m_bdc_obstacles = new BoundaryDataController *[m_number_of_obstacles];
    }
#ifndef BENCHMARKING
    m_logger->debug("finished parsing obstacle parameter");
#endif
}

BoundaryController::~BoundaryController() {
    delete m_multigrid;
}


BoundaryController *BoundaryController::getInstance(Settings::Settings const &settings) {
    if (singleton == nullptr) {
        singleton = new BoundaryController(settings);
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
    DomainData::getInstance()->print();
    m_bdc_boundary->print();
    for (size_t i = 0; i < m_number_of_obstacles; i++) {
        m_obstacle_list[i]->print();
        m_bdc_obstacles[i]->print();
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

// =========================================== Apply BCs ===========================================
// *************************************************************************************************
/// \brief  applies boundary conditions to all boundaries (domain, surfaces, obstacles)
/// \param  field   field
/// \param  sync    synchronization (default: false)
// *************************************************************************************************
void BoundaryController::apply_boundary(Field &field, bool sync) {
    m_multigrid->apply_boundary_condition(field, sync);
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
#ifndef BENCHMARKING
            m_logger->debug("scan for neighbouring cells for '{}' and '{}'", obstacle1->get_name(),
                            obstacle2->get_name());
#endif
            Obstacle::remove_circular_constraints(obstacle1, obstacle2);
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("finished detecting neighbouring obstacles");
#endif
}

std::vector<FieldType> BoundaryController::get_used_fields() const {
    return m_bdc_boundary->get_used_fields();
}
