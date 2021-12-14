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
    auto [surfaces, bdc_surfaces, obstacles, bdc_obstacles] = read_XML();
    size_t multigrid_level = DomainData::getInstance()->get_levels();
    m_multigrid = new Multigrid(surfaces, bdc_surfaces,
                                obstacles, bdc_obstacles,
                                m_bdc_boundary,
                                multigrid_level);
#ifndef BENCHMARKING
    print_boundaries(bdc_surfaces, bdc_obstacles);
#endif
}

// ================================= Read XML ======================================================
// *************************************************************************************************
/// \brief  Reads in all parameters of boundary, obstacles and surfaces
// *************************************************************************************************
return_xml_objects BoundaryController::read_XML() {
#ifndef BENCHMARKING
    m_logger->debug("start parsing XML");
#endif

    parse_boundary_parameter(m_settings.get_boundaries());
    auto [obstacles, bdc_obstacles] = parse_obstacle_parameter(m_settings.get_obstacles());
    detect_neighbouring_obstacles(obstacles);
    auto [surfaces, bdc_surfaces] = parse_surface_parameter(m_settings.get_surfaces());
#ifndef BENCHMARKING
    m_logger->debug("finished parsing XML");
#endif
    return {surfaces, bdc_surfaces, obstacles, bdc_obstacles};
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
    m_bdc_boundary = new BoundaryDataController(boundaries);
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
            bdc_surfaces.emplace_back(surface.get_boundaries());
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
return_obstacle BoundaryController::parse_obstacle_parameter(const std::vector<Settings::ObstacleSetting> &obstacle_setting) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing obstacle parameter");
#endif
// OBSTACLES
    m_has_obstacles = m_settings.get_bool("obstacles/enabled");
    std::vector<Obstacle> obstacles;
    std::vector<BoundaryDataController> bdc_obstacles;
    if (m_has_obstacles) {
        obstacles.reserve(obstacle_setting.size());
        bdc_obstacles.reserve(obstacle_setting.size());
        for (const Settings::ObstacleSetting &obstacle: obstacle_setting) {
            std::string name = obstacle.get_name();
#ifndef BENCHMARKING
            m_logger->debug("read obstacle '{}'", name);
#endif
            real ox1 = obstacle.get_ox1();
            real ox2 = obstacle.get_ox2();
            real oy1 = obstacle.get_oy1();
            real oy2 = obstacle.get_oy2();
            real oz1 = obstacle.get_oz1();
            real oz2 = obstacle.get_oz2();

            obstacles.emplace_back(ox1, ox2, oy1, oy2, oz1, oz2, name);
            bdc_obstacles.emplace_back(obstacle.get_boundaries());
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("finished parsing obstacle parameter");
#endif
    return {obstacles, bdc_obstacles};
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
void BoundaryController::print_boundaries(const std::vector<BoundaryDataController> &bdc_surfaces, const std::vector<BoundaryDataController> &bdc_obstacles) {
#ifndef BENCHMARKING
    m_logger->info("-- Info summary");
    DomainData::getInstance()->print();
    m_bdc_boundary->print();
    for (const auto & bdc_obstacle : bdc_obstacles) {
        bdc_obstacle.print();
    }
    for (const auto & bdc_surface : bdc_surfaces) {
        bdc_surface.print();
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
void BoundaryController::detect_neighbouring_obstacles(std::vector<Obstacle> &obstacle_list) {
#ifndef BENCHMARKING
    m_logger->debug("start detecting neighbouring obstacles");
#endif
    size_t number_of_obstacles = obstacle_list.size();
    for (size_t id1 = 0; id1 < number_of_obstacles; id1++) {
        Obstacle &obstacle1 = obstacle_list[id1];
        for (size_t id2 = id1 + 1; id2 < number_of_obstacles; id2++) {
            Obstacle &obstacle2 = obstacle_list[id2];
#ifndef BENCHMARKING
            m_logger->debug("scan for neighbouring cells for '{}' and '{}'", obstacle1.get_name(),
                            obstacle2.get_name());
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
