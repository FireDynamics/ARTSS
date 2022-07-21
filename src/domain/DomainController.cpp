/// \file       DomainController.cpp
/// \brief      Controller class for domain, grants access to indices lists
/// \date       Oct 01, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "DomainController.h"
#include <string>

std::unique_ptr<DomainController> DomainController::single{};

DomainController::DomainController(const Settings::Settings &settings) :
        m_settings(settings) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    auto [bdc_domain, surfaces, bdc_surfaces, obstacles, bdc_obstacles] = read_XML();
    size_t multigrid_level = DomainData::getInstance()->get_levels();
    m_multigrid = new Multigrid(surfaces, bdc_surfaces,
                                obstacles, bdc_obstacles,
                                bdc_domain,
                                multigrid_level);
#ifndef BENCHMARKING
    print_boundaries(bdc_domain, bdc_surfaces, bdc_obstacles);
#endif
}

// ================================= Read XML ======================================================
// *************************************************************************************************
/// \brief  Reads in all parameters of boundary, obstacles and surfaces
// *************************************************************************************************
return_xml_objects DomainController::read_XML() {
#ifndef BENCHMARKING
    m_logger->debug("start parsing XML");
    m_logger->debug("start parsing boundary parameter");
#endif
   BoundaryDataController bdc_domain(m_settings.boundary_parameters.boundaries);
#ifndef BENCHMARKING
    m_logger->debug("finished parsing boundary parameter");
#endif
    auto [obstacles, bdc_obstacles] = parse_obstacle_parameter(m_settings.obstacles_parameters);
    detect_neighbouring_obstacles(obstacles);
    auto [surfaces, bdc_surfaces] = parse_surface_parameter(m_settings.surfaces_parameters);
#ifndef BENCHMARKING
    m_logger->debug("finished parsing XML");
#endif
    return {bdc_domain, surfaces, bdc_surfaces, obstacles, bdc_obstacles};
}

// ================================= Parser ========================================================
// *************************************************************************************************
/// \brief  parses surfaces from XML file
/// \param  xmlParameter pointer to XMLElement to start with
// *************************************************************************************************
return_surface DomainController::parse_surface_parameter(const Settings::surfaces_parameters &surface_settings) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing surface parameter");
#endif
    // SURFACES
    // TODO(issue 5): surfaces
    std::vector<Surface> surfaces;
    std::vector<BoundaryDataController> bdc_surfaces;
    if (surface_settings.enabled) {
        surfaces.reserve(surface_settings.surfaces.size());
        bdc_surfaces.reserve(surface_settings.surfaces.size());
        for (const auto &surface: surface_settings.surfaces) {
            surfaces.emplace_back(surface.start_coords, surface.end_coords, surface.name, surface.patch);
            bdc_surfaces.emplace_back(surface.boundaries);
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
return_obstacle DomainController::parse_obstacle_parameter(const Settings::obstacles_parameters &obstacle_settings) {
#ifndef BENCHMARKING
    m_logger->debug("start parsing obstacle parameter");
#endif
    std::vector<Obstacle> obstacles;
    std::vector<BoundaryDataController> bdc_obstacles;
    if (obstacle_settings.enabled) {
       obstacles.reserve(obstacle_settings.obstacles.size());
       bdc_obstacles.reserve(obstacle_settings.obstacles.size());
       for (const Settings::obstacle &o: obstacle_settings.obstacles) {
           obstacles.emplace_back(o.start_coords, o.end_coords, o.name);
           bdc_obstacles.emplace_back(o.boundaries);
       }
    }
#ifndef BENCHMARKING
    m_logger->debug("finished parsing obstacle parameter");
#endif
    return {obstacles, bdc_obstacles};
}

DomainController::~DomainController() {
    delete m_multigrid;
}

// ================================= Printer =======================================================
// *************************************************************************************************
/// \brief  prints boundaries (outer, inner, surfaces)
// *************************************************************************************************
void DomainController::print_boundaries(const BoundaryDataController &bdc_domain, const std::vector<BoundaryDataController> &bdc_surfaces, const std::vector<BoundaryDataController> &bdc_obstacles) {
#ifndef BENCHMARKING
    m_logger->info("-- Info summary");
    DomainData::getInstance()->print();
    bdc_domain.print();
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
void DomainController::detect_neighbouring_obstacles(std::vector<Obstacle> &obstacle_list) {
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

std::vector<FieldType> DomainController::get_used_fields() const {
    return m_multigrid->get_used_fields();
}

bool DomainController::is_blocked_by_obstacle(Coordinate<size_t> &start, Coordinate<size_t> &end) {
    return m_multigrid->is_blocked_by_obstacle(start, end);
}
