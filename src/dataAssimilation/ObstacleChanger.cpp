/// \file      ObstacleChanger.cpp
/// \brief
/// \date      Nov 23, 2022
/// \author    My Linh Wuerzburger
/// \copyright <2015-2022> Forschungszentrum Juelich. All rights reserved

#include <chrono>
#include "ObstacleChanger.h"
#include "../domain/DomainController.h"

return_parameter_reader ObstacleChanger::read_config(const std::string &filename) {
    try {
        m_logger->debug("parse file to string {}", filename);
        auto file_content = Settings::parse_settings_from_file(filename);
        m_logger->debug("parse document from {} to XMLTree {}", filename, file_content);
        tinyxml2::XMLDocument doc;
        doc.Parse(file_content.c_str());
        m_logger->debug("parse obstacles {}", static_cast<void *>(doc.RootElement()));
        auto obstacle_parameters = Settings::parse_obstacles_parameters(doc.RootElement());
        bool parameter_changes = false;
        for (const auto& obstacle: obstacle_parameters.obstacles) {
            if (obstacle.state == State::MODIFIED || obstacle.state == State::NEW || obstacle.state == State::DELETED) {
                parameter_changes = true;
                break;
            }
        }
        if (parameter_changes) {
            m_logger->debug("apply obstacle changes");
            std::chrono::time_point<std::chrono::system_clock> start, end;
            start = std::chrono::system_clock::now();
            DomainController::getInstance()->replace_obstacles(obstacle_parameters);
            end = std::chrono::system_clock::now();
            long ms = std::chrono::duration_cast < std::chrono::milliseconds > (end - start).count();
            m_logger->debug("replace obstacles time: {}", ms);
#ifndef BENCHMARKING
            int counter_unmodified = 0;
            int counter_deleted = 0;
            int counter_rest = 0;
            for (const auto& obstacle: obstacle_parameters.obstacles) {
                if (obstacle.state == State::UNMODIFIED) {
                    counter_unmodified++;
                } else if (obstacle.state == State::DELETED) {
                    counter_deleted++;
                } else if (obstacle.state == State::MODIFIED || obstacle.state == State::NEW) {
                    counter_rest++;
                } else {
                    m_logger->warn("wrong state: {} ({})", obstacle.state, obstacle.name);
                }
            }
            m_logger->debug("changed obstacles, unmodified: {}, modified/new: {}, deleted {}", counter_unmodified, counter_deleted, counter_rest);
#endif
            m_solver_controller.update_sight();
            start = std::chrono::system_clock::now();
            m_solver_controller.m_solver->update_obstacle_change();
            end = std::chrono::system_clock::now();
            ms = std::chrono::duration_cast < std::chrono::milliseconds > (end - start).count();
            m_logger->debug("update source: {}", ms);
#ifndef BENCHMARKING
        } else {
            m_logger->debug("no obstacle changes");
#endif
        }
        m_logger->debug("parse field changes");

        auto field_changes = Settings::parse_field_changes(doc.RootElement(), "TemperatureSourceChanger");
        return {parameter_changes, field_changes};
    } catch (const std::exception &ex) {
        std::cerr << ex.what() << std::endl;
    }
    Settings::data_assimilation::field_changes field_changes;
    field_changes.changed = false;
    return {false, field_changes};
}
