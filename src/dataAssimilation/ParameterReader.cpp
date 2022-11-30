/// \file      ParameterReader.cpp
/// \brief     read file with only field changes.
/// \date      Jan 03, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include <chrono>
#include "ParameterReader.h"
#include "DataAssimilation.h"
#include "../domain/DomainController.h"

return_parameter_reader ParameterReader::read_config(const std::string &file_name) {
    bool parameter_changes = false;
    try {
        m_logger->debug("parse file to string {}", file_name);
        auto file_content = Settings::parse_settings_from_file(file_name);
        m_logger->debug("parse document to XMLTree\n{}", file_content);
        tinyxml2::XMLDocument doc;
        doc.Parse(file_content.c_str());
        m_logger->debug("parse field changes");
        auto field_changes = Settings::parse_field_changes(doc.RootElement());
        auto da_methods = Settings::parse_data_assimilation_methods(doc.RootElement());
        for (std::tuple<std::string, std::string> tuple: da_methods) {
            std::string name = std::get<0>(tuple);
            std::string tag = std::get<1>(tuple);
            if (name == AssimilationMethods::temperature_source) {
                auto subsection = Settings::get_subsection(tag, doc.RootElement());
                parameter_changes = parameter_changes || temperature_source_changer(subsection, tag);
            } else if (name == AssimilationMethods::obstacle_changer) {
                auto subsection = Settings::get_subsection(tag, doc.RootElement());
                parameter_changes = parameter_changes || obstacle_changer(subsection, tag);
            } else {
                m_logger->debug("Unknown Assimilation method: {}", name);
            }
        }
        return {true, field_changes};
    } catch (const std::exception &ex) {
        std::cerr << ex.what() << std::endl;
    }
    Settings::data_assimilation::field_changes field_changes;
    field_changes.changed = false;
    return {parameter_changes, field_changes};
}

bool ParameterReader::temperature_source_changer(const tinyxml2::XMLElement *head, const std::string &context) {
    auto temperature_source = Settings::solver::parse_temperature_source(head, context);
    bool parameter_changes = true;
    // TODO (c++20)
    // bool parameter_changes = temperature_source != m_temperature_source;
    if (parameter_changes) {
        m_logger->debug("apply heat source changes");
        m_solver_controller.m_solver->replace_heat_source(temperature_source);
    }
    return parameter_changes;
}

bool ParameterReader::obstacle_changer(const tinyxml2::XMLElement *head, const std::string &context) {
    Settings::obstacles_parameters obstacle_parameters{ };
    for (const auto *i = head->FirstChildElement(); i; i = i->NextSiblingElement()) {
        obstacle_parameters.obstacles.emplace_back(Settings::parse_obstacle(i, context));
    }
    int counter_unmodified = 0;
    int counter_deleted = 0;
    int counter_new = 0;
    int counter_modified = 0;
    int counter_xml = 0;
    for (const auto &obstacle: obstacle_parameters.obstacles) {
        if (obstacle.state == State::UNMODIFIED) {
            counter_unmodified++;
        } else if (obstacle.state == State::DELETED) {
            counter_deleted++;
        } else if (obstacle.state == State::MODIFIED) {
            counter_modified++;
        } else if (obstacle.state == State::NEW) {
            counter_new++;
        } else if (obstacle.state == State::XML) {
            counter_xml++;
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("obstacles read in:\n  "
                    "unmodified: {}\n  "
                    "modified: {}\n  "
                    "deleted: {}\n  "
                    "new: {}\n  "
                    "XML: {}",
                    counter_unmodified,
                    counter_modified,
                    counter_deleted,
                    counter_new,
                    counter_xml);
#endif
    bool parameter_changes = counter_deleted + counter_new + counter_modified > 0;
    obstacle_parameters.enabled = counter_new + counter_modified + counter_unmodified + counter_xml > 0;
    if (parameter_changes) {
#ifndef BENCHMARKING
        m_logger->debug("apply obstacle changes");
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
#endif
        DomainController::getInstance()->replace_obstacles(obstacle_parameters);
#ifndef BENCHMARKING
        end = std::chrono::system_clock::now();
        long ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        m_logger->debug("replace obstacles time: {}", ms);
#endif
        m_solver_controller.update_sight();
#ifndef BENCHMARKING
        start = std::chrono::system_clock::now();
        m_solver_controller.m_solver->update_obstacle_change();
        end = std::chrono::system_clock::now();
        ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        m_logger->debug("update source: {}", ms);
    } else {
        m_logger->debug("no obstacle changes");
#endif
    }
    return parameter_changes;
}
