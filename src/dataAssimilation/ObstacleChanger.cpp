/// \file      ObstacleChanger.cpp
/// \brief
/// \date      Nov 23, 2022
/// \author    My Linh Wuerzburger
/// \copyright <2015-2022> Forschungszentrum Juelich. All rights reserved

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
        bool parameter_changes = true;
        // TODO (c++20)
        // bool parameter_changes = new struct obstacle_parameters != old struct obstacle_parameters;
        if (parameter_changes) {
            m_logger->debug("apply heat source changes");
            DomainController::getInstance()->replace_obstacles(obstacle_parameters);
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
