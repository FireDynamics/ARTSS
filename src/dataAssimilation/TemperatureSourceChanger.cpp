/// \file       TemperatureSourceChanger.cpp
/// \brief      Class for changing temperature source
/// \date       Jan 05, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2022> Forschungszentrum Juelich All rights reserved.

#include "TemperatureSourceChanger.h"
#include "../solver/SolverSelection.h"

return_parameter_reader TemperatureSourceChanger::read_config(const std::string &filename) {
    try {
        m_logger->debug("parse file to string");
        auto file_content = Settings::parse_settings_from_file(filename);
        m_logger->debug("parse document from {} to XMLTree {}", filename, file_content);
        tinyxml2::XMLDocument doc;
        doc.Parse(file_content.c_str());
        m_logger->debug("parse heat source changes {}", static_cast<void *>(doc.RootElement()));
        auto temperature_source = Settings::solver::parse_temperature_source(doc.RootElement(), "temperature");
        bool parameter_changes = true;
        // TODO
        // bool parameter_changes = temperature_source != m_temperature_source;
        if (parameter_changes) {
            m_logger->debug("apply heat source changes");
            m_solver_controller.m_solver->replace_heat_source(temperature_source);
        }
        m_logger->debug("parse field changes");
        auto field_changes = Settings::parse_field_changes(doc.RootElement(), "TemperatureSourceChanger");
        return {parameter_changes, field_changes};
    } catch (const std::exception &ex) {
        std::cerr << ex.what() << std::endl;
    }
    Settings::data_assimilation::field_changes field_changes;
    field_changes.changed = false;
    return {false, field_changes} ;
}

TemperatureSourceChanger::TemperatureSourceChanger(const SolverController &solver_controller,
                                                   const Settings::solver::temperature_source &temperature_source) :
        m_solver_controller(solver_controller),
        m_temperature_source(temperature_source) {
    m_logger = Utility::create_logger(typeid(this).name());
}
