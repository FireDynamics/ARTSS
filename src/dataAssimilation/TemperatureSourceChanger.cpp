/// \file       TemperatureSourceChanger.cpp
/// \brief      Class for changing temperature source
/// \date       Jan 05, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2022> Forschungszentrum Juelich All rights reserved.

#include "TemperatureSourceChanger.h"

Settings::data_assimilation::field_changes TemperatureSourceChanger::read_config(const std::string &filename, const real t_cur) {
    m_logger->debug("parse file to string");
    auto file_content = Settings::parse_settings_from_file(filename);
    m_logger->debug("parse document to XMLTree {}", file_content);
    tinyxml2::XMLDocument doc;
    doc.Parse(file_content.c_str());
    m_logger->debug("parse heat source changes {}", static_cast<void *>(doc.RootElement()));
    // TODO maybe replace the required section with optional with old values, see m_temperature_source
    auto temperature_source = Settings::solver::parse_temperature_source(doc.RootElement(), "temperature");
    m_logger->debug("apply heat source changes");
    m_solver_controller.m_solver->replace_heat_source(temperature_source, t_cur);
    m_logger->debug("parse field changes");
    auto changes = Settings::parse_field_changes(doc.RootElement(), "field_changes");
    return changes;
}

TemperatureSourceChanger::TemperatureSourceChanger(const SolverController &solver_controller,
                                                   const Settings::solver::temperature_source &temperature_source) :
        m_solver_controller(solver_controller),
        m_temperature_source(temperature_source) {
    m_logger = Utility::create_logger(typeid(this).name());
}
