/// \file       TemperatureSourceChanger.cpp
/// \brief      Class for changing temperatur source
/// \date       Jan 05, 2021
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2022> Forschungszentrum Juelich All rights reserved.

#include "TemperatureSourceChanger.h"

Settings::data_assimilation::field_changes TemperatureSourceChanger::read_config(const std::string &filename) {
    auto file_content = Settings::parse_settings_from_file(filename);
    tinyxml2::XMLDocument doc;
    doc.Parse(file_content.c_str());
    auto temperature_source = Settings::solver::parse_temperature_source(doc.RootElement(), "temperature");
    m_solver_controller.m_solver->replace_heat_source(temperature_source);
    auto changes = Settings::parse_field_changes(doc.RootElement(), "field_changes");
    return changes;
}
