//
// Created by linh on 05.01.22.
//

#include "TemperatureSourceChanger.h"

Settings::data_assimilation::field_changes TemperatureSourceChanger::read_config(const std::string &filename) {
    auto file_content = Settings::parse_settings_from_file(filename);
    auto root = Settings::parse_file_content(file_content);
    auto temperature_source = Settings::parse_temperature_source(root, "temperature");
    // TODO  replace temperature source
    auto changes = Settings::parse_field_changes(root, "field_changes");
    return changes;
}
