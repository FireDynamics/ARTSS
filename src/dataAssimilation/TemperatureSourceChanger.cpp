//
// Created by linh on 05.01.22.
//

#include "TemperatureSourceChanger.h"

std::vector<FieldType> TemperatureSourceChanger::read_config(const std::string &filename) {
    auto file_content = Settings::parse_settings_from_file(filename);
    auto root = Settings::parse_file_content(file_content);
    auto temperature_source = Settings::parse_temperature_source(root, "temperature");

    auto changes = Settings::parse_field_changes(root, "field_changes");
    std::vector<FieldType> fields;
    return fields;
}
