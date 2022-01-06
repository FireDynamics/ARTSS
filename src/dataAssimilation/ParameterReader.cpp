/// \file      ParameterReader.cpp
/// \brief
/// \date      Jan 03, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include "ParameterReader.h"

Settings::data_assimilation::field_changes ParameterReader::read_config(const std::string &filename) {
    auto file_content = Settings::parse_settings_from_file(filename);
    auto root = Settings::parse_file_content(file_content);
    return Settings::parse_field_changes(root, "field_changes");
}
