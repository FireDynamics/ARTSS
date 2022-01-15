/// \file      ParameterReader.cpp
/// \brief     read file with only field changes.
/// \date      Jan 03, 2021
/// \author    My Linh Wuerzburger
/// \copyright <2015-2021> Forschungszentrum Juelich. All rights reserved

#include "ParameterReader.h"

ParameterReader::ParameterReader() {
    m_logger = Utility::create_logger(typeid(this).name());
}

return_parameter_reader ParameterReader::read_config(const std::string &filename) {
    try {
        m_logger->debug("parse file to string");
        auto file_content = Settings::parse_settings_from_file(filename);
        m_logger->debug("parse document to XMLTree");
        tinyxml2::XMLDocument doc;
        doc.Parse(file_content.c_str());
        m_logger->debug("parse field changes");
        auto field_changes = Settings::parse_field_changes(doc.RootElement(), "field_changes");
        return {true, field_changes};
    } catch (const std::exception &ex) {
        std::cerr << ex.what() << std::endl;
    }
    Settings::data_assimilation::field_changes field_changes;
    field_changes.changed = false;
    return {false, field_changes} ;
}
