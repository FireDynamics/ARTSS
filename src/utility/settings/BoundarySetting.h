/// \file       BoundarySetting.h
/// \brief      Model of Boundaries in the XML
/// \date       Dec 01, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_SETTINGS_BOUNDARYSETTING_H
#define ARTSS_UTILITY_SETTINGS_BOUNDARYSETTING_H

#include <string>

#include "tinyxml2.h"
#include "../GlobalMacrosTypes.h"

#ifndef BENCHMARKING
#include <spdlog/logger.h>
#endif


namespace Settings {
class BoundarySetting {
 public:
    explicit BoundarySetting(tinyxml2::XMLElement *xml_element);

    std::string get_field() const { return field; }
    std::string get_patch() const { return patch; }
    std::string get_type() const { return type; }
    real get_value() const { return value; }

#ifndef BENCHMARKING
    void print(spdlog::logger logger) const;
#endif

 private:
    std::string field;
    std::string patch;
    std::string type;
    real value;
};
}

#endif
