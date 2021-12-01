/// \file       ObstacleSetting.h
/// \brief      Model of Obstacles in the XML
/// \date       Dec 01, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_SETTINGS_OBSTACLESETTING_H
#define ARTSS_UTILITY_SETTINGS_OBSTACLESETTING_H

#include <string>
#include <vector>

#include "tinyxml2.h"
#include "BoundarySetting.h"
#include "../GlobalMacrosTypes.h"

#ifndef BENCHMARKING
#include <spdlog/logger.h>
#endif


namespace Settings {
class ObstacleSetting {
 public:
    explicit ObstacleSetting(tinyxml2::XMLElement *xml_element);

    std::string get_name() const { return name; }
    real get_ox1() const { return ox1; }
    real get_ox2() const { return ox2; }
    real get_oy1() const { return oy1; }
    real get_oy2() const { return oy2; }
    real get_oz1() const { return oz1; }
    real get_oz2() const { return oz2; }

#ifndef BENCHMARKING
    void print(spdlog::logger logger) const;
#endif

    std::vector<BoundarySetting> get_boundaries() const { return boundaries; }

 private:
    std::string name;
    real ox1, ox2;
    real oy1, oy2;
    real oz1, oz2;
    std::vector<BoundarySetting> boundaries;
};
}

#endif
