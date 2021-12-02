/// \file       SurfaceSetting.h
/// \brief      Model of Surfaces in the XML
/// \date       Dec 01, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_SETTINGS_SURFACESETTING_H
#define ARTSS_UTILITY_SETTINGS_SURFACESETTING_H

#include <string>

#include "BoundarySetting.h"
#include "tinyxml2.h"
#include "../GlobalMacrosTypes.h"

#ifndef BENCHMARKING
#include <spdlog/logger.h>
#endif


namespace Settings {
class SurfaceSetting {
 public:
    explicit SurfaceSetting(tinyxml2::XMLElement *xml_element);

    std::string get_name() const { return name; }
    real get_sx1() const { return sx1; }
    real get_sx2() const { return sx2; }
    real get_sy1() const { return sy1; }
    real get_sy2() const { return sy2; }
    real get_sz1() const { return sz1; }
    real get_sz2() const { return sz2; }

#ifndef BENCHMARKING
    void print(spdlog::logger logger) const;
#endif

    std::vector<BoundarySetting> get_boundaries() const { return boundaries; }

 private:
    std::string name;
    std::string patch;
    real sx1, sx2;
    real sy1, sy2;
    real sz1, sz2;

    std::vector<BoundarySetting> boundaries;
};
}

#endif
