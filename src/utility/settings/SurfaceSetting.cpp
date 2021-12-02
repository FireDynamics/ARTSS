/// \file       SurfaceSetting.cpp
/// \brief      Model of Surfaces in the XML
/// \date       Dec 01, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.

#include "SurfaceSetting.h"

namespace Settings {
SurfaceSetting::SurfaceSetting(tinyxml2::XMLElement *xml_element) :
    name(xml_element->Attribute("name")) {
    for (auto i = xml_element->FirstChildElement(); i; i = i->NextSiblingElement()) {
        if (i->Name() == std::string("geometry")) {
            patch = i->Attribute("patch");
            sx1 = std::atof(i->Attribute("sx1"));
            sx2 = std::atof(i->Attribute("sx2"));
            sy1 = std::atof(i->Attribute("sy1"));
            sy2 = std::atof(i->Attribute("sy2"));
            sz1 = std::atof(i->Attribute("sz1"));
            sz2 = std::atof(i->Attribute("sz2"));
        } else if (i->Name() == std::string("boundary")) {
            boundaries.push_back(BoundarySetting(i));
        }
    }
}

#ifndef BENCHMARKING
void SurfaceSetting::print(spdlog::logger logger) const {
    logger.debug("name={} patch={} sx=({}, {}) sy=({}, {}) oz=({}, {})",
            name,
            patch,
            sx1, sx2,
            sy1, sy2,
            sz1, sz2);
    for (const auto &i: boundaries) {
        i.print(logger);
    }
}
#endif
}
