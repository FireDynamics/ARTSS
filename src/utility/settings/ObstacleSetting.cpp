
#include "ObstacleSetting.h"

namespace Settings {
ObstacleSetting::ObstacleSetting(tinyxml2::XMLElement *xml_element) :
    name(xml_element->Attribute("name")) {
    for (auto i=xml_element->FirstChildElement(); i; i = i->NextSiblingElement()) {
        if (i->Name() == std::string("geometry")) {
            ox1 = std::atof(i->Attribute("ox1"));
            ox2 = std::atof(i->Attribute("ox2"));
            oy1 = std::atof(i->Attribute("oy1"));
            oy2 = std::atof(i->Attribute("oy2"));
            oz1 = std::atof(i->Attribute("oz1"));
            oz2 = std::atof(i->Attribute("oz2"));
        } else if (i->Name() == std::string("boundary")) {
            boundaries.push_back(BoundarySetting(i));
        }
    }
}

#ifndef BENCHMARKING
void ObstacleSetting::print(spdlog::logger logger) const {
    logger.debug("ox=({}, {}) oy=({}, {}) oz=({}, {})", ox1, ox2, oy1, oy2, oz1, oz2);
    for(auto i : boundaries) {
        i.print(logger);
    }
}
#endif
}
