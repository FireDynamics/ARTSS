
#include "SurfaceSetting.h"

namespace Settings {
SurfaceSetting::SurfaceSetting(tinyxml2::XMLElement *xml_element) :
    id(std::atoi(xml_element->Attribute("ID"))),
    sx1(std::atof(xml_element->Attribute("sx1"))), sx2(std::atof(xml_element->Attribute("sx2"))),
    sy1(std::atof(xml_element->Attribute("sy1"))), sy2(std::atof(xml_element->Attribute("sy2"))),
    sz1(std::atof(xml_element->Attribute("sz1"))), sz2(std::atof(xml_element->Attribute("sz2"))) {
}

#ifndef BENCHMARKING
void SurfaceSetting::print(spdlog::logger logger) const {
    logger.debug("id={} sx=({}, {}) sy=({}, {}) oz=({}, {})",
            id,
            sx1, sx2,
            sy1, sy2,
            sz1, sz2);
}
#endif
}
