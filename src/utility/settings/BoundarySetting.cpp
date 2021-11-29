
#include "BoundarySetting.h"

namespace Settings {
BoundarySetting::BoundarySetting(tinyxml2::XMLElement *xml_element) :
    field(xml_element->Attribute("field")),
    patch(xml_element->Attribute("patch")),
    type(xml_element->Attribute("type")),
    value(std::atof(xml_element->Attribute("value"))) {
}

#ifndef BENCHMARKING
void BoundarySetting::print(spdlog::logger logger) const {
    logger.debug("field={} patch={} type={} value={}", field, patch, type, value);
}
#endif
}
