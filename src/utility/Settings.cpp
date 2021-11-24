
#include "Settings.h"

#include "Utility.h"

Settings::Settings(std::string path) :
        filename(path) {
    tinyxml2::XMLDocument doc;
    doc.LoadFile(path.c_str());

    for (auto i = doc.RootElement()->FirstChildElement(); i; i = i->NextSiblingElement()) {
        std::cout << i->Name();
        if (i->Name() == std::string("boundaries")) {
            read_boundaries(i);
        } else if (i->Name() == std::string("obstacles")) {
            read_obstacles(i);
        } else if (i->Name() == std::string("surfaces")) {
            read_surfaces(i);
        } else {
            read_config("", i);
        }
    }

    m_logger = Utility::create_logger(*this, "XMLFile");
    m_logger->debug("start the simulation of \"{}\"", path);
}


void Settings::read_boundaries(tinyxml2::XMLElement *elem) {
    for (auto i = elem->FirstChildElement(); i; i = i->NextSiblingElement()) {
        m_boundaries.push_back(BoundarySetting(i));
    }
}

void Settings::read_obstacles(tinyxml2::XMLElement *elem) {
    for (auto i = elem->FirstChildElement(); i; i = i->NextSiblingElement()) {
        m_obstacles.push_back(ObstacleSetting(i));
    }
}

void Settings::read_surfaces(tinyxml2::XMLElement *elem) {
    for (auto i = elem->FirstChildElement(); i; i = i->NextSiblingElement()) {
        m_surfaces.push_back(SurfaceSetting(i));
    }
}


void Settings::read_config(std::string prefix, tinyxml2::XMLElement *elem) {
    prefix += elem->Name();
    prefix += "/";
    for (auto i = elem->FirstAttribute(); i; i = i->Next()) {
        sset(prefix + i->Name(), i->Value());
    }

    for (auto i = elem->FirstChildElement(); i; i = i->NextSiblingElement()) {
        if (i->GetText()) {
            sset(prefix + i->Name(), i->GetText());
        }

        read_config(prefix, i);
    }
}

void Settings::print_config() const {
#ifndef BENCHMARKING
    std::map<std::string, std::string> ordered(m_proxy.begin(), m_proxy.end());
    for(auto i = ordered.begin(); i != ordered.end(); ++i) {
        m_logger->debug("{} = \"{}\"", i->first, i->second);
    }

    for(auto i : get_boundaries()) {
        i.print(*m_logger);
    }

    for(auto i : get_obstacles()) {
        i.print(*m_logger);
    }

    for(auto i : get_surfaces()) {
        i.print(*m_logger);
    }
#endif
}

std::string Settings::sget(std::string path) const {
    auto iter = m_proxy.find(path);
    return iter->second;
}

std::string Settings::get(std::string path) const {
    auto iter = m_proxy.find(path);
    if (iter == m_proxy.end()) {
#ifndef BENCHMARKING
        m_logger->error("didn't found \"{}\" in settings", path);
        print_config();
        throw std::invalid_argument("\"" + path + "\" is an unkown setting");
#endif
    }

#ifndef BENCHMARKING
     m_logger->debug("reading: \"{}\" with value: \"{}\"", path, iter->second);
#endif

    return iter->second;
}

BoundarySetting::BoundarySetting(tinyxml2::XMLElement *xml_element) :
    field(xml_element->Attribute("field")),
    patch(xml_element->Attribute("patch")),
    type(xml_element->Attribute("type")),
    value(std::atof(xml_element->Attribute("value"))) {
        assert(xml_element->Name() == std::string("boundary"));
}

void BoundarySetting::print(spdlog::logger logger) const {
    logger.debug("field={} patch={} type={} value={}", field, patch, type, value);
}

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
        } else {
            assert(1 == 0);
        }
    }
}

void ObstacleSetting::print(spdlog::logger logger) const {
    logger.debug("ox=({}, {}) oy=({}, {}) oz=({}, {})", ox1, ox2, oy1, oy2, oz1, oz2);
    for(auto i : boundaries) {
        i.print(logger);
    }
}

SurfaceSetting::SurfaceSetting(tinyxml2::XMLElement *xml_element) :
    id(std::atoi(xml_element->Attribute("ID"))),
    sx1(std::atof(xml_element->Attribute("sx1"))), sx2(std::atof(xml_element->Attribute("sx2"))),
    sy1(std::atof(xml_element->Attribute("sy1"))), sy2(std::atof(xml_element->Attribute("sy2"))),
    sz1(std::atof(xml_element->Attribute("sz1"))), sz2(std::atof(xml_element->Attribute("sz2"))) {
        assert(xml_element->Name() == std::string("surface"));
}

void SurfaceSetting::print(spdlog::logger logger) const {
    logger.debug("id={} sx=({}, {}) sy=({}, {}) oz=({}, {})",
            id,
            sx1, sx2,
            sy1, sy2,
            sz1, sz2);
}
