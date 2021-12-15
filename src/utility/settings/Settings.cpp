/// \file       Settings.cpp
/// \brief      Access parameters of XML file
/// \date       Nov 17, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Settings.h"

#include "../Utility.h"

namespace Settings {
Settings::Settings(const std::string &path) :
        filename(path) {
    tinyxml2::XMLDocument doc;
    doc.LoadFile(path.c_str());

    for (auto i = doc.RootElement()->FirstChildElement(); i; i = i->NextSiblingElement()) {
        if (i->Name() == std::string("boundaries")) {
            read_boundaries(i);
        } else if (i->Name() == std::string("obstacles")) {
            read_obstacles(i);
        } else if (i->Name() == std::string("surfaces")) {
            read_surfaces(i);
        }

        read_config("", i);
    }

#ifndef BENCHMARKING
    Utility::create_logger(sget("logging/level"), sget("logging/file"));  // create global logger
    m_logger = Utility::create_logger("XMLFile");
    m_logger->debug("start the simulation of \"{}\"", path);
    print_config();
#endif
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
    for(auto &i : ordered) {
        m_logger->debug(R"({} = "{}")", i.first, i.second);
    }

    for(const auto &i : get_boundaries()) {
        i.print(*m_logger);
    }

    for(const auto &i : get_obstacles()) {
        i.print(*m_logger);
    }

    for(const auto &i : get_surfaces()) {
        i.print(*m_logger);
    }
#endif
}

std::string Settings::sget(const std::string &path) const {
    auto iter = m_proxy.find(path);
    return iter->second;
}

std::string Settings::get(const std::string &path) const {
    auto iter = m_proxy.find(path);
    if (iter == m_proxy.end()) {
#ifndef BENCHMARKING
        m_logger->error(R"(didn't found "{}" in settings)", path);
        print_config();
        throw std::invalid_argument("\"" + path + "\" is an unknown setting");
#endif
    }

#ifndef BENCHMARKING
     m_logger->debug(R"(reading: "{}" with value: "{}")", path, iter->second);
#endif

    return iter->second;
}
}
