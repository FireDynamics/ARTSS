
#include "Settings.h"
#include "Utility.h"

Settings::Settings(std::string path) :
        filename(path) {
    tinyxml2::XMLDocument doc;
    doc.LoadFile(path.c_str());

    for (auto i = doc.RootElement()->FirstChildElement(); i; i = i->NextSiblingElement()) {
        read_config("", i);
    }

    m_logger = Utility::create_logger(*this, "XMLFile");
    m_logger->debug("start the simulation of \"{}\"", path);
    print_config();
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

std::vector<BoundarySetting> Settings::get_boundaries() const {
    return {};
}

std::vector<ObstacleSetting> Settings::get_obstacles() const {
    return {};
}

std::vector<SurfaceSetting> Settings::get_surfaces() const {
    return {};
}
