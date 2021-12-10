/// \file       Settings.cpp
/// \brief      Access parameters of XML file
/// \date       Nov 17, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.

#include <fstream>
#include <filesystem>
#include "Settings.h"

#include "../Utility.h"
#include <fmt/format.h>

namespace Settings {
using Map = std::map<std::string, std::string>;
static const std::string xml_true = "Yes";

Map map_parameters(tinyxml2::XMLDocument &doc, const std::string &context) {
    Map values;
    for (auto i = doc.RootElement()->FirstChildElement(); i; i = i->NextSiblingElement()) {
        if (i->Name() == context) {
            for (auto e = i->FirstAttribute(); e; e = e->Next()) {
                values[e->Name()] = e->Value();
            }
            for (auto e = i->FirstChildElement(); e; e = e->NextSiblingElement()) {
                values[e->Name()] = e->GetText();
            }
            break;
        }
    }
    return values;
}

size_t get_optional_size_t(const Map &map, const std::string &key, size_t default_value) {
    if (auto iter = map.find(key); iter != map.end()) {
        return std::stoi(iter->second);
    }
    return default_value;
}
size_t get_required_size_t(const Map &map, const std::string &key, const std::string &context) {
    try {
        size_t value = std::stoi(map.at(key));
        return value;
    } catch (const std::exception &e) {
        throw config_error(fmt::format("Value {} is for {} required.", key, context));
    }
}
real get_optional_real(const Map &map, const std::string &key, real default_value) {
    if (auto iter = map.find(key); iter != map.end()) {
        return std::stod(iter->second);
    }
    return default_value;
}
real get_required_real(const Map &map, const std::string &key, const std::string &context) {
    try {
        real value = std::stod(map.at(key));
        return value;
    } catch (const std::exception &e) {
        throw config_error(fmt::format("Value {} is for {} required.", key, context));
    }
}
std::string get_optional_string(const Map &map, const std::string &key, std::string default_value) {
    if (auto iter = map.find(key); iter != map.end()) {
        return iter->second;
    }
    return default_value;
}
std::string get_required_string(const Map &map, const std::string &key, const std::string &context) {
    try {
        return map.at(key);
    } catch (const std::exception &e) {
        throw config_error(fmt::format("Value {} is for {} required.", key, context));
    }
}
bool get_required_bool(const Map &map, const std::string &key, const std::string &context) {
    try {
        return map.at(key) == xml_true;
    } catch (const std::exception &e) {
        throw config_error(fmt::format("Value {} is for {} required.", key, context));
    }
}

domain_parameters parse_domain_parameters(tinyxml2::XMLDocument &doc) {
    std::string context = "domain_parameters";
    Map values = map_parameters(doc, context);
    domain_parameters dp{};
    dp.X1 = get_required_real(values, "X1", context);
    dp.X2 = get_required_real(values, "X2", context);
    dp.Y1 = get_required_real(values, "Y1", context);
    dp.Y2 = get_required_real(values, "Y2", context);
    dp.Z1 = get_required_real(values, "Z1", context);
    dp.Z2 = get_required_real(values, "Z2", context);

    dp.x1 = get_optional_real(values, "x1", 0);
    dp.x2 = get_optional_real(values, "x2", 0);
    dp.y1 = get_optional_real(values, "y1", 0);
    dp.y2 = get_optional_real(values, "y2", 0);
    dp.z1 = get_optional_real(values, "z1", 0);
    dp.z2 = get_optional_real(values, "z2", 0);

    dp.nx = get_required_size_t(values, "nx", context);
    dp.ny = get_required_size_t(values, "ny", context);
    dp.nz = get_required_size_t(values, "nz", context);

    dp.enable_computational_domain = get_required_bool(values, "enable_computational_domain", context);
    return dp;
}
visualisation_parameters parse_visualisation_parameters(tinyxml2::XMLDocument &doc) {
    std::string context = "visualisation";
    Map values = map_parameters(doc, context);
    visualisation_parameters vp{};
    vp.csv_nth_plot = get_optional_size_t(values, "csv_nth_plot", 1);
    vp.vtk_nth_plot = get_optional_size_t(values, "vtk_nth_plot", 1);
    vp.save_csv = get_required_bool(values, "save_csv", context);
    vp.save_vtk = get_required_bool(values, "save_vtk", context);
    return vp;
}
initial_conditions_parameters parse_initial_conditions_parameters(tinyxml2::XMLDocument &doc) {
    std::string context = "initial_conditions";
    Map values = map_parameters(doc, context);
    initial_conditions_parameters icp{};
    icp.random = get_required_bool(values, "random", context);
    icp.usr_fct = get_required_string(values, "usr_fct", context);
    return icp;
}
obstacles_parameters parse_obstacles_parameters(tinyxml2::XMLDocument &doc) {
    std::string context = "obstacles";
    Map values = map_parameters(doc, context);
    obstacles_parameters op{};
    op.enabled = get_required_bool(values, "enabled", context);
    return op;
}
surfaces_parameters parse_surfaces_parameters(tinyxml2::XMLDocument &doc) {
    std::string context = "surfaces";
    Map values = map_parameters(doc, context);
    surfaces_parameters sp{};
    sp.enabled = get_required_bool(values, "enabled", context);
    return sp;
}
logging_parameters parse_logging_parameters(tinyxml2::XMLDocument &doc) {
    std::string context = "logging";
    Map values = map_parameters(doc, context);
    logging_parameters lp{};
    lp.file = get_required_string(values, "file", context);
    lp.level = get_required_string(values, "level", context);
    return lp;
}
physical_parameters parse_physical_parameters(tinyxml2::XMLDocument &doc) {
    std::string context = "physical_parameters";
    Map values = map_parameters(doc, context);
    physical_parameters pp{};
    pp.t_end = get_required_real(values, "t_end", context);
    pp.dt = get_required_real(values, "dt", context);
    pp.nu = get_optional_real(values, "nu", 0.000031);
    pp.beta = get_optional_real(values, "beta", 0.00334);
    pp.g = get_optional_real(values, "g", -9.81);
    pp.kappa = get_optional_real(values, "kappa", 0.0000425);
    return pp;
}
Settings_new parse_settings(const std::string &file_content) {
    tinyxml2::XMLDocument doc;
    doc.Parse(file_content.c_str());
    return {parse_physical_parameters(doc),
            parse_domain_parameters(doc),
            parse_obstacles_parameters(doc),
            parse_surfaces_parameters(doc),
            parse_initial_conditions_parameters(doc),
            parse_visualisation_parameters(doc),
            parse_logging_parameters(doc)};
}

Settings_new parse_settings_from_file(const std::filesystem::path &path) {
    std::ifstream in(path);
    std::stringstream sstr;
    sstr << in.rdbuf();
    return parse_settings(sstr.str());
}


Settings::Settings(std::string path) :
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
    m_logger = Utility::create_logger(*this, "XMLFile");
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
}
