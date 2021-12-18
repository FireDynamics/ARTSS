/// \file       Settings.cpp
/// \brief      Access parameters of XML file
/// \date       Nov 17, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.

#include <fstream>
#include <filesystem>
#include "Settings.h"

#include "../Utility.h"
#include "../Mapping.h"
#include <fmt/format.h>

namespace Settings {
using Map = std::map<std::string, std::string>;
static const std::string xml_true = "Yes";
static const std::string xml_false = "No";

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
bool get_optional_bool(const Map &map, const std::string &key, const bool default_value) {
    if (auto iter = map.find(key); iter != map.end()) {
        return iter->second == xml_true;
    }
    return default_value;
}
bool get_required_bool(const Map &map, const std::string &key, const std::string &context) {
    try {
        return map.at(key) == xml_true;
    } catch (const std::exception &e) {
        throw config_error(fmt::format("Value {} is for {} required.", key, context));
    }
}
std::string create_context(const std::string &parent_context, const std::string &own_context) {
    return parent_context + "/" + own_context;
}

domain_parameters parse_domain_parameters(tinyxml2::XMLDocument &doc) {
    std::string context = "domain_parameters";
    Map values = map_parameters(doc, context);
    domain_parameters dp{};

    dp.enable_computational_domain = get_optional_bool(values, "enable_computational_domain", false);
    for (size_t a = 0; a < number_of_axes; a++) {
        auto axis = CoordinateAxis(a);
        std::string axis_name = Mapping::get_axis_name(axis);
        std::string axis_name_low = Utility::to_lower(axis_name);
        dp.start_coords_PD[axis] = get_required_real(values, axis_name + "1", context);
        dp.end_coords_PD[axis] = get_required_real(values, axis_name + "2", context);
        dp.number_of_inner_cells[axis] = get_required_size_t(values, "n" +axis_name_low, context);
        if (dp.enable_computational_domain) {
            dp.start_coords_PD[axis] = get_required_real(values, axis_name_low + "1", context);
            dp.end_coords_PD[axis] = get_required_real(values, axis_name_low + "2", context);
        } else {
           dp.start_coords_PD[axis] = dp.start_coords_PD[axis];
           dp.end_coords_PD[axis] =  dp.end_coords_PD[axis];
        }
    }
    return dp;
}
visualisation_parameters parse_visualisation_parameters(tinyxml2::XMLDocument &doc) {
    std::string context = "visualisation";
    Map values = map_parameters(doc, context);
    visualisation_parameters vp{};
    vp.csv_nth_plot = get_optional_size_t(values, "csv_nth_plot", 1);
    if (vp.csv_nth_plot) {
        vp.save_csv = get_required_bool(values, "save_csv", context);
    } else {
        vp.save_csv = get_optional_bool(values, "csv_nth_plot", 0);
    }
    vp.save_vtk = get_required_bool(values, "save_vtk", context);
    if (vp.vtk_nth_plot) {
        vp.save_vtk = get_required_bool(values, "save_vtk", context);
    } else {
        vp.save_vtk = get_optional_bool(values, "vtk_nth_plot", 0);
    }
    return vp;
}
random_parameters parse_random_parameters(tinyxml2::XMLDocument &doc, const std::string &parent_context, bool is_random) {
    std::string context = create_context(parent_context, "random");
    Map values = map_parameters(doc, context);
    random_parameters rp{};
    if (is_random) {
        rp.absolute = get_optional_bool(values, "absolute", false);
        rp.custom_seed = get_required_bool(values, "custom_seed", context);
        if (rp.custom_seed) {
            rp.seed = get_required_size_t(values, "seed", context);
        } else {
            rp.seed = get_optional_size_t(values, "seed", 0);
        }
        rp.custom_steps = get_required_bool(values, "custom_steps", context);
        if (rp.custom_steps) {
            rp.step_size = get_required_real(values, "step_size", context);
        } else {
            rp.step_size = get_optional_real(values, "step_size", 1);
        }
    }
    return rp;
}
initial_conditions_parameters parse_initial_conditions_parameters(tinyxml2::XMLDocument &doc) {
    std::string context = "initial_conditions";
    Map values = map_parameters(doc, context);
    initial_conditions_parameters icp{};

    icp.random = get_required_bool(values, "random", context);
    icp.random_parameters = parse_random_parameters(doc, context, icp.random);

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
    pp.nu = get_optional_real(values, "nu", 3.1e-5);
    pp.beta = get_optional_real(values, "beta", 3.34e-3);
    pp.g = get_optional_real(values, "g", -9.81);
    pp.kappa = get_optional_real(values, "kappa", 4.25e-5);
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
