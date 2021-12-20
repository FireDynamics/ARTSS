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
#include "../../Functions.h"
#include <fmt/format.h>
#include <iostream>

namespace Settings {
    static const std::string xml_true = "Yes";
    static const char delimiter = ',';
    Map map_parameter_line(const tinyxml2::XMLElement *line, const std::string &context) {
        Map values;
        for (auto e = line->FirstAttribute(); e; e = e->Next()) {
            values[e->Name()] = e->Value();
        }
        for (auto e = line->FirstChildElement(); e; e = e->NextSiblingElement()) {
            if (e->GetText()) {
                values[e->Name()] = e->GetText();
            }
        }
        return values;
    }
    return_xml_data map_parameter_section(const tinyxml2::XMLElement *head, const std::string &context) {
        Map values;
        const tinyxml2::XMLElement *subsection;
        for (auto i = head->FirstChildElement(); i; i = i->NextSiblingElement()) {
            if (i->Name() == context) {
                subsection = i;
                for (auto e = i->FirstAttribute(); e; e = e->Next()) {
                    values[e->Name()] = e->Value();
                }
                for (auto e = i->FirstChildElement(); e; e = e->NextSiblingElement()) {
                    if (e->GetText()) {
                        values[e->Name()] = e->GetText();
                    }
                }
                break;
            }
        }
        return {subsection, values};
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

    domain_parameters parse_domain_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "domain_parameters";
        auto[subsection, values] = map_parameter_section(root, context);
        domain_parameters dp{};

        dp.enable_computational_domain = get_optional_bool(values, "enable_computational_domain", false);
        for (size_t a = 0; a < number_of_axes; a++) {
            auto axis = CoordinateAxis(a);
            std::string axis_name = Mapping::get_axis_name(axis);
            std::string axis_name_low = Utility::to_lower(axis_name);
            dp.start_coords_PD[axis] = get_required_real(values, axis_name + "1", context);
            dp.end_coords_PD[axis] = get_required_real(values, axis_name + "2", context);
            dp.number_of_inner_cells[axis] = get_required_size_t(values, "n" + axis_name_low, context);
            if (dp.enable_computational_domain) {
                dp.start_coords_PD[axis] = get_required_real(values, axis_name_low + "1", context);
                dp.end_coords_PD[axis] = get_required_real(values, axis_name_low + "2", context);
            } else {
                dp.start_coords_PD[axis] = dp.start_coords_PD[axis];
                dp.end_coords_PD[axis] = dp.end_coords_PD[axis];
            }
        }
        return dp;
    }

    visualisation_parameters parse_visualisation_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "visualisation";
        auto[subsection, values] = map_parameter_section(root, context);
        visualisation_parameters vp{};
        vp.save_csv = get_optional_bool(values, "save_csv", false);
        if (vp.save_csv) {
            vp.csv_nth_plot = get_required_size_t(values, "csv_nth_plot", context);
        } else {
            vp.csv_nth_plot = get_optional_size_t(values, "csv_nth_plot", 0);
        }
        vp.save_vtk = get_optional_bool(values, "save_vtk", false);
        if (vp.save_vtk) {
            vp.vtk_nth_plot = get_required_size_t(values, "vtk_nth_plot", context);
        } else {
            vp.vtk_nth_plot = get_optional_size_t(values, "vtk_nth_plot", 0);
        }
        return vp;
    }

    random_parameters parse_random_parameters(const tinyxml2::XMLElement *head,
                                              const std::string &parent_context,
                                              bool is_random) {
        std::string own_context = "random";
        std::string context = create_context(parent_context, own_context);
        auto[subsection, values] = map_parameter_section(head, own_context);
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
            rp.range = get_required_real(values, "range", context);
        }
        return rp;
    }

    namespace initial_conditions {
        uniform parse_uniform_parameters(const Map &values,
                                         const std::string &parent_context) {
            std::string own_context = FunctionNames::uniform;
            std::string context = create_context(parent_context, own_context);
            uniform uniform{};

            uniform.value = get_required_real(values, "val", context);
            return uniform;
        }

        gauss_bubble parse_gauss_bubble_parameters(const Map &values,
                                                   const std::string &parent_context) {
            std::string own_context = FunctionNames::gauss_bubble;
            std::string context = create_context(parent_context, own_context);
            gauss_bubble gauss_bubble{};

            gauss_bubble.l = get_required_real(values, "l", context);
            gauss_bubble.u_lin = get_required_real(values, "u_lin", context);
            gauss_bubble.v_lin = get_required_real(values, "v_lin", context);
            gauss_bubble.w_lin = get_required_real(values, "w_lin", context);
            gauss_bubble.x_shift = get_required_real(values, "x_shift", context);
            gauss_bubble.y_shift = get_required_real(values, "y_shift", context);
            gauss_bubble.z_shift = get_required_real(values, "z_shift", context);
            return gauss_bubble;
        }
    }

    initial_conditions_parameters parse_initial_conditions_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "initial_conditions";
        auto[subsection, values] = map_parameter_section(root, context);
        initial_conditions_parameters icp{};

        icp.random = get_required_bool(values, "random", context);
        icp.random_parameters = parse_random_parameters(subsection, context, icp.random);

        icp.usr_fct = get_required_string(values, "usr_fct", context);
        if (icp.usr_fct == FunctionNames::uniform) {
            icp.ic = initial_conditions::parse_uniform_parameters(values, context);
        } else if (icp.usr_fct == FunctionNames::gauss_bubble) {
            icp.ic = initial_conditions::parse_gauss_bubble_parameters(values, context);
        } else {
            throw config_error(fmt::format("{} has no parsing implementation.", icp.usr_fct));
        }
        return icp;
    }

    boundary parse_boundary(const tinyxml2::XMLElement *head, const std::string &parent_context) {
        std::string own_context = "boundary";
        std::string context = create_context(parent_context, own_context);
        auto values = map_parameter_line(head, own_context);
        boundary boundary{};

        boundary.value = get_required_real(values, "value", context);
        boundary.boundary_condition = Mapping::match_boundary_condition(get_required_string(values, "type", context));
        auto fields = Utility::split(get_required_string(values, "field", context), delimiter);
        for (const std::string &string: fields) {
            boundary.field_type.emplace_back(Mapping::match_field(string));
        }
        auto patches = Utility::split(get_required_string(values, "patch", context), delimiter);
        for (const std::string &string: patches) {
            boundary.patch.emplace_back(Mapping::match_patch(string));
        }
        return boundary;
    }

    boundaries_parameters parse_boundaries_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "boundaries";
        auto[subsection, values] = map_parameter_section(root, context);
        boundaries_parameters bp{};

        for (auto i = subsection->FirstChildElement(); i; i = subsection->NextSiblingElement()) {
            if (i->Name() == std::string("boundary")) {
                bp.boundaries.emplace_back(parse_boundary(i, context));
            }
        }
        return bp;
    }

    obstacle parse_obstacle(const tinyxml2::XMLElement *head, const std::string &parent_context) {
        std::string own_context = "obstacle";
        std::string context = create_context(parent_context, own_context);
        auto[subsection, values] = map_parameter_section(head, own_context);

        obstacle obstacle{};
        obstacle.name = get_required_string(values, "name", context);

        std::string context_geometry = create_context(context, fmt::format("geometry ({})", obstacle.name));
        auto[subsection_geometry, values_geometry] = map_parameter_section(subsection, "geometry");
        for (size_t a = 0; a < number_of_axes; a++) {
            auto axis = CoordinateAxis(a);
            std::string axis_name = Mapping::get_axis_name(axis);
            obstacle.start_coords[axis] = get_required_real(values_geometry, "o" + axis_name + "1", context);
            obstacle.end_coords[axis] = get_required_real(values_geometry, "o" + axis_name + "2", context);
        }

        std::string context_boundaries = create_context(context, fmt::format("boundaries ({})", obstacle.name));
        for (auto i = subsection->FirstChildElement(); i; i = subsection->NextSiblingElement()) {
            if (subsection->Name() == std::string("boundary")) {
                obstacle.boundaries.emplace_back(parse_boundary(i, context_boundaries));
            }
        }
        return obstacle;
    }

    obstacles_parameters parse_obstacles_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "obstacles";
        auto[subsection, values] = map_parameter_section(root, context);
        obstacles_parameters op{};
        op.enabled = get_required_bool(values, "enabled", context);
        if (op.enabled) {
            for (auto i = subsection->FirstChildElement(); i; i = subsection->NextSiblingElement()) {
                op.obstacles.emplace_back(parse_obstacle(i, context));
            }
        }
        op.obstacles.shrink_to_fit();
        return op;
    }

    surfaces_parameters parse_surfaces_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "surfaces";
        auto[subsection, values] = map_parameter_section(root, context);
        surfaces_parameters sp{};

        sp.enabled = get_required_bool(values, "enabled", context);
        if (sp.enabled) {
            throw config_error(
                    fmt::format("handling of surfaces are not fully implemented yet. Please check issue #5."));
        }
        return sp;
    }

    adaption_parameters parse_adaption_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "adaption";
        auto[subsection, values] = map_parameter_section(root, context);
        adaption_parameters ap{};

        ap.enabled = get_required_bool(values, "enabled", context);
        if (ap.enabled) {
            //TODO (issue 178) parse adaption parameters
            throw config_error(
                    fmt::format("handling of adaption is not fully implemented yet. Please check issue #178"));
        }
        return ap;
    }

    logging_parameters parse_logging_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "logging";
        auto[subsection, values] = map_parameter_section(root, context);
        logging_parameters lp{};
        lp.file = get_required_string(values, "file", context);
        lp.level = get_required_string(values, "level", context);
        return lp;
    }

    physical_parameters parse_physical_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "physical_parameters";
        auto[subsection, values] = map_parameter_section(root, context);
        physical_parameters pp{};
        pp.t_end = get_required_real(values, "t_end", context);
        pp.dt = get_required_real(values, "dt", context);
        pp.nu = get_optional_real(values, "nu", 3.1e-5);
        pp.beta = get_optional_real(values, "beta", 3.34e-3);
        pp.g = get_optional_real(values, "g", -9.81);
        pp.kappa = get_optional_real(values, "kappa", 4.25e-5);
        return pp;
    }

    solver_parameters parse_solver_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "solver";
        auto[subsection, values] = map_parameter_section(root, context);
        solver_parameters sp{};

        sp.description = get_required_string(values, "description", context);


        return sp;
    }

    Settings_new parse_settings(const std::string &file_content) {
        tinyxml2::XMLDocument doc;
        doc.Parse(file_content.c_str());
        tinyxml2::XMLElement *root = doc.RootElement();
        return {parse_physical_parameters(root),
                parse_solver_parameters(root),
                parse_domain_parameters(root),
                parse_adaption_parameters(root),
                parse_boundaries_parameters(root),
                parse_obstacles_parameters(root),
                parse_surfaces_parameters(root),
                parse_initial_conditions_parameters(root),
                parse_visualisation_parameters(root),
                parse_logging_parameters(root)};
    }

    Settings_new parse_settings_from_file(const std::filesystem::path &path) {
        std::ifstream in(path);
        std::stringstream sstr;
        sstr << in.rdbuf();
        return parse_settings(sstr.str());
    }


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
        for (auto &i: ordered) {
            m_logger->debug(R"({} = "{}")", i.first, i.second);
        }

        for (const auto &i: get_boundaries()) {
            i.print(*m_logger);
        }

        for (const auto &i: get_obstacles()) {
            i.print(*m_logger);
        }

        for (const auto &i: get_surfaces()) {
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
