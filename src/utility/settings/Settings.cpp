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
#include "../../interfaces/ISolver.h"
#include "../../solver/SolverSelection.h"
#include <fmt/format.h>
#include <iostream>

namespace Settings {
    static const std::string xml_true = "Yes";
    static const char delimiter = ',';

    Map map_parameter_line(const tinyxml2::XMLElement *line) {
        Map values;
        for (const auto *e = line->FirstAttribute(); e; e = e->Next()) {
            values[e->Name()] = e->Value();
        }
        for (const auto *e = line->FirstChildElement(); e; e = e->NextSiblingElement()) {
            if (e->GetText()) {
                values[e->Name()] = e->GetText();
            }
        }
        return values;
    }

    return_xml_data map_parameter_section(const tinyxml2::XMLElement *head, const std::string &context) {
        Map values;
        const tinyxml2::XMLElement *subsection;
        for (const auto *i = head->FirstChildElement(); i; i = i->NextSiblingElement()) {
            if (i->Name() == context) {
                subsection = i;
                for (const auto *e = i->FirstAttribute(); e; e = e->Next()) {
                    values[e->Name()] = e->Value();
                }
                for (const auto *e = i->FirstChildElement(); e; e = e->NextSiblingElement()) {
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
        const auto iter = map.find(key);
        if (iter != map.end()) {
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
        const auto iter = map.find(key);
        if (iter != map.end()) {
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
        const auto iter = map.find(key);
        if (iter != map.end()) {
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
        const auto iter = map.find(key);
        if (iter != map.end()) {
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
                dp.start_coords_CD[axis] = get_required_real(values, axis_name_low + "1", context);
                dp.end_coords_CD[axis] = get_required_real(values, axis_name_low + "2", context);
            } else {
                dp.start_coords_CD[axis] = dp.start_coords_PD[axis];
                dp.end_coords_CD[axis] = dp.end_coords_PD[axis];
            }
        }
        return dp;
    }

    visualisation_parameters parse_visualisation_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "visualisation";
        auto[subsection, values] = map_parameter_section(root, context);
        visualisation_parameters vp{};
        vp.save_csv = get_required_bool(values, "save_csv", context);
        if (vp.save_csv) {
            vp.csv_nth_plot = get_required_size_t(values, "csv_nth_plot", context);
        }
        vp.save_vtk = get_required_bool(values, "save_vtk", context);
        if (vp.save_vtk) {
            vp.vtk_nth_plot = get_required_size_t(values, "vtk_nth_plot", context);
        }
        return vp;
    }

    random_parameters parse_random_parameters(const tinyxml2::XMLElement *head,
                                              const std::string &parent_context) {
        std::string own_context = "random";
        std::string context = create_context(parent_context, own_context);
        auto[subsection, values] = map_parameter_section(head, own_context);
        random_parameters rp{};
        rp.absolute = get_required_bool(values, "absolute", context);
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

        exp_sinus_prod parse_exp_sinus_prod(const Map &values,
                                            const std::string &parent_context) {
            std::string own_context = FunctionNames::exp_sinus_prod;
            std::string context = create_context(parent_context, own_context);
            exp_sinus_prod exp_sinus_prod{};

            exp_sinus_prod.l = get_required_real(values, "l", context);
            return exp_sinus_prod;
        }

        gauss_bubble parse_gauss_bubble_parameters(const Map &values,
                                                   const std::string &parent_context) {
            std::string own_context = FunctionNames::gauss_bubble;
            std::string context = create_context(parent_context, own_context);
            gauss_bubble gauss_bubble{};

            gauss_bubble.l = get_required_real(values, "l", context);
            gauss_bubble.velocity_lin[CoordinateAxis::X] = get_required_real(values, "u_lin", context);
            gauss_bubble.velocity_lin[CoordinateAxis::Y] = get_required_real(values, "v_lin", context);
            gauss_bubble.velocity_lin[CoordinateAxis::Z] = get_required_real(values, "w_lin", context);
            gauss_bubble.shift[CoordinateAxis::X] = get_required_real(values, "x_shift", context);
            gauss_bubble.shift[CoordinateAxis::Y] = get_required_real(values, "y_shift", context);
            gauss_bubble.shift[CoordinateAxis::Z] = get_required_real(values, "z_shift", context);
            return gauss_bubble;
        }

        hat parse_hat_parameters(const Map &values, const std::string &parent_context) {
            std::string own_context = FunctionNames::hat;
            std::string context = create_context(parent_context, own_context);
            hat hat{};

            hat.start_coords[CoordinateAxis::X] = get_required_real(values, "x1", context);
            hat.start_coords[CoordinateAxis::Y] = get_required_real(values, "y1", context);
            hat.start_coords[CoordinateAxis::Z] = get_required_real(values, "z1", context);
            hat.end_coords[CoordinateAxis::X] = get_required_real(values, "x2", context);
            hat.end_coords[CoordinateAxis::Y] = get_required_real(values, "y2", context);
            hat.end_coords[CoordinateAxis::Z] = get_required_real(values, "z2", context);
            hat.val_out = get_required_real(values, "val_out", context);
            hat.val_in = get_required_real(values, "val_in", context);
            return hat;
        }

        drift parse_drift_parameters(const Map &values, const std::string &parent_context) {
            std::string own_context = FunctionNames::drift;
            std::string context = create_context(parent_context, own_context);
            drift drift{};

            drift.velocity_lin[CoordinateAxis::X] = get_required_real(values, "u_lin", context);
            drift.velocity_lin[CoordinateAxis::Y] = get_required_real(values, "v_lin", context);
            drift.velocity_lin[CoordinateAxis::Z] = get_required_real(values, "w_lin", context);
            drift.pa = get_required_real(values, "pa", context);
            return drift;
        }

        vortex parse_vortex_parameters(const Map &values, const std::string &parent_context) {
            std::string own_context = FunctionNames::vortex;
            std::string context = create_context(parent_context, own_context);
            vortex vortex{};

            vortex.velocity_lin[CoordinateAxis::X] = get_required_real(values, "u_lin", context);
            vortex.velocity_lin[CoordinateAxis::Y] = get_required_real(values, "v_lin", context);
            vortex.velocity_lin[CoordinateAxis::Z] = get_required_real(values, "w_lin", context);
            vortex.pa = get_required_real(values, "pa", context);
            vortex.rhoa = get_required_real(values, "rhoa", context);
            return vortex;
        }

        mc_dermott parse_mc_dermott_parameters(const Map &values, const std::string &parent_context) {
            std::string own_context = FunctionNames::mcdermott;
            std::string context = create_context(parent_context, own_context);
            mc_dermott mc_dermott{};

            mc_dermott.A = get_required_real(values, "A", context);
            return mc_dermott;
        }

        jet parse_jet_parameters(const Map &values, const std::string &parent_context) {
            std::string own_context = FunctionNames::jet;
            std::string context = create_context(parent_context, own_context);
            jet jet{};

            jet.dir = Mapping::match_axis(get_required_string(values, "dir", context));
            if (jet.dir == CoordinateAxis::X) {
                jet.start_coords[CoordinateAxis::Y] = get_required_real(values, "y1", context);
                jet.start_coords[CoordinateAxis::Z] = get_required_real(values, "z1", context);
                jet.end_coords[CoordinateAxis::Y] = get_required_real(values, "y2", context);
                jet.end_coords[CoordinateAxis::Z] = get_required_real(values, "z2", context);
            } else if (jet.dir == CoordinateAxis::Y) {
                jet.start_coords[CoordinateAxis::X] = get_required_real(values, "x1", context);
                jet.start_coords[CoordinateAxis::Z] = get_required_real(values, "z1", context);
                jet.end_coords[CoordinateAxis::X] = get_required_real(values, "x2", context);
                jet.end_coords[CoordinateAxis::Z] = get_required_real(values, "z2", context);
            } else {
                jet.start_coords[CoordinateAxis::X] = get_required_real(values, "x1", context);
                jet.start_coords[CoordinateAxis::Y] = get_required_real(values, "y1", context);
                jet.end_coords[CoordinateAxis::X] = get_required_real(values, "x2", context);
                jet.end_coords[CoordinateAxis::Y] = get_required_real(values, "y2", context);
            }
            jet.value = get_required_real(values, "value", context);
            return jet;
        }
        layers_temperature parse_layers_parameters(const Map &values, const std::string &parent_context) {
            std::string own_context = FunctionNames::layers;
            std::string context = create_context(parent_context, own_context);
            layers_temperature layers{};

            layers.number_of_layers = get_required_size_t(values, "n_layers", context);
            if (layers.number_of_layers == 0) {
                throw config_error(fmt::format("Numbers of layers has to be at least 1. Current value: {}",
                                               layers.number_of_layers));
            }
            layers.dir = Mapping::match_axis(get_required_string(values, "dir", context));
            layers.values.reserve(layers.number_of_layers);
            layers.borders.reserve(layers.number_of_layers - 1);
            for (size_t i = 1; i < layers.number_of_layers; i++) {
                layers.values.emplace_back(get_required_real(values, "value_" + std::to_string(i), context));
                layers.borders.emplace_back(get_required_real(values, "border_" + std::to_string(i), context));
            }
            layers.values.emplace_back(
                    get_required_real(values, "value_" + std::to_string(layers.number_of_layers), context));
            return layers;
        }

        sin_sin_sin parse_sin_sin_sin_parameters(const Map &values, const std::string &parent_context) {
            std::string own_context = FunctionNames::sin_sin_sin;
            std::string context = create_context(parent_context, own_context);
            sin_sin_sin sin{};
            sin.l = get_required_real(values, "l", context);
            return sin;
        }
    }

    initial_conditions_parameters parse_initial_conditions_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "initial_conditions";
        auto[subsection, values] = map_parameter_section(root, context);
        initial_conditions_parameters icp{};

        icp.random = get_required_bool(values, "random", context);
        if (icp.random) {
            icp.random_parameters = parse_random_parameters(subsection, context);
        }

        icp.usr_fct = get_required_string(values, "usr_fct", context);
        if (icp.usr_fct == FunctionNames::uniform) {
            icp.ic = initial_conditions::parse_uniform_parameters(values, context);
        } else if (icp.usr_fct == FunctionNames::drift) {
            icp.ic = initial_conditions::parse_drift_parameters(values, context);
        } else if (icp.usr_fct == FunctionNames::gauss_bubble) {
            icp.ic = initial_conditions::parse_gauss_bubble_parameters(values, context);
        } else if (icp.usr_fct == FunctionNames::hat) {
            icp.ic = initial_conditions::parse_hat_parameters(values, context);
        } else if (icp.usr_fct == FunctionNames::exp_sinus_prod) {
            icp.ic = initial_conditions::parse_exp_sinus_prod(values, context);
        } else if (icp.usr_fct == FunctionNames::exp_sinus_sum) {
            // no values
        } else if (icp.usr_fct == FunctionNames::jet) {
            icp.ic = initial_conditions::parse_jet_parameters(values, context);
        } else if (icp.usr_fct == FunctionNames::layers) {
            icp.ic = initial_conditions::parse_layers_parameters(values, context);
        } else if (icp.usr_fct == FunctionNames::mcdermott) {
            icp.ic = initial_conditions::parse_mc_dermott_parameters(values, context);
        } else if (icp.usr_fct == FunctionNames::sin_sin_sin) {
            icp.ic = initial_conditions::parse_sin_sin_sin_parameters(values, context);
        } else if (icp.usr_fct == FunctionNames::vortex) {
            icp.ic = initial_conditions::parse_vortex_parameters(values, context);
        } else if (icp.usr_fct == FunctionNames::zero) {
            // do nothing
        } else {
            throw config_error(fmt::format("{} has no parsing implementation.", icp.usr_fct));
        }
        return icp;
    }

    boundary parse_boundary(const tinyxml2::XMLElement *head, const std::string &parent_context) {
        std::string own_context = "boundary";
        std::string context = create_context(parent_context, own_context);
        auto values = map_parameter_line(head);
        boundary boundary{};

        boundary.boundary_condition = Mapping::match_boundary_condition(get_required_string(values, "type", context));
        if (boundary.boundary_condition != BoundaryCondition::PERIODIC) {
            boundary.value = get_required_real(values, "value", context);
        }
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

    boundary_parameters parse_boundaries_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "boundaries";
        auto[subsection, values] = map_parameter_section(root, context);
        boundary_parameters bp{};

        for (auto i = subsection->FirstChildElement(); i; i = i->NextSiblingElement()) {
            if (i->Name() == std::string("boundary")) {
                bp.boundaries.emplace_back(parse_boundary(i, context));
            }
        }
        return bp;
    }

    obstacle parse_obstacle(const tinyxml2::XMLElement *head, const std::string &parent_context) {
        std::string own_context = "obstacle";
        std::string context = create_context(parent_context, own_context);
        auto values = map_parameter_line(head);

        obstacle obstacle{};
        obstacle.name = get_required_string(values, "name", context);

        std::string context_geometry = create_context(context, fmt::format("geometry ({})", obstacle.name));
        auto[head_geometry, values_geometry] = map_parameter_section(head, "geometry");
        for (size_t a = 0; a < number_of_axes; a++) {
            auto axis = CoordinateAxis(a);
            std::string axis_name = Utility::to_lower(Mapping::get_axis_name(axis));
            obstacle.start_coords[axis] = get_required_real(values_geometry, "o" + axis_name + "1", context);
            obstacle.end_coords[axis] = get_required_real(values_geometry, "o" + axis_name + "2", context);
        }

        std::string context_boundaries = create_context(context, fmt::format("boundaries ({})", obstacle.name));
        for (const auto *i = head->FirstChildElement(); i; i = i->NextSiblingElement()) {
            if (i->Name() == std::string("boundary")) {
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
            for (const auto *i = subsection->FirstChildElement(); i; i = i->NextSiblingElement()) {
                op.obstacles.emplace_back(parse_obstacle(i, context));
            }
        }
        op.obstacles.shrink_to_fit();
        return op;
    }

    boundary parse_surface_boundary(const tinyxml2::XMLElement *head, const std::string &parent_context) {
        std::string own_context = "boundary";
        std::string context = create_context(parent_context, own_context);
        auto values = map_parameter_line(head);
        boundary boundary{};

        boundary.boundary_condition = Mapping::match_boundary_condition(get_required_string(values, "type", context));
        if (boundary.boundary_condition != BoundaryCondition::PERIODIC) {
            boundary.value = get_required_real(values, "value", context);
        }
        auto fields = Utility::split(get_required_string(values, "field", context), delimiter);
        for (const std::string &string: fields) {
            boundary.field_type.emplace_back(Mapping::match_field(string));
        }
        return boundary;
    }

    surface parse_surface(const tinyxml2::XMLElement *head, const std::string &parent_context) {
        std::string own_context = "surface";
        std::string context = create_context(parent_context, own_context);
        auto values = map_parameter_line(head);

        surface surface{};
        surface.name = get_required_string(values, "name", context);
        surface.patch = Mapping::match_patch(get_required_string(values, "patch", context));

        std::string context_geometry = create_context(context, fmt::format("geometry ({})", surface.name));
        auto[head_geometry, values_geometry] = map_parameter_section(head, "geometry");
        std::vector<CoordinateAxis> axes = {CoordinateAxis::X, CoordinateAxis::Y, CoordinateAxis::Z};
        CoordinateAxis coord_axis = Mapping::to_axis(surface.patch);
        auto tmp = std::remove(axes.begin(), axes.end(), coord_axis);
        for (auto it = axes.begin(); it != tmp; ++it) {
            auto axis = *it;
            std::string axis_name = Utility::to_lower(Mapping::get_axis_name(axis));
            surface.start_coords[axis] = get_required_real(values_geometry, "s" + axis_name + "1", context);
            surface.end_coords[axis] = get_required_real(values_geometry, "s" + axis_name + "2", context);
        }

        std::string context_boundaries = create_context(context, fmt::format("boundaries ({})", surface.name));
        for (const auto *i = head->FirstChildElement(); i; i = i->NextSiblingElement()) {
            if (i->Name() == std::string("boundary")) {
                surface.boundaries.emplace_back(parse_surface_boundary(i, context_boundaries));
                surface.boundaries.back().patch.emplace_back(surface.patch);
            }
        }
        return surface;
    }
    surfaces_parameters parse_surfaces_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "surfaces";
        auto[subsection, values] = map_parameter_section(root, context);
        surfaces_parameters sp{};

        sp.enabled = get_required_bool(values, "enabled", context);
        if (sp.enabled) {
            for (const auto *i = subsection->FirstChildElement(); i; i = i->NextSiblingElement()) {
                sp.surfaces.emplace_back(parse_surface(i, context));
            }
        }
        sp.surfaces.shrink_to_fit();
        return sp;
    }

    adaption_parameters parse_adaption_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "adaption";
        auto[subsection, values] = map_parameter_section(root, context);
        adaption_parameters ap{};

        ap.enabled = get_required_bool(values, "dynamic", context);
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

    physical_parameters
    parse_physical_parameters(const tinyxml2::XMLElement *root, const std::string &solver_description) {
        std::string context = "physical_parameters";
        auto[subsection, values] = map_parameter_section(root, context);
        physical_parameters pp{};
        pp.t_end = get_required_real(values, "t_end", context);
        pp.dt = get_required_real(values, "dt", context);
        if (solver_description != SolverTypes::AdvectionSolver &&
            solver_description != SolverTypes::PressureSolver) {
            pp.nu = get_required_real(values, "nu", context);
        }
        if (solver_description == SolverTypes::NSSolver ||
            solver_description == SolverTypes::NSTempSolver ||
            solver_description == SolverTypes::NSTurbSolver ||
            solver_description == SolverTypes::NSTempTurbSolver ||
            solver_description == SolverTypes::NSTempTurbConSolver ||
            solver_description == SolverTypes::NSTempConSolver) {
            // navier stokes
            pp.rhoa = get_optional_real(values, "rhoa", 1);
        }
        pp.beta = get_optional_real(values, "beta", 3.34e-3);  // for buoyancy
        pp.g = get_optional_real(values, "g", -9.81);  // for buoyancy
        if (solver_description == SolverTypes::NSTempSolver ||
            solver_description == SolverTypes::NSTempTurbSolver ||
            solver_description == SolverTypes::NSTempTurbConSolver ||
            solver_description == SolverTypes::NSTempConSolver) {
            // temperature
            pp.kappa = get_required_real(values, "kappa", context);
        }
        if (solver_description == SolverTypes::NSTempTurbConSolver ||
            solver_description == SolverTypes::NSTempConSolver) {
            // concentration
            pp.gamma = get_required_real(values, "gamma", context);
        }
        return pp;
    }

    namespace solver {
        advection_solver parse_advection_solver(const tinyxml2::XMLElement *head, const std::string &parent_context) {
            std::string own_context = "advection";
            std::string context = create_context(parent_context, own_context);
            auto[subsection, values] = map_parameter_section(head, own_context);

            advection_solver solver_advection{};
            solver_advection.type = get_required_string(values, "type", context);
            if (solver_advection.type == AdvectionMethods::SemiLagrangian) {
                // no values
            } else {
                throw config_error(
                        fmt::format("advection solver '{}' has no parsing implementation.", solver_advection.type));
            }

            auto fields = Utility::split(get_required_string(values, "field", context), delimiter);
            for (const std::string &string: fields) {
                solver_advection.fields.emplace_back(Mapping::match_field(string));
            }
            return solver_advection;
        }

        namespace diffusion_solvers {
            colored_gauss_seidel parse_colored_gauss_seidel(const Map &values, const std::string &parent_context) {
                std::string own_context = DiffusionMethods::ColoredGaussSeidel;
                std::string context = create_context(parent_context, own_context);

                colored_gauss_seidel colored_gauss_seidel{};
                colored_gauss_seidel.w = get_required_real(values, "w", context);
                colored_gauss_seidel.tol_res = get_required_real(values, "tol_res", context);
                colored_gauss_seidel.max_iter = get_required_size_t(values, "max_iter", context);
                return colored_gauss_seidel;
            }

            jacobi parse_jacobi(const Map &values, const std::string &parent_context) {
                std::string own_context = DiffusionMethods::Jacobi;
                std::string context = create_context(parent_context, own_context);

                jacobi jacobi{};
                jacobi.w = get_required_real(values, "w", context);
                jacobi.tol_res = get_required_real(values, "tol_res", context);
                jacobi.max_iter = get_required_size_t(values, "max_iter", context);
                return jacobi;
            }
        }

        diffusion_solver parse_diffusion_solver(const tinyxml2::XMLElement *head, const std::string &parent_context) {
            std::string own_context = "diffusion";
            std::string context = create_context(parent_context, own_context);
            auto[subsection, values] = map_parameter_section(head, own_context);

            diffusion_solver solver_diffusion{};
            solver_diffusion.type = get_required_string(values, "type", context);
            if (solver_diffusion.type == DiffusionMethods::Jacobi) {
                solver_diffusion.solver = diffusion_solvers::parse_jacobi(values, context);
            } else if (solver_diffusion.type == DiffusionMethods::ColoredGaussSeidel) {
                solver_diffusion.solver = diffusion_solvers::parse_colored_gauss_seidel(values, context);
            } else if (solver_diffusion.type == DiffusionMethods::Explicit) {
                // no values
            } else {
                throw config_error(
                        fmt::format("diffusion solver '{}' has no parsing implementation.", solver_diffusion.type));
            }

            auto fields = Utility::split(get_required_string(values, "field", context), delimiter);
            for (const std::string &string: fields) {
                solver_diffusion.fields.emplace_back(Mapping::match_field(string));
            }
            return solver_diffusion;
        }

        namespace turbulence_solvers {
            const_smagorinsky parse_const_smagorinsky(const Map &values, const std::string &parent_context) {
                std::string own_context = "ConstSmagorinsky";
                std::string context = create_context(parent_context, own_context);

                const_smagorinsky cs{};
                cs.cs = get_required_real(values, "Cs", context);
                return cs;
            }
        }

        turbulence_solver parse_turbulence_solver(const tinyxml2::XMLElement *head, const std::string &parent_context) {
            std::string own_context = "turbulence";
            std::string context = create_context(parent_context, own_context);
            auto[subsection, values] = map_parameter_section(head, own_context);

            turbulence_solver solver_turbulence{};
            solver_turbulence.type = get_required_string(values, "type", context);
            if (solver_turbulence.type == TurbulenceMethods::ConstSmagorinsky) {
                solver_turbulence.solver = turbulence_solvers::parse_const_smagorinsky(values, context);
            } else if (solver_turbulence.type == TurbulenceMethods::DynamicSmagorinsky) {
                // no values
            } else {
                throw config_error(
                        fmt::format("turbulence solver '{}' has no parsing implementation.", solver_turbulence.type));
            }

            return solver_turbulence;
        }

        namespace source_solvers {
            buoyancy parse_buoyancy(const Map &values, const std::string &parent_context) {
                std::string own_context = SourceMethods::Buoyancy;
                std::string context = create_context(parent_context, own_context);

                buoyancy buoyancy{};
                buoyancy.use_init_values = get_required_bool(values, "use_init_values", context);
                if (!buoyancy.use_init_values) {
                    buoyancy.ambient_temperature_value = get_required_real(values, "ambient_temperature_value",
                                                                           context);
                }
                return buoyancy;
            }

            uniform parse_uniform(const Map &values, const std::string &parent_context) {
                std::string own_context = SourceMethods::Uniform;
                std::string context = create_context(parent_context, own_context);

                uniform uniform{};
                uniform.velocity_value[CoordinateAxis::X] = get_required_real(values, "val_x", context);
                uniform.velocity_value[CoordinateAxis::Y] = get_required_real(values, "val_y", context);
                uniform.velocity_value[CoordinateAxis::Z] = get_required_real(values, "val_z", context);
                return uniform;
            }
        }

        source_solver parse_source_solver(const tinyxml2::XMLElement *head, const std::string &parent_context) {
            std::string own_context = "source";
            std::string context = create_context(parent_context, own_context);
            auto[subsection, values] = map_parameter_section(head, own_context);

            source_solver solver_source{};
            solver_source.type = get_required_string(values, "type", context);
            solver_source.force_fct = get_required_string(values, "force_fct", context);
            if (solver_source.force_fct == SourceMethods::Buoyancy) {
                solver_source.force_function = source_solvers::parse_buoyancy(values, context);
            } else if (solver_source.force_fct == SourceMethods::Uniform) {
                solver_source.force_function = source_solvers::parse_uniform(values, context);
            } else if (solver_source.force_fct == SourceMethods::Zero) {
                // do nothing
            } else {
                throw config_error(fmt::format("source type '{}' has no parsing implementation.", solver_source.type));
            }
            auto fields = Utility::split(get_required_string(values, "dir", context), delimiter);
            for (const std::string &string: fields) {
                solver_source.direction.emplace_back(Mapping::match_axis(string));
            }
            return solver_source;
        }

        namespace pressure_solvers {
            vcycle_mg
            parse_vcycle_mg(const Map &values, const tinyxml2::XMLElement *head, const std::string &parent_context) {
                std::string own_context = "vcycle";
                std::string context = create_context(parent_context, own_context);
                vcycle_mg mg{};
                mg.n_level = get_required_size_t(values, "n_level", context);
                mg.n_cycle = get_required_size_t(values, "n_cycle", context);
                mg.max_cycle = get_required_size_t(values, "max_cycle", context);
                mg.n_relax = get_required_size_t(values, "n_relax", context);
                mg.tol_res = get_required_real(values, "tol_res", context);
                mg.smoother = parse_diffusion_solver(head, context);
                return mg;
            }
        }

        pressure_solver parse_pressure_solver(const tinyxml2::XMLElement *head, const std::string &parent_context) {
            std::string own_context = "pressure";
            std::string context = create_context(parent_context, own_context);
            auto[subsection, values] = map_parameter_section(head, own_context);

            pressure_solver solver_pressure{};
            solver_pressure.type = get_required_string(values, "type", context);
            solver_pressure.field = Mapping::match_field(get_required_string(values, "field", context));
            if (solver_pressure.type == PressureMethods::VCycleMG) {
                solver_pressure.solver = pressure_solvers::parse_vcycle_mg(values, subsection, context);
            } else {
                throw config_error(
                        fmt::format("pressure solver '{}' has no parsing implementation.", solver_pressure.type));
            }
            return solver_pressure;
        }

        namespace sources {
            gauss parse_gauss(const Map &values, const std::string &parent_context) {
                std::string context = create_context(parent_context, SourceMethods::Gauss);
                gauss gauss{};
                gauss.heat_release_rate = get_required_real(values, "HRR", context);
                gauss.heat_capacity = get_required_real(values, "cp", context);
                for (const auto axis: {CoordinateAxis::X, CoordinateAxis::Y, CoordinateAxis::Z}) {
                    std::string axis_name = Utility::to_lower(Mapping::get_axis_name(axis));
                    gauss.position[axis] = get_required_real(values, axis_name + "0", context);
                    gauss.dimension[axis] = get_required_real(values, "sigma_" + axis_name, context);
                }
                gauss.tau = get_required_real(values, "tau", context);
                return gauss;
            }
            cube parse_cube(const Map &values, const std::string &parent_context) {
                std::string context = create_context(parent_context, SourceMethods::Cube);
                cube cube{};
                for (const auto axis: {CoordinateAxis::X, CoordinateAxis::Y, CoordinateAxis::Z}) {
                    std::string axis_name = Utility::to_lower(Mapping::get_axis_name(axis));
                    cube.coords_start[axis] = get_required_real(values, axis_name + "_start", context);
                    cube.coords_end[axis] = get_required_real(values, axis_name + "_end", context);
                }
                cube.value = get_required_real(values, "val", context);
                return cube;
            }
        }

        temperature_solver parse_temperature_solver(const tinyxml2::XMLElement *head, const std::string &parent_context,
                                                    bool has_turbulence) {
            std::string own_context = "temperature";
            std::string context = create_context(parent_context, own_context);
            auto[subsection, values] = map_parameter_section(head, own_context);

            temperature_solver solver_temperature{};
            solver_temperature.advection = parse_advection_solver(subsection, context);
            solver_temperature.diffusion = parse_diffusion_solver(subsection, context);

            if (has_turbulence) {
                std::string context_turb = "turbulence";
                auto[subsection_turb, values_turb] = map_parameter_section(subsection, context_turb);
                solver_temperature.has_turbulence = get_required_bool(values_turb, "include",
                                                                      create_context(context, context_turb));
                if (solver_temperature.has_turbulence) {
                    solver_temperature.prandtl_number = get_required_real(values_turb, "Pr_T",
                                                                          create_context(context, context_turb));
                }
            } else {
                solver_temperature.has_turbulence = false;
            }

            std::string context_source = "source";
            auto[subsection_source, values_source] = map_parameter_section(subsection, context_source);
            solver_temperature.source.type = get_required_string(values_source, "type", context_source);
            auto fields = Utility::split(get_required_string(values_source, "dir", context_source), delimiter);
            for (const std::string &string: fields) {
                solver_temperature.source.dir.emplace_back(Mapping::match_axis(string));
            }
            solver_temperature.source.dissipation = get_required_bool(values_source, "dissipation", context_source);
            solver_temperature.source.random = get_required_bool(values_source, "random", context_source);
            if (solver_temperature.source.random) {
                solver_temperature.source.random_parameters = parse_random_parameters(subsection_source,
                                                                                      create_context(context,
                                                                                                     context_source));
            }
            solver_temperature.source.temp_fct = get_required_string(values_source, "temp_fct", context_source);
            if (solver_temperature.source.temp_fct == SourceMethods::Gauss) {
                solver_temperature.source.temp_function = sources::parse_gauss(values_source, context_source);
            } else if (solver_temperature.source.temp_fct == SourceMethods::Cube) {
                solver_temperature.source.temp_function = sources::parse_cube(values_source, context_source);
            } else if (solver_temperature.source.temp_fct == SourceMethods::Buoyancy || solver_temperature.source.temp_fct == SourceMethods::Zero) {
                // do nothing
            } else {
                throw config_error(fmt::format("temperature source function '{}' has no parsing implementation.",
                                               solver_temperature.source.temp_fct));
            }
            return solver_temperature;
        }

        concentration_solver
        parse_concentration_solver(const tinyxml2::XMLElement *head, const std::string &parent_context,
                                   bool has_turbulence) {
            std::string own_context = "concentration";
            std::string context = create_context(parent_context, own_context);
            auto[subsection, values] = map_parameter_section(head, own_context);

            concentration_solver solver_concentration{};
            solver_concentration.advection = parse_advection_solver(subsection, context);
            solver_concentration.diffusion = parse_diffusion_solver(subsection, context);

            if (has_turbulence) {
                std::string context_turb = "turbulence";
                auto[subsection_turb, values_turb] = map_parameter_section(subsection, context_turb);
                solver_concentration.has_turbulence =
                        get_required_bool(values_turb, "include",
                                          create_context(context, context_turb));
                if (solver_concentration.has_turbulence) {
                    solver_concentration.turbulent_schmidt_number =
                            get_required_real(values_turb, "Sc_T",
                                              create_context(context, context_turb));
                }
            } else {
                solver_concentration.has_turbulence = false;
            }

    std::string context_source = "source";
            auto[subsection_source, values_source] = map_parameter_section(subsection, context_source);
            solver_concentration.source.type = get_required_string(values_source, "type", context_source);
            auto fields = Utility::split(get_required_string(values_source, "dir", context_source), delimiter);
            for (const std::string &string: fields) {
                solver_concentration.source.dir.emplace_back(Mapping::match_axis(string));
            }
            solver_concentration.source.random = get_required_bool(values_source, "random", context_source);
            if (solver_concentration.source.random) {
                solver_concentration.source.random_parameters = parse_random_parameters(subsection_source,
                                                                                        create_context(context,
                                                                                                       context_source));
            }
            solver_concentration.source.con_fct = get_required_string(values_source, "con_fct", context_source);
            if (solver_concentration.source.con_fct == SourceMethods::Gauss) {
                solver_concentration.source.con_function = sources::parse_gauss(values_source, context_source);
            } else if (solver_concentration.source.con_fct == SourceMethods::Cube) {
                solver_concentration.source.con_function = sources::parse_cube(values_source, context_source);
            } else {
                throw config_error(fmt::format("concentration source function '{}' has no parsing implementation.",
                                               solver_concentration.source.con_fct));
            }
            return solver_concentration;
        }

        solution parse_solution(const tinyxml2::XMLElement *head, const std::string &parent_context) {
            std::string own_context = "solution";
            std::string context = create_context(parent_context, own_context);
            auto[subsection, values] = map_parameter_section(head, own_context);

            solution sol{};
            sol.analytical_solution = get_required_bool(values, "available", context);
            if (sol.analytical_solution) {
                sol.solution_tolerance = get_optional_real(values, "tol", 1e-3);
            }
            return sol;
        }
    }

    solver_parameters parse_solver_parameters(const tinyxml2::XMLElement *root) {
        std::string context = "solver";
        auto[subsection, values] = map_parameter_section(root, context);
        solver_parameters sp{};

        sp.description = get_required_string(values, "description", context);

        if (sp.description == SolverTypes::AdvectionSolver) {
            sp.advection = solver::parse_advection_solver(subsection, context);
        } else if (sp.description == SolverTypes::AdvectionDiffusionSolver) {
            sp.advection = solver::parse_advection_solver(subsection, context);
            sp.diffusion = solver::parse_diffusion_solver(subsection, context);
        } else if (sp.description == SolverTypes::DiffusionSolver) {
            sp.diffusion = solver::parse_diffusion_solver(subsection, context);
        } else if (sp.description == SolverTypes::DiffusionTurbSolver) {
            sp.diffusion = solver::parse_diffusion_solver(subsection, context);
            sp.turbulence = solver::parse_turbulence_solver(subsection, context);
        } else if (sp.description == SolverTypes::NSSolver) {
            sp.advection = solver::parse_advection_solver(subsection, context);
            sp.diffusion = solver::parse_diffusion_solver(subsection, context);
            sp.pressure = solver::parse_pressure_solver(subsection, context);
            sp.source = solver::parse_source_solver(subsection, context);
        } else if (sp.description == SolverTypes::NSTempSolver) {
            sp.advection = solver::parse_advection_solver(subsection, context);
            sp.diffusion = solver::parse_diffusion_solver(subsection, context);
            sp.pressure = solver::parse_pressure_solver(subsection, context);
            sp.temperature = solver::parse_temperature_solver(subsection, context, false);
            sp.source = solver::parse_source_solver(subsection, context);
        } else if (sp.description == SolverTypes::NSTempConSolver) {
            sp.advection = solver::parse_advection_solver(subsection, context);
            sp.diffusion = solver::parse_diffusion_solver(subsection, context);
            sp.pressure = solver::parse_pressure_solver(subsection, context);
            sp.temperature = solver::parse_temperature_solver(subsection, context, false);
            sp.concentration = solver::parse_concentration_solver(subsection, context, false);
            sp.source = solver::parse_source_solver(subsection, context);
        } else if (sp.description == SolverTypes::NSTempTurbSolver) {
            sp.advection = solver::parse_advection_solver(subsection, context);
            sp.diffusion = solver::parse_diffusion_solver(subsection, context);
            sp.pressure = solver::parse_pressure_solver(subsection, context);
            sp.temperature = solver::parse_temperature_solver(subsection, context, true);
            sp.turbulence = solver::parse_turbulence_solver(subsection, context);
            sp.source = solver::parse_source_solver(subsection, context);
        } else if (sp.description == SolverTypes::NSTempTurbConSolver) {
            sp.advection = solver::parse_advection_solver(subsection, context);
            sp.diffusion = solver::parse_diffusion_solver(subsection, context);
            sp.pressure = solver::parse_pressure_solver(subsection, context);
            sp.temperature = solver::parse_temperature_solver(subsection, context, true);
            sp.turbulence = solver::parse_turbulence_solver(subsection, context);
            sp.concentration = solver::parse_concentration_solver(subsection, context, true);
            sp.source = solver::parse_source_solver(subsection, context);
        } else if (sp.description == SolverTypes::NSTurbSolver) {
            sp.advection = solver::parse_advection_solver(subsection, context);
            sp.diffusion = solver::parse_diffusion_solver(subsection, context);
            sp.pressure = solver::parse_pressure_solver(subsection, context);
            sp.turbulence = solver::parse_turbulence_solver(subsection, context);
            sp.source = solver::parse_source_solver(subsection, context);
        } else if (sp.description == SolverTypes::PressureSolver) {
            sp.pressure = solver::parse_pressure_solver(subsection, context);
        } else {
            throw config_error(fmt::format("solver '{}' has no parsing implementation.", sp.description));
        }

        sp.solution = solver::parse_solution(subsection, context);
        return sp;
    }

    Settings_new parse_settings(const std::string &filename, const std::string &file_content) {
        tinyxml2::XMLDocument doc;
        doc.Parse(file_content.c_str());
        tinyxml2::XMLElement *root = doc.RootElement();
        auto solver_params = parse_solver_parameters(root);
        return {filename,
                parse_physical_parameters(root, solver_params.description),
                solver_params,
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
        auto settings = parse_settings(path.filename(), sstr.str());
#ifndef BENCHMARKING
        Utility::create_logger(settings.logging_parameters.level, settings.logging_parameters.file);  // create global logger
        auto logger = Utility::create_logger("XML File");
        logger->debug(sstr.str());
#endif
        return settings;
    }
}
