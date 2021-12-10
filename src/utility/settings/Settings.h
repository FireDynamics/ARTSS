/// \file       Settings.h
/// \brief      Access parameters of XML file
/// \date       Nov 17, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_SETTINGS_H
#define ARTSS_UTILITY_SETTINGS_H

#include "tinyxml2.h"
#include "BoundarySetting.h"
#include "ObstacleSetting.h"
#include "SurfaceSetting.h"
#include "../GlobalMacrosTypes.h"

#ifndef BENCHMARKING
#include <spdlog/logger.h>
#include <memory>
#endif

#include <map>
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <variant>


namespace Settings {
    class config_error : std::runtime_error { ;
    public:
        explicit config_error(const std::string &message) : std::runtime_error(message) {}
        explicit config_error(const char *message) : std::runtime_error(message) {};
    };
    struct domain_parameters {
        bool enable_computational_domain;
        real X1;
        real X2;
        real Y1;
        real Y2;
        real Z1;
        real Z2;
        real x1;
        real x2;
        real y1;
        real y2;
        real z1;
        real z2;
        size_t nx;
        size_t ny;
        size_t nz;
    };
    struct physical_parameters {
        real t_end;
        real dt;
        real nu;
        real beta;
        real g;
        real kappa;
    };
    struct visualisation_parameters {
        bool save_vtk;
        bool save_csv;
        size_t vtk_nth_plot;
        size_t csv_nth_plot;
    };
    struct logging_parameters {
        std::string file;
        std::string level;
    };
    struct Uniform {
        real value;
    };
    struct Zero {};
    struct GaussBubble {
        real u_lin;
        real v_lin;
        real w_lin;
        real x_shift;
        real y_shift;
        real z_shift;
        real l;
    };
    struct random_parameters {
        bool absolute;
        bool custom_seed;
        bool custom_steps;
        size_t seed;
        real step_size;
        real range;
    };
    struct initial_conditions_parameters {
        std::string usr_fct;
        bool random;
        std::variant<Uniform,Zero,GaussBubble> ic;
        struct random_parameters random_parameters;
    };
    struct boundary {
        std::string patch;
        std::string field_type;
        std::string boundary_condition;
        real value;
    };
    struct surface {
        real sx1;
        real sy1;
        real sz1;
        real sx2;
        real sy2;
        real sz2;
        struct boundary boundary_parameters;
    };
    struct surfaces_parameters {
        bool enabled;
        std::vector<struct surface> surfaces;
    };
    struct boundaries {
        std::vector<struct boundary> boundaries;
    };
    struct obstacle {
        std::string name;
        real ox1;
        real oy1;
        real oz1;
        real ox2;
        real oy2;
        real oz2;
        std::vector<struct boundary> boundaries;
    };
    struct obstacles_parameters {
        bool enabled;
        std::vector<struct obstacle> obstacles;
    };
    struct Settings_new {
        struct physical_parameters physical_parameters;
        //struct solver_parameters solver__parameters;
        struct domain_parameters domain_parameters;
        //struct adaption_parameters adaption_parameters;
        //struct boundaries_parameters boundaries_parameters;
        struct obstacles_parameters obstacles_parameters;
        struct surfaces_parameters surfaces_parameters;
        struct initial_conditions_parameters initial_conditions_parameters;
        struct visualisation_parameters visualisation_parameters;
        struct logging_parameters logging_parameters;
    };
    surfaces_parameters parse_surfaces_parameters(tinyxml2::XMLDocument &doc);
    obstacles_parameters parse_obstacles_parameters(tinyxml2::XMLDocument &doc);
    initial_conditions_parameters parse_initial_conditions_parameters(tinyxml2::XMLDocument &doc);
    visualisation_parameters parse_visualisation_parameters(tinyxml2::XMLDocument &doc);
    logging_parameters parse_logging_parameters(tinyxml2::XMLDocument &doc);
    domain_parameters parse_domain_parameters(tinyxml2::XMLDocument &doc);
    physical_parameters parse_physical_parameters(tinyxml2::XMLDocument &doc);
    Settings_new parse_settings(const std::string &file_content);
    Settings_new parse_settings_from_file(const std::filesystem::path &path);

class Settings {
public:
     explicit Settings(std::string path);
     void print_config() const;

     std::string get(std::string path) const;
     std::string sget(std::string path) const;
     void set(std::string path, std::string val) {
         sset(path, val);
#ifndef BENCHMARKING
         m_logger->debug(R"(set: "{}" to value: "{}")", path, val);
#endif
     }
     void sset(std::string path, std::string val) { m_proxy.insert({path, val}); }

     bool get_bool(std::string path) const { return get(path) == "Yes"; }
     int get_int(std::string path) const { return std::stoi(get(path)); }
     int get_size_t(std::string path) const { return std::stol(get(path)); }
     real get_real(std::string path) const { return real(std::stod(get(path))); }

     std::string get_filename() const { return filename; }

     std::vector<BoundarySetting> get_boundaries() const { return m_boundaries; }
     std::vector<ObstacleSetting> get_obstacles() const { return m_obstacles; }
     std::vector<SurfaceSetting> get_surfaces() const { return m_surfaces; }

 private:
     void read_config(std::string prefix, tinyxml2::XMLElement *elem);
     void read_boundaries(tinyxml2::XMLElement *elem);
     void read_obstacles(tinyxml2::XMLElement *elem);
     void read_surfaces(tinyxml2::XMLElement *elem);

     std::vector<BoundarySetting> m_boundaries;
     std::vector<ObstacleSetting> m_obstacles;
     std::vector<SurfaceSetting> m_surfaces;

     std::unordered_multimap<std::string, std::string> m_proxy;
     std::string filename;

#ifndef BENCHMARKING
     std::shared_ptr<spdlog::logger> m_logger;
#endif
};
}
#endif
