/// \file       Settings.h
/// \brief      Access parameters of XML file
/// \date       Nov 17, 2021
/// \author     von Mach
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.


#ifndef ARTSS_UTILITY_SETTINGS_H
#define ARTSS_UTILITY_SETTINGS_H

#include <map>
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include <variant>

#ifndef BENCHMARKING
#include <spdlog/logger.h>
#include <memory>
#endif

#include "tinyxml2.h"
#include "BoundarySetting.h"
#include "ObstacleSetting.h"
#include "SurfaceSetting.h"
#include "../GlobalMacrosTypes.h"
#include "../../domain/Coordinate.h"

using Map = std::map<std::string, std::string>;
using return_xml_data = std::tuple<const tinyxml2::XMLElement*,Map>;
namespace Settings {
    class config_error : std::runtime_error { ;
    public:
        explicit config_error(const std::string &message) : std::runtime_error(message) {}
        explicit config_error(const char *message) : std::runtime_error(message) {};
    };
    struct domain_parameters {
        bool enable_computational_domain;
        Coordinate<real> start_coords_CD;  // x1/y1/z1
        Coordinate<real> end_coords_CD;  // x2/y2/z2
        Coordinate<real> start_coords_PD;  // X1/Y1/Z1
        Coordinate<real> end_coords_PD;  // X2/Y2/Z2
        Coordinate<size_t> number_of_inner_cells;  // nx/ny/nz
    };
    struct physical_parameters {
        real t_end;
        real dt;
        std::optional<real> nu;
        real beta;
        real g;
        std::optional<real> kappa;  // thermal diffusion
        std::optional<real> gamma;  // concentration diffusion
        std::optional<real> rhoa;  // fluid density for buoyancy
    };
    struct visualisation_parameters {
        bool save_vtk;
        bool save_csv;
        std::optional<size_t> vtk_nth_plot;
        std::optional<size_t> csv_nth_plot;
    };
    struct logging_parameters {
        std::string file;
        std::string level;
    };
    namespace initial_conditions {
        struct uniform {
            real value;
        };
        struct drift {
            Coordinate<real> velocity_lin;
            real pa;
        };
        struct gauss_bubble {
            Coordinate<real> velocity_lin;
            Coordinate<real> shift;
            real l;
        };
        struct exp_sinus_prod {
            real l;
        };
        struct hat {
            Coordinate<real> start_coords;
            Coordinate<real> end_coords;
            real val_in;
            real val_out;
        };
        struct jet {
            Coordinate<real> start_coords;
            Coordinate<real> end_coords;
            real value;
            CoordinateAxis dir;
        };
        struct vortex {
            Coordinate<real> velocity_lin;
            real pa;
            real rhoa;
        };
        struct mc_dermott {
            real A;
        };
        struct layers_temperature {
            size_t number_of_layers;
            std::vector<real> borders;
            std::vector<real> values;
            CoordinateAxis dir;
        };
        struct sin_sin_sin {
            real l;
        };
        struct beltrami {
            real a;
            real d;
        };
    }
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
        std::optional<std::variant<initial_conditions::uniform,
                     initial_conditions::beltrami,
                     initial_conditions::drift,
                     initial_conditions::exp_sinus_prod,
                     initial_conditions::gauss_bubble,
                     initial_conditions::hat,
                     initial_conditions::jet,
                     initial_conditions::layers_temperature,
                     initial_conditions::mc_dermott,
                     initial_conditions::sin_sin_sin,
                     initial_conditions::vortex>> ic;
        std::optional<struct random_parameters> random_parameters;
    };
    struct boundary {
        std::vector<Patch> patch;
        std::vector<FieldType> field_type;
        BoundaryCondition boundary_condition;
        std::optional<real> value;
    };
    struct surface {
        Coordinate<real> start_coords;
        Coordinate<real> end_coords;
        struct boundary boundary_parameters;
    };
    struct surfaces_parameters {
        bool enabled;
        std::vector<struct surface> surfaces;
    };
    struct boundaries_parameters {
        std::vector<struct boundary> boundaries;
    };
    struct obstacle {
        std::string name;
        Coordinate<real> start_coords;
        Coordinate<real> end_coords;
        std::vector<struct boundary> boundaries;
    };
    struct obstacles_parameters {
        bool enabled;
        std::vector<struct obstacle> obstacles;
    };
    namespace adaption_classes {
        struct layers {
            size_t number_of_buffer_cells;
            real check_value;
            size_t time_step;
            size_t expansion_size;
        };
        struct vortex {
            Coordinate<real> velocity;
            bool reduction;
            std::vector<CoordinateAxis> dir;
            size_t buffer;
            size_t threshold;
        };
        struct data_extraction {
            bool has_time_measuring;
            bool has_data_extraction_endresult;
            bool has_data_extraction_after;
            bool has_data_extraction_before;
        };
    }
    struct adaption_parameters {
        bool enabled;
        std::optional<std::string> class_name;
        std::optional<std::variant<adaption_classes::layers, adaption_classes::vortex>> adaption_class;
        bool has_data_extraction;
        std::optional<adaption_classes::data_extraction> data_extraction;
    };
    namespace solver {
        struct advection_solver {
            std::string type;
            std::vector<FieldType> fields;
        };
        namespace diffusion_solvers {
            struct jacobi {
                size_t max_iter;
                real tol_res;
                real w;
            };
            struct colored_gauss_seidel {
                size_t max_iter;
                real tol_res;
                real w;
            };
        }
        struct diffusion_solver {
            std::string type;
            std::vector<FieldType> fields;
            std::optional<std::variant<diffusion_solvers::jacobi, diffusion_solvers::colored_gauss_seidel>> solver;
        };
        namespace turbulence_solvers {
            struct const_smagorinsky {
                real cs;
            };
        }
        struct turbulence_solver {
            std::string type;
            std::optional<turbulence_solvers::const_smagorinsky> solver;
        };
        namespace source_solvers {
            struct buoyancy {
                bool use_init_values;
                std::optional<real> ambient_temperature_value;
            };
            struct uniform {
                Coordinate<real> velocity_value;
            };
        }
        struct source_solver {
            std::string type;
            std::string force_fct;
            std::vector<CoordinateAxis> direction;
            std::variant<source_solvers::buoyancy,source_solvers::uniform> force_function;
        };
        namespace pressure_solvers {
            struct vcycle_mg {
                size_t n_level;
                size_t n_cycle;
                size_t max_cycle;
                size_t n_relax;
                real tol_res;
                struct diffusion_solver smoother;
            };
        }
        struct pressure_solver {
            std::string type;
            FieldType field;
            struct pressure_solvers::vcycle_mg solver;
        };
        namespace sources {
            struct gauss {
                real heat_release_rate;
                real heat_capacity;
                Coordinate<real> position;
                Coordinate<real> dimension;
                real tau;
            };
            struct cube {
                Coordinate<real> coords_start;
                Coordinate<real> coords_end;
                real value;
            };
        }
        struct temperature_source {
            std::string type;
            std::vector<CoordinateAxis> dir;
            std::string temp_fct;
            bool dissipation;
            bool random;
            struct random_parameters random_parameters;
            std::variant<sources::gauss, sources::cube> temp_function;
        };
        struct temperature_solver {
            advection_solver advection;
            diffusion_solver diffusion;
            bool has_turbulence;
            std::optional<real> prandtl_number;
            temperature_source source;
        };
        struct concentration_source {
            std::string type;
            std::vector<CoordinateAxis> dir;
            std::string con_fct;
            bool random;
            struct random_parameters random_parameters;
            std::variant<sources::gauss, sources::cube> con_function;
        };
        struct concentration_solver {
            advection_solver advection;
            diffusion_solver diffusion;
            bool has_turbulence;
            std::optional<real> turbulent_schmidt_number;
            concentration_source source;
        };
        struct solution {
            bool analytical_solution;
            std::optional<real> solution_tolerance;
        };
    }
    struct solver_parameters {
        std::string description;
        solver::advection_solver advection;
        solver::diffusion_solver diffusion;
        solver::turbulence_solver turbulence;
        solver::source_solver source;
        solver::pressure_solver pressure;
        solver::temperature_solver temperature;
        solver::concentration_solver concentration;
        solver::solution solution;
    };
    struct Settings_new {
        std::string filename;
        struct physical_parameters physical_parameters;
        struct solver_parameters solver_parameters;
        struct domain_parameters domain_parameters;
        struct adaption_parameters adaption_parameters;
        struct boundaries_parameters boundaries_parameters;
        struct obstacles_parameters obstacles_parameters;
        struct surfaces_parameters surfaces_parameters;
        struct initial_conditions_parameters initial_conditions_parameters;
        struct visualisation_parameters visualisation_parameters;
        struct logging_parameters logging_parameters;
    };
    random_parameters parse_random_parameters(const tinyxml2::XMLElement *head, const std::string &parent_context);
    solver_parameters parse_solver_parameters(const tinyxml2::XMLElement *root);
    surfaces_parameters parse_surfaces_parameters(const tinyxml2::XMLElement *root);
    obstacles_parameters parse_obstacles_parameters(const tinyxml2::XMLElement *root);
    adaption_parameters parse_adaption_parameters(const tinyxml2::XMLElement *root);
    boundaries_parameters parse_boundaries_parameters(const tinyxml2::XMLElement *root);
    initial_conditions_parameters parse_initial_conditions_parameters(const tinyxml2::XMLElement *root);
    visualisation_parameters parse_visualisation_parameters(const tinyxml2::XMLElement *root);
    logging_parameters parse_logging_parameters(const tinyxml2::XMLElement *root);
    domain_parameters parse_domain_parameters(const tinyxml2::XMLElement *root);
    physical_parameters parse_physical_parameters(const tinyxml2::XMLElement *root, const std::string &solver_description);
    Settings_new parse_settings(const std::string &filename, const std::string &file_content);
    Settings_new parse_settings_from_file(const std::filesystem::path &path);

class Settings {
 public:
     Settings();
     explicit Settings(const std::string& path);
     void print_config() const;

     std::string get(const std::string& path) const;
     std::string sget(const std::string& path) const;
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
