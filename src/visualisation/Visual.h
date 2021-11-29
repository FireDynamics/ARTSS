/// \file       Visual.h
/// \brief      coordinator of different visualisation methods
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_VISUALISATION_VISUAL_H_
#define ARTSS_VISUALISATION_VISUAL_H_

#include "../analysis/Solution.h"
#include "../interfaces/ISolver.h"
#include "../field/FieldController.h"
#include "../utility/Utility.h"
#include "../utility/settings/Settings.h"

class Visual {
public:
    Visual(Settings::Settings const &settings, Solution const &solution, bool has_analytical_solution);

    void visualise(const FieldController &field_controller, real t);

    static void initialise_grid(real *x_coords, real *y_coords, real *z_coords, int Nx, int Ny, int Nz, real dx, real dy, real dz);

    static void write_csv(FieldController &field_controller, const std::string& filename);
    static void write_vtk(FieldController &field_controller, const std::string& filename);
    static void write_vtk_debug(FieldController &field_controller, const std::string& filename);

private:
    Settings::Settings const &m_settings;

    std::string m_filename;
    const Solution &m_solution;
    bool m_save_csv = false;
    int m_csv_plots = 0;
    bool m_save_vtk = false;
    int m_vtk_plots = 0;

    static std::string create_filename(const std::string &filename, int counter, bool analytical);

    bool m_has_analytical_solution = false;
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};

#endif /* ARTSS_VISUALISATION_VISUAL_H_ */
