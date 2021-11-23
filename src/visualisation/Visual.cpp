/// \file       Visual.cpp
/// \brief      coordinator of different visualisation methods
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>
#include <iomanip>
#include <string>
#include <utility>

#include "Visual.h"
#include "../Domain.h"
#include "CSVWriter.h"
#include "VTKWriter.h"

Visual::Visual(Settings const &settings, Solution const &solution, bool has_analytical_solution) :
        m_solution(solution),
        m_has_analytical_solution(has_analytical_solution) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(settings, typeid(this).name());
#endif
    m_filename = Utility::remove_extension(settings.get_filename());

    m_save_csv = settings.get_bool("visualisation/save_csv");
    m_save_vtk = settings.get_bool("visualisation/save_vtk");

    m_dt = settings.get_real("physical_parameters/dt");
    m_t_end = settings.get_real("physical_parameters/t_end");
    if (m_save_csv) {
        m_csv_plots = settings.get_int("visualisation/csv_nth_plot");
    }
    if (m_save_vtk) {
        m_vtk_plots = settings.get_int("visualisation/vtk_nth_plot");
    }
}

void Visual::visualise(const FieldController &field_controller, real t) {
#ifndef BENCHMARKING
    m_logger->info("Visualise ...");
#endif
    int n = static_cast<int> (std::round(t / m_dt));
    std::string filename = create_filename(m_filename, n, false);
    if (m_save_vtk) {
        if (fmod(n, m_vtk_plots) == 0 || t >= m_t_end) {
            VTKWriter::write_numerical(field_controller, filename);
            if (m_has_analytical_solution) {
                VTKWriter::write_analytical(m_solution, filename);
            }
        }
    }

    if (m_save_csv) {
        if (fmod(n, m_csv_plots) == 0 || t >= m_t_end) {
            CSVWriter::write_numerical(field_controller, filename);
            if (m_has_analytical_solution) {
                CSVWriter::write_analytical(m_solution, filename);
            }
        }
    }
}

void Visual::write_csv(FieldController &field_controller, const std::string &filename){
    //TODO method update host/device whatever
    field_controller.get_field_u().update_host();
    field_controller.get_field_v().update_host();
    field_controller.get_field_w().update_host();
    field_controller.get_field_p().update_host();
    field_controller.get_field_rhs().update_host();
    field_controller.get_field_T().update_host();
    field_controller.get_field_concentration().update_host();
    field_controller.get_field_source_T().update_host();
    field_controller.get_field_source_concentration().update_host();
    field_controller.get_field_nu_t().update_host();
    CSVWriter::write_numerical(field_controller, filename);
}

void Visual::write_vtk(FieldController &field_controller, const std::string &filename){
    //TODO method update host/device whatever
    field_controller.get_field_u().update_host();
    field_controller.get_field_v().update_host();
    field_controller.get_field_w().update_host();
    field_controller.get_field_p().update_host();
    field_controller.get_field_rhs().update_host();
    field_controller.get_field_T().update_host();
    field_controller.get_field_concentration().update_host();
    field_controller.get_field_source_T().update_host();
    field_controller.get_field_source_concentration().update_host();
    field_controller.get_field_nu_t().update_host();
    VTKWriter::write_numerical(field_controller, filename);
}

void Visual::write_vtk_debug(FieldController &field_controller, const std::string &filename){
    //TODO method update host/device whatever
    field_controller.get_field_u().update_host();
    field_controller.get_field_v().update_host();
    field_controller.get_field_w().update_host();
    field_controller.get_field_p().update_host();
    field_controller.get_field_rhs().update_host();
    field_controller.get_field_T().update_host();
    field_controller.get_field_concentration().update_host();
    field_controller.get_field_source_T().update_host();
    field_controller.get_field_source_concentration().update_host();
    field_controller.get_field_nu_t().update_host();
    field_controller.get_field_force_x().update_host();
    field_controller.get_field_force_y().update_host();
    field_controller.get_field_force_z().update_host();
    field_controller.get_field_kappa().update_host();
    field_controller.get_field_gamma().update_host();
    VTKWriter::write_numerical_debug(field_controller, filename);
}

void Visual::initialise_grid(real *x_coords, real *y_coords, real *z_coords,
                             int Nx, int Ny, int Nz, real dx, real dy, real dz) {
    Domain *domain = Domain::getInstance();
    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    // Initialize grid
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                x_coords[index] = (X1 + (i - 0.5) * dx);
                y_coords[index] = (Y1 + (j - 0.5) * dy);
                z_coords[index] = (Z1 + (k - 0.5) * dz);
            }
        }
    }
}

std::string Visual::create_filename(const std::string &filename,
                                    int counter, bool analytical) {
    std::string fname = std::move(filename);
    if (analytical) {
        fname.append("_ana_");
    } else {
        fname.append("_num_");
    }
    std::ostringstream tstep;
    tstep << std::setw(6) << std::setfill('0') << counter;
    fname.append(tstep.str());
    return fname;
}
