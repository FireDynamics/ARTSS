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
#include "../utility/Parameters.h"
#include "../DomainData.h"
#include "CSVWriter.h"
#include "VTKWriter.h"

Visual::Visual(const Solution &solution, bool has_analytical_solution) :
        m_solution(solution),
        m_has_analytical_solution(has_analytical_solution) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    auto params = Parameters::getInstance();
    m_filename = Utility::remove_extension(params->get_filename());

    m_save_csv = (params->get("visualisation/save_csv") == XML_TRUE);
    m_save_vtk = (params->get("visualisation/save_vtk") == XML_TRUE);

    m_dt = params->get_real("physical_parameters/dt");
    m_t_end = params->get_real("physical_parameters/t_end");
    if (m_save_csv) {
        m_csv_plots = params->get_int("visualisation/csv_nth_plot");
    }
    if (m_save_vtk) {
        m_vtk_plots = params->get_int("visualisation/vtk_nth_plot");
    }
}

void Visual::visualise(const FieldController &field_controller, real t) {
#ifndef BENCHMARKING
    m_logger->info("Visualise ...");
#endif
    int n = static_cast<int> (std::round(t / m_dt));
    std::string filename_numerical = create_filename(m_filename, n, false);
    std::string filename_analytical = create_filename(m_filename, n, true);
    if (m_save_vtk) {
        if (fmod(n, m_vtk_plots) == 0 || t >= m_t_end) {
            VTKWriter::write_numerical(field_controller, filename_numerical);
            if (m_has_analytical_solution) {
                VTKWriter::write_analytical(m_solution, filename_analytical);
            }
        }
    }

    if (m_save_csv) {
        if (fmod(n, m_csv_plots) == 0 || t >= m_t_end) {
            CSVWriter::write_numerical(field_controller, filename_numerical);
            if (m_has_analytical_solution) {
                CSVWriter::write_analytical(m_solution, filename_analytical);
            }
        }
    }
}

void Visual::write_csv(FieldController &field_controller, const std::string &filename){
    field_controller.update_host();
    CSVWriter::write_numerical(field_controller, filename);
}

void Visual::write_vtk(FieldController &field_controller, const std::string &filename){
    field_controller.update_host();
    VTKWriter::write_numerical(field_controller, filename);
}

void Visual::write_vtk_debug(FieldController &field_controller, const std::string &filename){
    field_controller.update_host_debug();
    VTKWriter::write_numerical_debug(field_controller, filename);
}

void Visual::initialise_grid(real *x_coords, real *y_coords, real *z_coords,
                             int Nx, int Ny, int Nz, real dx, real dy, real dz) {
    DomainData *domain = DomainData::getInstance();
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
