/// \file       Visual.cpp
/// \brief      coordinator of different visualisation methods
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Visual.h"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>

#include "../domain/DomainData.h"
#include "CSVWriter.h"
#include "VTKWriter.h"

Visual::Visual(const Settings::visualisation_parameters &settings, const Solution &solution, const std::string &file_name) :
        m_settings(settings),
        m_file_name(Utility::remove_extension(file_name)),
        m_solution(solution),
        m_has_analytical_solution(m_solution.has_analytical_solution()) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
}

void Visual::visualise(const FieldController &field_controller, real t) {
#ifndef BENCHMARKING
    m_logger->info("Visualise ...");
#endif
    auto domain_data = DomainData::getInstance();
    real dt = domain_data->get_physical_parameters().dt;
    real t_end = domain_data->get_physical_parameters().t_end;

    int n = static_cast<int> (std::round(t / dt));
#ifndef BENCHMARKING
    m_logger->debug("{}th time step", n);
#endif

    std::string file_name_numerical = create_file_name(m_file_name, n, false);
    std::string file_name_analytical = create_file_name(m_file_name, n, true);
    if (m_settings.save_vtk) {
        real vtk_plot = static_cast<real>(m_settings.vtk_nth_plot.value());
        if (fmod(n, vtk_plot) == 0 || t >= t_end) {
#ifndef BENCHMARKING
            m_logger->debug("write vtk file {}", n);
#endif
            VTKWriter::write_numerical(field_controller, file_name_numerical);
            if (m_has_analytical_solution) {
                VTKWriter::write_analytical(m_solution, file_name_analytical);
            }
        }
    }

    if (m_settings.save_csv) {
        real csv_plot = static_cast<real>(m_settings.csv_nth_plot.value());
        if (fmod(n, csv_plot) == 0 || t >= t_end) {
#ifndef BENCHMARKING
            m_logger->debug("write csv file {}", n);
#endif
            CSVWriter::write_numerical(field_controller, file_name_numerical);
            if (m_has_analytical_solution) {
                CSVWriter::write_analytical(m_solution, file_name_analytical);
            }
        }
    }
}

void Visual::write_csv(FieldController &field_controller, const std::string &file_name){
    field_controller.update_host();
    CSVWriter::write_numerical(field_controller, file_name);
}

void Visual::write_vtk(FieldController &field_controller, const std::string &file_name){
    field_controller.update_host();
    VTKWriter::write_numerical(field_controller, file_name);
}

void Visual::write_vtk_debug(FieldController &field_controller, const std::string &file_name){
    field_controller.update_host_debug();
    VTKWriter::write_numerical_debug(field_controller, file_name);
}

void Visual::initialise_grid(real *x_coords, real *y_coords, real *z_coords,
                             int Nx, int Ny, int Nz, real dx, real dy, real dz) {
    auto domain_data = DomainData::getInstance();
    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

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

std::string Visual::create_file_name(const std::string &file_name,
                                    int counter, bool analytical) {
    std::string fname = file_name;
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
