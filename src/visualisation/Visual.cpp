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
#include "../Domain.h"
#include "CSVWriter.h"
#include "VTKWriter.h"

Visual::Visual(Solution *solution) {
    auto params = Parameters::getInstance();
    m_filename = remove_extension(params->get_filename());

    m_save_csv = (params->get("visualisation/save_csv") == "Yes");
    m_save_vtk = (params->get("visualisation/save_vtk") == "Yes");

    m_dt = params->get_real("physical_parameters/dt");
    m_t_end = params->get_real("physical_parameters/t_end");
    if (m_save_csv) {
        m_csv_plots = params->get_int("visualisation/csv_nth_plot");
    }
    if (m_save_vtk) {
        m_vtk_plots = params->get_int("visualisation/vtk_nth_plot");
    }
    m_has_analytical_solution = (params->get("solver/solution/available") == "Yes");
    m_solution = solution;
}

void Visual::visualise(ISolver *solver, const real t) {
    int n = static_cast<int> (std::round(t / m_dt));
    if (m_has_analytical_solution) {
        m_solution->calc_analytical_solution(t);
    }

    std::string filename = create_filename(m_filename, static_cast<int>(std::round(t / m_dt)), false);
    if (m_save_vtk) {
        if (fmod(n, m_vtk_plots) == 0 || t >= m_t_end) {
            VTKWriter::write_numerical(solver, filename);
            if (m_has_analytical_solution) {
                VTKWriter::write_analytical(m_solution, create_filename(m_filename, static_cast<int>(t/m_dt), true));
            }
        }
    }

    if (m_save_csv) {
        if (fmod(n, m_csv_plots) == 0 || t >= m_t_end) {
            CSVWriter::write_numerical(solver, filename);
            if (m_has_analytical_solution) {
                CSVWriter::write_analytical(m_solution, create_filename(m_filename, static_cast<int>(t/m_dt), true));
            }
        }
    }
}

void Visual::initialise_grid(float *x_coords, float *y_coords, float *z_coords, int Nx, int Ny, int Nz, real dx, real dy, real dz) {
    Domain *domain = Domain::getInstance();
    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    // Initialize grid
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                x_coords[index] = static_cast<float> (X1 + (i - 0.5) * dx);
                y_coords[index] = static_cast<float> (Y1 + (j - 0.5) * dy);
                z_coords[index] = static_cast<float> (Z1 + (k - 0.5) * dz);
            }
        }
    }
}

void Visual::prepare_fields(read_ptr *fields, float **vars, int size) {
    Domain *domain = Domain::getInstance();

    int domain_size = static_cast<int>(domain->get_size());

    // Cast variables to floats
    for (int index = 0; index < domain_size; index++) {
        for (int v = 0; v < size; v++) {
            vars[v][index] = static_cast<float>(fields[v][index]);
        }
    }
}

std::string Visual::create_filename(std::string filename, int counter, bool analytical) {
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

//================================= Remove extension ==================================
// ***************************************************************************************
/// \brief  Removes extension from filename
/// \param  filename    xml-file name (via argument)
// ***************************************************************************************
std::string Visual::remove_extension(const std::string &filename) {
    size_t lastdot = filename.find_last_of('.');
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}
