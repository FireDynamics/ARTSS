/// \file       CSVWriter.cpp
/// \brief      class to write out csv files
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <fstream>
#include <iomanip>
#include "CSVWriter.h"
#include "../Domain.h"
#include "Visual.h"

static std::string ending = ".csv";
const static char delimiter = ',';

void CSVWriter::write_numerical(FieldController *field_controller, const std::string& filename) {
    auto u = field_controller->get_field_u_data();
    auto v = field_controller->get_field_v_data();
    auto w = field_controller->get_field_w_data();
    auto p = field_controller->get_field_p_data();
    auto div = field_controller->get_field_rhs_data();
    auto T = field_controller->get_field_T_data();
    auto C = field_controller->get_field_concentration_data();
    auto s = field_controller->get_field_sight_data();
    auto nu_t = field_controller->get_field_nu_t_data();
    auto S_T = field_controller->get_field_source_T_data();
    CSVWriter::csvPrepareAndWrite((filename + ending).c_str(), u, v, w, p, div, T, C, s, nu_t, S_T);
}

void CSVWriter::write_analytical(Solution *solution, const std::string& filename) {
    auto u = solution->GetU_data();
    auto v = solution->GetV_data();
    auto w = solution->GetW_data();
    auto p = solution->GetP_data();
    auto T = solution->GetT_data();
    CSVWriter::csvPrepareAndWrite((filename + ending).c_str(), u, v, w, p, T);
}

void CSVWriter::csvPrepareAndWrite(const char *filename, real *u, real *v, real *w, real *p, real *div, real *T, real *C, real *s, real *nu_t, real *S_T) {
    Domain *domain = Domain::getInstance();
    int size = static_cast<int>(domain->get_size());

    // Initialize variables
    int size_vars = 10; // Number of variables
    const char *var_names[] = {"x-velocity (m/s)", "y-velocity (m/s)", "z-velocity (m/s)",
                               "pressure (kg/(m s^2))",
                               "divergence (1/s)",
                               "temperature (Celsius)",
                               "concentration (g/m^3)",
                               "sight",
                               "turb viscosity (m^2/s)",
                               "temperature source (K/s)"};

    // velocities
    auto *u_vel = new real[size];
    auto *v_vel = new real[size];
    auto *w_vel = new real[size];
    // pressure
    auto *pres = new real[size];
    // divergence
    auto *vel_div = new real[size];
    // temperature
    auto *Temp = new real[size];
    // smoke concentration
    auto *Con = new real[size];
    // boundary sight
    auto *Sight = new real[size];
    // turbulent viscosity
    auto *turb_visc = new real[size];
    // energy source
    auto *source_T = new real[size];

    // Summarize pointers to variables in an array
    real* fields[] = {u, v, w, p, div, T, C, s, nu_t, S_T};
    CSVWriter::csv_write(filename, fields, size_vars, var_names);

    // Clean up
    delete[] (u_vel);
    delete[] (v_vel);
    delete[] (w_vel);
    delete[] (pres);
    delete[] (vel_div);
    delete[] (Temp);
    delete[] (Con);
    delete[] (Sight);
    delete[] (turb_visc);
    delete[] (source_T);
}

void CSVWriter::csvPrepareAndWrite(const char *filename, real *u, real *v, real *w, real *p, real *T) {
    Domain *domain = Domain::getInstance();
    int size = static_cast<int>(domain->get_size());

    // Initialize variables
    int size_vars = 5; // Number of variables
    const char *var_names[] = {"x-velocity (m/s)", "y-velocity (m/s)", "z-velocity (m/s)",
                               "pressure (kg/(m s^2))",
                               "temperature (Celsius)"};

    // velocities
    auto *u_vel = new real[size];
    auto *v_vel = new real[size];
    auto *w_vel = new real[size];
    // pressure
    auto *pres = new real[size];
    // temperature
    auto *Temp = new real[size];

    // Summarize pointers to variables in an array
    real* fields[] = {u, v, w, p, T};

    CSVWriter::csv_write(filename, fields, size_vars, var_names);

    // Clean up
    delete[] (u_vel);
    delete[] (v_vel);
    delete[] (w_vel);
    delete[] (pres);
    delete[] (Temp);
}

void CSVWriter::csv_write(const char *filename, real **vars, int size_vars, const char **var_names) {
    Domain *domain = Domain::getInstance();

    int Nx = static_cast<int>(domain->get_Nx());
    int Ny = static_cast<int>(domain->get_Ny());
    int Nz = static_cast<int>(domain->get_Nz());

    int size = static_cast<int>(domain->get_size());

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    const char *coord_names[] = {"i", "j", "k", "index", "x-coords (m)", "y-coords (m)", "z-coords (m)"};

    auto *x_centres = new float[size];
    auto *y_centres = new float[size];
    auto *z_centres = new float[size];
    float *coords[] = {x_centres, y_centres, z_centres};
    Visual::initialise_grid(x_centres, y_centres, z_centres, Nx, Ny, Nz, dx, dy, dz);

    // write data to csv
    std::ofstream outputFile;
    outputFile.open(filename, std::ofstream::out);

    // var_names as column titles
    for (auto & coord_name : coord_names) {
        outputFile << coord_name << delimiter;
    }
    for (int i = 0; i < size_vars - 1; i++) {
        outputFile << var_names[i] << delimiter;
    }
    // last column
    outputFile << var_names[size_vars - 1] << std::endl;

    // write variables to csv
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                outputFile << i << delimiter
                           << j << delimiter
                           << k << delimiter
                           << index << delimiter
                           << std::setprecision(24);

                for (auto & coord : coords) {
                    outputFile << coord[index] << delimiter;
                }
                for (int v = 0; v < size_vars - 1; v++) {
                    outputFile << std::setprecision(24) << vars[v][index] << delimiter;
                }
                outputFile << std::setprecision(24) << vars[size_vars - 1][index] << std::endl;
            }
        }
    }
    outputFile.close();

    delete[] (x_centres);
    delete[] (y_centres);
    delete[] (z_centres);
}

void CSVWriter::write_data(std::string *data_titles,
        real **data, size_t size_data, const std::string& filename) {
    const char *var_names[size_data];
    for (size_t i = 0; i < size_data; i++){
        var_names[i] = data_titles[i].c_str();
    }

    CSVWriter::csv_write((filename+ending).c_str(), data, size_data, var_names);
}

