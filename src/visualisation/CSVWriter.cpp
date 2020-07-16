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

void CSVWriter::write_numerical(ISolver *solver, const std::string& filename) {
    auto u = solver->get_u();
    auto v = solver->get_v();
    auto w = solver->get_w();
    auto p = solver->get_p();
    auto div = solver->get_rhs();
    auto T = solver->get_T();
    auto C = solver->get_concentration();
    auto s = solver->get_sight();
    auto nu_t = solver->get_nu_t();
    auto S_T = solver->get_S_T();
    CSVWriter::csvPrepareAndWrite((filename + ending).c_str(), u, v, w, p, div, T, C, s, nu_t, S_T);
}

void CSVWriter::write_analytical(Solution *solution, const std::string& filename) {
    auto u = solution->GetU();
    auto v = solution->GetV();
    auto w = solution->GetW();
    auto p = solution->GetP();
    auto T = solution->GetT();
    CSVWriter::csvPrepareAndWrite((filename + ending).c_str(), u, v, w, p, T);
}

void CSVWriter::csvPrepareAndWrite(const char *filename, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr div, read_ptr T, read_ptr C, read_ptr s, read_ptr nu_t, read_ptr S_T) {
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
    auto *u_vel = new float[size];
    auto *v_vel = new float[size];
    auto *w_vel = new float[size];
    // pressure
    auto *pres = new float[size];
    // divergence
    auto *vel_div = new float[size];
    // temperature
    auto *Temp = new float[size];
    // smoke concentration
    auto *Con = new float[size];
    // boundary sight
    auto *Sight = new float[size];
    // turbulent viscosity
    auto *turb_visc = new float[size];
    // energy source
    auto *source_T = new float[size];

    // Summarize pointers to variables in an array
    float *vars[] = {static_cast<float *> (u_vel),
                     static_cast<float *> (v_vel),
                     static_cast<float *> (w_vel),
                     static_cast<float *> (pres),
                     static_cast<float *> (vel_div),
                     static_cast<float *> (Temp),
                     static_cast<float *> (Con),
                     static_cast<float *> (Sight),
                     static_cast<float *> (turb_visc),
                     static_cast<float *> (source_T)};
    read_ptr fields[] = {u, v, w, p, div, T, C, s, nu_t, S_T};
    Visual::prepare_fields(fields, vars, size_vars);

    CSVWriter::csv_write(filename, vars, size_vars, var_names);

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

void CSVWriter::csvPrepareAndWrite(const char *filename, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr T) {
    Domain *domain = Domain::getInstance();
    int size = static_cast<int>(domain->get_size());

    // Initialize variables
    int size_vars = 5; // Number of variables
    const char *var_names[] = {"x-velocity (m/s)", "y-velocity (m/s)", "z-velocity (m/s)",
                               "pressure (kg/(m s^2))",
                               "temperature (Celsius)"};

    // velocities
    auto *u_vel = new float[size];
    auto *v_vel = new float[size];
    auto *w_vel = new float[size];
    // pressure
    auto *pres = new float[size];
    // temperature
    auto *Temp = new float[size];

    // Summarize pointers to variables in an array
    float *vars[] = {static_cast<float *> (u_vel),
                     static_cast<float *> (v_vel),
                     static_cast<float *> (w_vel),
                     static_cast<float *> (pres),
                     static_cast<float *> (Temp)};
    read_ptr fields[] = {u, v, w, p, T};

    Visual::prepare_fields(fields, vars, size_vars);
    CSVWriter::csv_write(filename, vars, size_vars, var_names);

    // Clean up
    delete[] (u_vel);
    delete[] (v_vel);
    delete[] (w_vel);
    delete[] (pres);
    delete[] (Temp);
}

void CSVWriter::csv_write(const char *filename, float **vars, int size_vars, const char **var_names) {
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
                           << std::setprecision(16);

                for (auto & coord : coords) {
                    outputFile << coord[index] << delimiter;
                }
                for (int v = 0; v < size_vars - 1; v++) {
                    outputFile << vars[v][index] << delimiter;
                }
                outputFile << vars[size_vars - 1][index] << std::endl;
            }
        }
    }
    outputFile.close();

    delete[] (x_centres);
    delete[] (y_centres);
    delete[] (z_centres);
}

