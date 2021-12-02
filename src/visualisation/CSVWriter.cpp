/// \file       CSVWriter.cpp
/// \brief      class to write out csv files
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef BENCHMARKING
#define FMT_STRING_ALIAS 1
#include <fmt/os.h>
#include <fmt/compile.h>
#endif

#include <fstream>

#include "CSVWriter.h"

static std::string ending = ".csv";
const static std::string_view delimiter = ",";

void CSVWriter::write_numerical(const FieldController &field_controller, const std::string &filename) {
    auto u = field_controller.get_field_u_data();
    auto v = field_controller.get_field_v_data();
    auto w = field_controller.get_field_w_data();
    auto p = field_controller.get_field_p_data();
    auto div = field_controller.get_field_rhs_data();
    auto T = field_controller.get_field_T_data();
    auto C = field_controller.get_field_concentration_data();
    auto sight = field_controller.get_field_sight_data();
    auto nu_t = field_controller.get_field_nu_t_data();
    auto source_T = field_controller.get_field_source_T_data();
    CSVWriter::csv_prepare_and_write(filename + ending, u, v, w, p, div, T, C, sight, nu_t, source_T);
}

void CSVWriter::write_analytical(const Solution &solution, const std::string &filename) {
    auto u = solution.get_return_ptr_data_u();
    auto v = solution.get_return_ptr_data_v();
    auto w = solution.get_return_ptr_data_w();
    auto p = solution.get_return_ptr_data_p();
    auto T = solution.get_return_ptr_data_T();
    CSVWriter::csv_prepare_and_write(filename + ending, u, v, w, p, T);
}

void CSVWriter::csv_prepare_and_write(const std::string &filename,
                                      return_ptr u, return_ptr v, return_ptr w,
                                      return_ptr p,
                                      return_ptr div,
                                      return_ptr T,
                                      return_ptr C,
                                      return_ptr sight,
                                      return_ptr nu_t,
                                      return_ptr source_T) {
    // Initialize variables
    int size_vars = 10; // Number of variables
    std::vector<std::string> var_names = {"x-velocity (m/s)", "y-velocity (m/s)", "z-velocity (m/s)",
                                           "pressure (kg/(m s^2))",
                                           "divergence (1/s)",
                                           "temperature (Celsius)",
                                           "concentration (g/m^3)",
                                           "sight",
                                           "turb viscosity (m^2/s)",
                                           "temperature source (K/s)"};

    return_ptr fields[] = {u, v, w, p, div, T, C, sight, nu_t, source_T};
    CSVWriter::csv_write(filename, fields, size_vars, var_names);
}

void CSVWriter::csv_prepare_and_write(const std::string &filename, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr T) {
    // Initialize variables
    int size_vars = 5;  // Number of variables
    std::vector<std::string> var_names = {"x-velocity (m/s)", "y-velocity (m/s)", "z-velocity (m/s)",
                                           "pressure (kg/(m s^2))",
                                           "temperature (Celsius)"};

    // Summarize pointers to variables in an array
    return_ptr fields[] = {u, v, w, p, T};

    CSVWriter::csv_write(filename, fields, size_vars, var_names);
}

void CSVWriter::csv_write(const std::string &filename, return_ptr *vars, int size_vars, const std::vector<std::string> &var_names) {
#ifndef BENCHMARKING
    DomainData *domain = DomainData::getInstance();

    int Nx = static_cast<int>(domain->get_Nx());
    int Ny = static_cast<int>(domain->get_Ny());
    int Nz = static_cast<int>(domain->get_Nz());

    int size = static_cast<int>(domain->get_size());

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    std::vector<std::string> coord_names = {"i", "j", "k", "index",
                                            "x-coords (m)", "y-coords (m)", "z-coords (m)"};

    auto *x_centres = new real[size];
    auto *y_centres = new real[size];
    auto *z_centres = new real[size];
    real *coords[] = {x_centres, y_centres, z_centres};
    Visual::initialise_grid(x_centres, y_centres, z_centres,
                            Nx, Ny, Nz,
                            dx, dy, dz);

    // write data to csv
    std::ofstream output_file(filename, std::ofstream::binary);

    // var_names as column titles
    output_file << fmt::format("{},", fmt::join(coord_names, delimiter));
    output_file << fmt::format("{}\n", fmt::join(var_names, delimiter));

    // write variables to csv
    std::string result;
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                result = fmt::format(FMT_COMPILE("{},{},{},{},"),
                                                i, j, k, index);

                result.append(fmt::format(FMT_COMPILE("{},{},{},"),
                                                coords[0][index],
                                                coords[1][index],
                                                coords[2][index]));

                for (int v = 0; v < size_vars - 1; v++) {
                    result.append(fmt::format(FMT_COMPILE("{},"),
                                                vars[v][index]));
                }
                result.append(fmt::format(FMT_COMPILE("{}\n"),
                                                vars[size_vars - 1][index]));

                output_file << result;
            }
        }
    }

    delete[] x_centres;
    delete[] y_centres;
    delete[] z_centres;
#endif
}

