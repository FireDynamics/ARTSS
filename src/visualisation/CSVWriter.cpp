/// \file       CSVWriter.cpp
/// \brief      class to write out csv files
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <fmt/os.h>
#include "CSVWriter.h"

static std::string ending = ".csv";
const static std::string_view delimiter = ",";

void CSVWriter::write_numerical(FieldController *field_controller, const std::string& filename) {
    auto u = field_controller->field_u->data;
    auto v = field_controller->field_v->data;
    auto w = field_controller->field_w->data;
    auto p = field_controller->field_p->data;
    auto div = field_controller->field_rhs->data;
    auto T = field_controller->field_T->data;
    auto C = field_controller->field_concentration->data;
    auto s = field_controller->sight->data;
    auto nu_t = field_controller->field_nu_t->data;
    auto S_T = field_controller->field_source_T->data;
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

    real* fields[] = {u, v, w, p, div, T, C, s, nu_t, S_T};
    CSVWriter::csv_write(filename, fields, size_vars, var_names);
}

void CSVWriter::csvPrepareAndWrite(const char *filename, real *u, real *v, real *w, real *p, real *T) {
    // Initialize variables
    int size_vars = 5; // Number of variables
    std::vector<std::string> var_names = {"x-velocity (m/s)", "y-velocity (m/s)", "z-velocity (m/s)",
                                           "pressure (kg/(m s^2))",
                                           "temperature (Celsius)"};

    // Summarize pointers to variables in an array
    real* fields[] = {u, v, w, p, T};

    CSVWriter::csv_write(filename, fields, size_vars, var_names);
}

void CSVWriter::csv_write(const char *filename, real **vars, int size_vars, std::vector<std::string> var_names) {
    Domain *domain = Domain::getInstance();

    int Nx = static_cast<int>(domain->get_Nx());
    int Ny = static_cast<int>(domain->get_Ny());
    int Nz = static_cast<int>(domain->get_Nz());

    int size = static_cast<int>(domain->get_size());

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    std::vector<std::string> coord_names = {"i", "j", "k", "index",
                                            "x-coords (m)", "y-coords (m)", "z-coords (m)"};

    auto *x_centres = new float[size];
    auto *y_centres = new float[size];
    auto *z_centres = new float[size];
    float *coords[] = {x_centres, y_centres, z_centres};
    Visual::initialise_grid(x_centres, y_centres, z_centres, Nx, Ny, Nz, dx, dy, dz);

    // write data to csv
    auto output_file = fmt::output_file(filename);

    // var_names as column titles
    output_file.print("{}", fmt::join(coord_names, ","));
    output_file.print("{}\n", fmt::join(var_names, delimiter));

    // write variables to csv
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                output_file.print("{0}{4}{1}{4}{2}{4}{3}{4}", i, j, k, index, delimiter);

                output_file.print("{0}{3}{1}{3}{2}{3}", coords[0][index],
                                                          coords[1][index],
                                                          coords[2][index],
                                                          delimiter);

                for (int v = 0; v < size_vars - 1; v++) {
                    output_file.print("{}{}", vars[v][index], delimiter);
                }
                output_file.print("{}\n", vars[size_vars - 1][index]);
            }
        }
    }
}

