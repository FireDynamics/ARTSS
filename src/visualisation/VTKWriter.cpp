/// \file       VTKWriter.cpp
/// \brief      class to write out vtk files
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "VTKWriter.h"
#include "../DomainData.h"
#include "visit_writer.h"  //( https://wci.llnl.gov/codes/visit/ )

static std::string ending = ".vtk";
int VTKWriter::vtk_counter = 0;

void VTKWriter::write_numerical_debug(const FieldController &field_controller, const std::string &filename) {
    auto u = field_controller.get_field_u_data();
    auto v = field_controller.get_field_v_data();
    auto w = field_controller.get_field_w_data();
    auto p = field_controller.get_field_p_data();
    auto div = field_controller.get_field_rhs_data();
    auto T = field_controller.get_field_T_data();
    auto C = field_controller.get_field_concentration_data();
    auto s = field_controller.get_field_sight_data();
    auto nu_t = field_controller.get_field_nu_t_data();
    auto S_T = field_controller.get_field_source_T_data();
    auto S_C = field_controller.get_field_source_concentration_data();
    auto f_x = field_controller.get_field_force_x_data();
    auto f_y = field_controller.get_field_force_y_data();
    auto f_z = field_controller.get_field_force_z_data();
    auto kappa = field_controller.get_field_kappa_data();
    auto gamma = field_controller.get_field_gamma_data();

    int size_vars = 16;
    read_ptr data[16] = {u, v, w, p, div, T, C, s, nu_t, S_T, S_C, f_x, f_y, f_z, kappa, gamma};
    // Initialise variables
    int var_dims[size_vars + 6]; // Dimensions of variables (x,y,z,u,v,w,p,div,T)
    int centering[size_vars + 6]; // Whether the variables are centered in a cell: 0 for zonal!
    for (int i = 0; i < size_vars + 6; i++) {
        var_dims[i] = 1;
        centering[i] = 0;
    }
    const char *var_names[] = {"x-coords", "y-coords", "z-coords",
                               "index_i", "index_j", "index_k",
                               "x-velocity", "y-velocity", "z-velocity",
                               "pressure",
                               "divergence",
                               "temperature",
                               "concentration",
                               "sight",
                               "nu_t",
                               "source_T", "source_C",
                               "force_x", "force_y", "force_z",
                               "kappa",
                               "gamma"};

    VTKWriter::vtk_prepare_and_write_debug((filename + ending).c_str(), data,
                                           size_vars, var_names, centering, var_dims);

}

void VTKWriter::write_numerical(const FieldController &field_controller, const std::string &filename) {
    return_ptr u = field_controller.get_field_u_data();
    return_ptr v = field_controller.get_field_v_data();
    return_ptr w = field_controller.get_field_w_data();
    return_ptr p = field_controller.get_field_p_data();
    return_ptr div = field_controller.get_field_rhs_data();
    return_ptr T = field_controller.get_field_T_data();
    return_ptr C = field_controller.get_field_concentration_data();
    return_ptr sight = field_controller.get_field_sight_data();
    return_ptr nu_t = field_controller.get_field_nu_t_data();
    return_ptr source_T = field_controller.get_field_source_T_data();
    VTKWriter::vtk_prepare_and_write((filename + ending).c_str(),
                                     u, v, w,
                                     p,
                                     div,
                                     T,
                                     C,
                                     sight,
                                     nu_t,
                                     source_T);
}

void VTKWriter::write_analytical(const Solution &solution, const std::string &filename) {
    auto u = solution.get_return_ptr_data_u();
    auto v = solution.get_return_ptr_data_v();
    auto w = solution.get_return_ptr_data_w();
    auto p = solution.get_return_ptr_data_p();
    auto T = solution.get_return_ptr_data_T();
    VTKWriter::vtk_prepare_and_write((filename + ending).c_str(), u, v, w, p, T);
}

//================================= Visualization (VTK) ==================================
// ***************************************************************************************
/// \brief  Prepares the (numerical) arrays in a correct format and writes the structured grid
///         and its variables
/// \param  fname xml-file name (via argument)
/// \param  u     constant input value (\a x -velocity)
/// \param  v     constant input value (\a y -velocity)
/// \param  w     constant input value (\a z -velocity)
/// \param  p     constant input value (pressure)
/// \param  div   constant input value (divergence)
/// \param  T     constant input value (temperature)
/// \param  C     constant input value (concentration)
/// \param  sight     constant input value (sight)
/// \param  nu_t    constant input value (turbulent viscosity)
/// \param  source_T   constant input values (energy source)
/// \author Severt
// ***************************************************************************************
void VTKWriter::vtk_prepare_and_write(const char *filename,
                                      read_ptr u, read_ptr v, read_ptr w,
                                      read_ptr p,
                                      read_ptr div,
                                      read_ptr T,
                                      read_ptr C,
                                      read_ptr sight,
                                      read_ptr nu_t,
                                      read_ptr source_T) {
    DomainData *domain = DomainData::getInstance();
    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    int Nx = static_cast<int>(domain->get_Nx());
    int Ny = static_cast<int>(domain->get_Ny());
    int Nz = static_cast<int>(domain->get_Nz());

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    int size = static_cast<int>(domain->get_size());

    // Initialise variables
    int size_vars = 13;  // Number of variables
    // Dimensions of variables (x,y,z,u,v,w,p,div,T,C,s,nu_t)
    int var_dims[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    // Whether the variables are centered in a cell: 0 for zonal!
    int centering[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const char *var_names[] = {"x-coords", "y-coords", "z-coords",
                               "x-velocity", "y-velocity", "z-velocity",
                               "pressure",
                               "divergence",
                               "temperature",
                               "concentration",
                               "sight",
                               "turb_visc",
                               "source_T"};

    // Dimensions of the rectilinear array (+1 for zonal values)
    int dims[] = {Nx + 1, Ny + 1, Nz + 1};

    auto *x_coords = new float[(Nx + 1)];
    auto *y_coords = new float[(Ny + 1)];
    auto *z_coords = new float[(Nz + 1)];

    // Initialise grid
    // faces of the grid cells
    for (int i = 0; i < Nx + 1; i++) {
        x_coords[i] = static_cast<float> (X1 + (i - 1) * dx);
    }

    for (int j = 0; j < Ny + 1; j++) {
        y_coords[j] = static_cast<float> (Y1 + (j - 1) * dy);
    }

    for (int k = 0; k < Nz + 1; k++) {
        z_coords[k] = static_cast<float> (Z1 + (k - 1) * dz);
    }

    // centers of the grid cells
    auto *x_centres = new float[size];
    auto *y_centres = new float[size];
    auto *z_centres = new float[size];

    // velocities
    auto *u_vel = new float[size];
    auto *v_vel = new float[size];
    auto *w_vel = new float[size];
    // pressure
    auto pres = new float[size];
    // divergence
    auto vel_div = new float[size];
    // temperature
    auto Temp = new float[size];
    // smoke concentration
    auto Con = new float[size];
    // boundary sight
    auto Sight = new float[size];
    // turbulent viscosity
    auto turb_visc = new float[size];
    // energy source
    auto Source_T = new float[size];

    // Cast variables to floats
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                x_centres[index] = x_coords[i] + static_cast<float> (0.5 * dx);
                y_centres[index] = y_coords[j] + static_cast<float> (0.5 * dy);
                z_centres[index] = z_coords[k] + static_cast<float> (0.5 * dz);
                u_vel[index] = static_cast<float>(u[index]);
                v_vel[index] = static_cast<float>(v[index]);
                w_vel[index] = static_cast<float>(w[index]);
                pres[index] = static_cast<float>(p[index]);
                vel_div[index] = static_cast<float>(div[index]);
                Temp[index] = static_cast<float>(T[index]);
                Con[index] = static_cast<float>(C[index]);
                Sight[index] = static_cast<float>(sight[index]);
                turb_visc[index] = static_cast<float>(nu_t[index]);
                Source_T[index] = static_cast<float>(source_T[index]);
            }
        }
    }

    // Summarise pointers to variables in an array
    float *vars[] = {static_cast<float *> (x_centres),
                     static_cast<float *> (y_centres),
                     static_cast<float *> (z_centres),
                     static_cast<float *> (u_vel),
                     static_cast<float *> (v_vel),
                     static_cast<float *> (w_vel),
                     static_cast<float *> (pres),
                     static_cast<float *> (vel_div),
                     static_cast<float *> (Temp),
                     static_cast<float *> (Con),
                     static_cast<float *> (Sight),
                     static_cast<float *> (turb_visc),
                     static_cast<float *> (Source_T)};

    // Use visit_writer to write data on mesh
    write_rectilinear_mesh(filename, 1, dims,
                           x_coords, y_coords, z_coords,
                           size_vars, var_dims, centering,
                           var_names, vars);

    // Clean up
    delete[] (x_coords);
    delete[] (y_coords);
    delete[] (z_coords);
    delete[] (x_centres);
    delete[] (y_centres);
    delete[] (z_centres);
    delete[] (u_vel);
    delete[] (v_vel);
    delete[] (w_vel);
    delete[] (pres);
    delete[] (vel_div);
    delete[] (Temp);
    delete[] (Con);
    delete[] (Sight);
    delete[] (turb_visc);
    delete[] (Source_T);
}

//================================= Visualization (VTK) ============================================
// *************************************************************************************************
/// \brief  Prepares the (analytical) arrays in a correct format and writes the structured grid
///         and its variables
/// \param  filename  xml filename (via argument)
/// \param  u     constant input value (\a x -velocity)
/// \param  v     constant input value (\a y -velocity)
/// \param  w     constant input value (\a z -velocity)
/// \param  p     constant input value (pressure)
/// \param  T     constant input value (temperature)
/// \param  s     constant input value (sight)
/// \author Severt
// *************************************************************************************************
void VTKWriter::vtk_prepare_and_write(const char *filename,
                                      read_ptr u, read_ptr v, read_ptr w,
                                      read_ptr p,
                                      read_ptr T) {
    DomainData *domain = DomainData::getInstance();
    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    int Nx = static_cast<int>(domain->get_Nx());
    int Ny = static_cast<int>(domain->get_Ny());
    int Nz = static_cast<int>(domain->get_Nz());

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    int size = static_cast<int>(domain->get_size());

    // Initialise variables
    int size_vars = 8; // Number of variables
    // Dimensions of variables (x,y,z,u,v,w,p,div,T)
    int var_dims[] = {1, 1, 1, 1, 1, 1, 1, 1};
    // Whether the variables are centered in a cell: 0 for zonal!
    int centering[] = {0, 0, 0, 0, 0, 0, 0, 0};
    const char *var_names[] = {"x-coords", "y-coords", "z-coords",
                               "x-velocity", "y-velocity", "z-velocity",
                               "pressure",
                               "temperature"};

    // Dimensions of the rectilinear array (+1 for zonal values)
    int dims[] = {Nx + 1, Ny + 1, Nz + 1};

    auto x_coords = new float[(Nx + 1)];
    auto y_coords = new float[(Ny + 1)];
    auto z_coords = new float[(Nz + 1)];

    // Initialise grid
    // faces of the grid cells
    for (int i = 0; i < Nx + 1; i++) {
        x_coords[i] = static_cast<float> (X1 + (i - 1) * dx);
    }

    for (int j = 0; j < Ny + 1; j++) {
        y_coords[j] = static_cast<float> (Y1 + (j - 1) * dy);
    }

    for (int k = 0; k < Nz + 1; k++) {
        z_coords[k] = static_cast<float> (Z1 + (k - 1) * dz);
    }

    // centers of the grid cells
    auto x_centres = new float[size];
    auto y_centres = new float[size];
    auto z_centres = new float[size];

    // velocities
    auto u_vel = new float[size];
    auto v_vel = new float[size];
    auto w_vel = new float[size];
    // pressure
    auto pres = new float[size];
    // temperature
    auto Temp = new float[size];

    // Cast variables to floats
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                x_centres[index] = x_coords[i] + static_cast<float> (0.5 * dx);
                y_centres[index] = y_coords[j] + static_cast<float> (0.5 * dy);
                z_centres[index] = z_coords[k] + static_cast<float> (0.5 * dz);
                u_vel[index] = static_cast<float>(u[index]);
                v_vel[index] = static_cast<float>(v[index]);
                w_vel[index] = static_cast<float>(w[index]);
                pres[index] = static_cast<float>(p[index]);
                Temp[index] = static_cast<float>(T[index]);
            }
        }
    }
    // Summarise pointers to variables in an array
    float *vars[] = {static_cast<float *> (x_centres),
                     static_cast<float *> (y_centres),
                     static_cast<float *> (z_centres),
                     static_cast<float *> (u_vel),
                     static_cast<float *> (v_vel),
                     static_cast<float *> (w_vel),
                     static_cast<float *> (pres),
                     static_cast<float *> (Temp)};

    // Use visit_writer to write data on mesh
    write_rectilinear_mesh(filename, 1, dims,
                           x_coords, y_coords, z_coords,
                           size_vars, var_dims, centering, var_names, vars);

    // Clean up
    delete[] (x_coords);
    delete[] (y_coords);
    delete[] (z_coords);
    delete[] (x_centres);
    delete[] (y_centres);
    delete[] (z_centres);
    delete[] (u_vel);
    delete[] (v_vel);
    delete[] (w_vel);
    delete[] (pres);
    delete[] (Temp);
}


//================================= Visualization (VTK) ============================================
// *************************************************************************************************
/// \brief  Prepares the arrays given in data in a correct format and writes the structured grid
/// and its variables
/// \param  filename  xml filename (via argument)
/// \param  data     containing all data fields to be written out
/// \param  size_vars length of data
/// \param  var_names title to be written out for the respective field in data
/// \param  centering  whether values are in the middle (0 for zonal value)
/// \param  var_dims  dimension of the respective field in data
/// \author Severt
// *************************************************************************************************
void VTKWriter::vtk_prepare_and_write_debug(const char *filename, read_ptr *data,
                                            int size_vars, const char * const *var_names,
                                            int *centering, int *var_dims) {
    DomainData *domain = DomainData::getInstance();
    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    int Nx = static_cast<int>(domain->get_Nx());
    int Ny = static_cast<int>(domain->get_Ny());
    int Nz = static_cast<int>(domain->get_Nz());

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    int size = static_cast<int>(domain->get_size());

    // Dimensions of the rectilinear array (+1 for zonal values)
    int dims[] = {Nx + 1, Ny + 1, Nz + 1};

    auto x_coords = new float[(Nx + 1)];
    auto y_coords = new float[(Ny + 1)];
    auto z_coords = new float[(Nz + 1)];

    // Initialise grid
    // faces of the grid cells
    for (int i = 0; i < Nx + 1; i++) {
        x_coords[i] = static_cast<float> (X1 + (i - 1) * dx);
    }

    for (int j = 0; j < Ny + 1; j++) {
        y_coords[j] = static_cast<float> (Y1 + (j - 1) * dy);
    }

    for (int k = 0; k < Nz + 1; k++) {
        z_coords[k] = static_cast<float> (Z1 + (k - 1) * dz);
    }

    float *write_out[size_vars + 6];
    for (int i = 0; i < size_vars + 6; i++){
        write_out[i] = new float[size];
    }

    // Cast variables to floats
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                write_out[0][index] = x_coords[i] + static_cast<float> (0.5 * dx);
                write_out[1][index] = y_coords[j] + static_cast<float> (0.5 * dy);
                write_out[2][index] = z_coords[k] + static_cast<float> (0.5 * dz);
                write_out[3][index] = static_cast<float> (i);
                write_out[4][index] = static_cast<float> (j);
                write_out[5][index] = static_cast<float> (k);

                for (int v = 0; v < size_vars; v++){
                    write_out[v + 6][index] = static_cast<float>(data[v][index]);
                }
            }
        }
    }
    // Summarise pointers to variables in an array
    float *vars[size_vars + 6];
    for (int i = 0; i < size_vars + 6; i++){
        vars[i] = static_cast<float *> (write_out[i]);
    }
    // Use visit_writer to write data on mesh
    write_rectilinear_mesh(filename, 1, dims, x_coords, y_coords, z_coords, size_vars + 6, var_dims, centering, var_names, vars);


    for (int i = 0; i < size_vars + 6; i++) {
        delete[] write_out[i];
    }
    // Clean up
    delete[] (x_coords);
    delete[] (y_coords);
    delete[] (z_coords);
}

void VTKWriter::write_field(const Field &field, const std::string &filename, const std::string &var_name) {
    /*
    std::string fname = std::to_string(VTKWriter::vtk_counter++) + "_" + filename;
    int size_vars = 1;
    const char *var_names[] = {"x-coords", "y-coords", "z-coords",
                               "index_i", "index_j", "index_k",
                               var_name.c_str()};
    DomainData *domain = DomainData::getInstance();
    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    int Nx = static_cast<int>(domain->get_Nx());
    int Ny = static_cast<int>(domain->get_Ny());
    int Nz = static_cast<int>(domain->get_Nz());

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    int size = static_cast<int>(domain->get_size());

    // Dimensions of the rectilinear array (+1 for zonal values)
    int dims[] = {Nx + 1, Ny + 1, Nz + 1};

    auto x_coords = new float[(Nx + 1)];
    auto y_coords = new float[(Ny + 1)];
    auto z_coords = new float[(Nz + 1)];

    // Initialise grid
    // faces of the grid cells
    for (int i = 0; i < Nx + 1; i++) {
        x_coords[i] = static_cast<float> (X1 + (i - 1) * dx);
    }

    for (int j = 0; j < Ny + 1; j++) {
        y_coords[j] = static_cast<float> (Y1 + (j - 1) * dy);
    }

    for (int k = 0; k < Nz + 1; k++) {
        z_coords[k] = static_cast<float> (Z1 + (k - 1) * dz);
    }

    float *write_out[size_vars + 6];
    for (int i = 0; i < size_vars + 6; i++){
        write_out[i] = new float[size];
    }

    // Cast variables to floats
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                write_out[0][index] = x_coords[i] + static_cast<float> (0.5 * dx);
                write_out[1][index] = y_coords[j] + static_cast<float> (0.5 * dy);
                write_out[2][index] = z_coords[k] + static_cast<float> (0.5 * dz);
                write_out[3][index] = static_cast<float> (i);
                write_out[4][index] = static_cast<float> (j);
                write_out[5][index] = static_cast<float> (k);
                write_out[6][index] = static_cast<float>(field.data[index]);
            }
        }
    }
    // Summarise pointers to variables in an array
    float *vars[size_vars + 6];
    for (int i = 0; i < size_vars + 6; i++){
        vars[i] = static_cast<float *> (write_out[i]);
    }
    int var_dims[] = {1, 1, 1, 1, 1, 1, 1};
    int centering[] = {0, 0, 0, 0, 0, 0, 0};
    write_rectilinear_mesh(fname.c_str(), 1, dims, x_coords, y_coords, z_coords, size_vars + 6, var_dims, centering, var_names, vars);


    for (int i = 0; i < size_vars + 6; i++) {
        delete[] write_out[i];
    }
    // Clean up
    delete[] (x_coords);
    delete[] (y_coords);
    delete[] (z_coords);
     */
}
