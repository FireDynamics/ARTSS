/// \file       VTKWriter.cpp
/// \brief      class to write out vtk files
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "VTKWriter.h"
#include "../Domain.h"
#include "visit_writer.h"  //( https://wci.llnl.gov/codes/visit/ )

static std::string ending = ".vtk";

void VTKWriter::write_numerical(ISolver *solver, const std::string& filename) {
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
    VTKWriter::vtkPrepareAndWrite((filename + ending).c_str(), u, v, w, p, div, T, C, s, nu_t, S_T);
}

void VTKWriter::write_analytical(Solution *solution, const std::string& filename) {
    auto u = solution->GetU();
    auto v = solution->GetV();
    auto w = solution->GetW();
    auto p = solution->GetP();
    auto T = solution->GetT();
    VTKWriter::vtkPrepareAndWrite((filename + ending).c_str(), u, v, w, p, T);
}

//================================= Visualization (VTK) ==================================
// ***************************************************************************************
/// \brief  Prepares the (numerical) arrays in a correct format and writes the structured grid and its variables
/// \param  fname xml-file name (via argument)
/// \param  u     constant input value (\a x -velocity)
/// \param  v     constant input value (\a y -velocity)
/// \param  w     constant input value (\a z -velocity)
/// \param  p     constant input value (pressure)
/// \param  div   constant input value (divergence)
/// \param  T     constant input value (temperature)
/// \param  C     constant input value (concentration)
/// \param  s     constant input value (sight)
/// \param  nu_t    constant input value (turbulent viscosity)
/// \param  S_T   constant input values (energy source)
/// \author Severt
// ***************************************************************************************
void VTKWriter::vtkPrepareAndWrite(const char *filename, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr div, read_ptr T, read_ptr C, read_ptr s, read_ptr nu_t, read_ptr S_T) {
    Domain *domain = Domain::getInstance();
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

// Initialize variables
    int size_vars = 13; // Number of variables
    int var_dims[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // Dimensions of variables (x,y,z,u,v,w,p,div,T,C,s,nu_t)
    int centering[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // Whether the variables are centered in a cell: 0 for zonal!
    const char *var_names[] = {"x-coords", "y-coords", "z-coords",
                              "x-velocity", "y-velocity", "z-velocity",
                              "pressure",
                              "divergence",
                              "temperature",
                              "concentration",
                              "sight",
                              "turb_visc",
                              "source_T"};

    int dims[] = {Nx + 1, Ny + 1, Nz + 1};            // Dimensions of the rectilinear array (+1 for zonal values)

    auto *x_coords = new float[(Nx + 1)];
    auto *y_coords = new float[(Ny + 1)];
    auto *z_coords = new float[(Nz + 1)];

    // Initialize grid
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
    auto source_T = new float[size];

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
                Sight[index] = static_cast<float>(s[index]);
                turb_visc[index] = static_cast<float>(nu_t[index]);
                source_T[index] = static_cast<float>(S_T[index]);
            }
        }
    }

    // Summarize pointers to variables in an array
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
                     static_cast<float *> (source_T)};

    // Use visit_writer to write data on mesh
    write_rectilinear_mesh(filename, 1, dims, x_coords, y_coords, z_coords, size_vars, var_dims, centering, var_names, vars);

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
    delete[] (source_T);
}

//================================= Visualization (VTK) ==================================
// ***************************************************************************************
/// \brief  Prepares the (analytical) arrays in a correct format and writes the structured grid and its variables
/// \param  filename  xml filename (via argument)
/// \param  u     constant input value (\a x -velocity)
/// \param  v     constant input value (\a y -velocity)
/// \param  w     constant input value (\a z -velocity)
/// \param  p     constant input value (pressure)
/// \param  T     constant input value (temperature)
/// \param  s     constant input value (sight)
/// \author Severt
// ***************************************************************************************
void VTKWriter::vtkPrepareAndWrite(const char *filename, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr T) {
    Domain *domain = Domain::getInstance();
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

    // Initialize variables
    int size_vars = 8; // Number of variables
    int var_dims[] = {1, 1, 1, 1, 1, 1, 1, 1}; // Dimensions of variables (x,y,z,u,v,w,p,div,T)
    int centering[] = {0, 0, 0, 0, 0, 0, 0, 0}; // Whether the variables are centered in a cell: 0 for zonal!
    const char *var_names[] = {"x-coords", "y-coords", "z-coords", "x-velocity", "y-velocity", "z-velocity", "pressure", "temperature"};

    int dims[] = {Nx + 1, Ny + 1, Nz + 1}; // Dimensions of the rectilinear array (+1 for zonal values)

    auto x_coords = new float[(Nx + 1)];
    auto y_coords = new float[(Ny + 1)];
    auto z_coords = new float[(Nz + 1)];

// Initialize grid
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
                z_centres[index] = z_coords[j] + static_cast<float> (0.5 * dz);
                u_vel[index] = static_cast<float>(u[index]);
                v_vel[index] = static_cast<float>(v[index]);
                w_vel[index] = static_cast<float>(w[index]);
                pres[index] = static_cast<float>(p[index]);
                Temp[index] = static_cast<float>(T[index]);
            }
        }
    }
    // Summarize pointers to variables in an array
    float *vars[] = {static_cast<float *> (x_centres),
                     static_cast<float *> (y_centres),
                     static_cast<float *> (z_centres),
                     static_cast<float *> (u_vel),
                     static_cast<float *> (v_vel),
                     static_cast<float *> (w_vel),
                     static_cast<float *> (pres),
                     static_cast<float *> (Temp)};

    // Use visit_writer to write data on mesh
    write_rectilinear_mesh(filename, 1, dims, x_coords, y_coords, z_coords, size_vars, var_dims, centering, var_names, vars);

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
