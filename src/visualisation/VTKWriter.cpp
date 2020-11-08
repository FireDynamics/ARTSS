/// \file       VTKWriter.cpp
/// \brief      class to write out vtk files
/// \date       Jun 25, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "VTKWriter.h"
#include "../Domain.h"
#include "visit_writer.h"  //( https://wci.llnl.gov/codes/visit/ )

#ifdef USEMPI
    #include "../utility/MPIHandler.h"
#endif

static std::string ending = ".vtk";

void VTKWriter::write_numerical(FieldController *field_controller, const std::string& filename) {
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
#ifdef USEMPI
    auto mpi_handler = MPIHandler::getInstance();
#endif
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
#ifdef USEMPI
    std::vector<int> rank_has_boundary {mpi_handler->get_mpi_neighbour()};
    int bx0{rank_has_boundary.at(Patch::LEFT)};
    int bx1{rank_has_boundary.at(Patch::RIGHT)};
    int by0{rank_has_boundary.at(Patch::BOTTOM)};
    int by1{rank_has_boundary.at(Patch::TOP)};
    int bz0{rank_has_boundary.at(Patch::BACK)};
    int bz1{rank_has_boundary.at(Patch::FRONT)};
#else
    int bz0{0};
    int bz1{0};
    int by0{0};
    int by1{0};
    int bx0{0};
    int bx1{0};
#endif

    int dims[] = {Nx + 1 - bx0 - bx1, Ny + 1 - by0 - by1, Nz + 1 - bz0 - bz1};            // Dimensions of the rectilinear array (+1 for zonal values)

    int size = static_cast<int>(dims[0]*dims[1]*dims[2]);

    auto *x_coords = new float[dims[0]];
    auto *y_coords = new float[dims[1]];
    auto *z_coords = new float[dims[2]];

    // Initialize grid
    // faces of the grid cells
    for (int i = bx0; i < Nx + 1 - bx1; i++) {
        x_coords[i-bx0] = static_cast<float> (X1 + (i - 1) * dx);
    }

    for (int j = by0; j < Ny + 1 - by1; j++) {
        y_coords[j-by0] = static_cast<float> (Y1 + (j - 1) * dy);
    }

    for (int k = bz0; k < Nz + 1 - bz1; k++) {
        z_coords[k-bz0] = static_cast<float> (Z1 + (k - 1) * dz);
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
    for (int k = bz0; k < Nz - bz1; k++) {
        for (int j = by0; j < Ny - by1; j++) {
            for (int i = bx0; i < Nx - bx1; i++) {
                size_t indexData = IX(i, j, k, Nx, Ny);
                size_t indexVTK = IX(i - bx0, j - by0, k - bz0, Nx - bx0 - bx1, Ny - by0 - by1);
                x_centres[indexVTK] = x_coords[i - bx0] + static_cast<float> (0.5 * dx);
                y_centres[indexVTK] = y_coords[j - by0] + static_cast<float> (0.5 * dy);
                z_centres[indexVTK] = z_coords[k - bz0] + static_cast<float> (0.5 * dz);
                u_vel[indexVTK] = static_cast<float>(u[indexData]);
                v_vel[indexVTK] = static_cast<float>(v[indexData]);
                w_vel[indexVTK] = static_cast<float>(w[indexData]);
                pres[indexVTK] = static_cast<float>(p[indexData]);
                vel_div[indexVTK] = static_cast<float>(div[indexData]);
                Temp[indexVTK] = static_cast<float>(T[indexData]);
                Con[indexVTK] = static_cast<float>(C[indexData]);
                Sight[indexVTK] = static_cast<float>(s[indexData]);
                turb_visc[indexVTK] = static_cast<float>(nu_t[indexData]);
                source_T[indexVTK] = static_cast<float>(S_T[indexData]);
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
#ifdef USEMPI
    auto mpi_handler = MPIHandler::getInstance();
#endif
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

    // Initialize variables
    int size_vars = 8; // Number of variables
    int var_dims[] = {1, 1, 1, 1, 1, 1, 1, 1}; // Dimensions of variables (x,y,z,u,v,w,p,div,T)
    int centering[] = {0, 0, 0, 0, 0, 0, 0, 0}; // Whether the variables are centered in a cell: 0 for zonal!
    const char *var_names[] = {"x-coords", "y-coords", "z-coords", "x-velocity", "y-velocity", "z-velocity", "pressure", "temperature"};

#ifdef USEMPI
    std::vector<int> rank_has_boundary {mpi_handler->get_mpi_neighbour()};
    int bx0{rank_has_boundary.at(Patch::LEFT)};
    int bx1{rank_has_boundary.at(Patch::RIGHT)};
    int by0{rank_has_boundary.at(Patch::BOTTOM)};
    int by1{rank_has_boundary.at(Patch::TOP)};
    int bz0{rank_has_boundary.at(Patch::BACK)};
    int bz1{rank_has_boundary.at(Patch::FRONT)};
#else
    int bz0{0};
    int bz1{0};
    int by0{0};
    int by1{0};
    int bx0{0};
    int bx1{0};
#endif

    int dims[] = {Nx + 1 - bx0 - bx1, Ny + 1 - by0 - by1, Nz + 1 - bz0 - bz1};            // Dimensions of the rectilinear array (+1 for zonal values)

    int size = static_cast<int>(dims[0]*dims[1]*dims[2]);

    auto x_coords = new float[dims[0]];
    auto y_coords = new float[dims[1]];
    auto z_coords = new float[dims[2]];

    // Initialize grid
    // faces of the grid cells
    for (int i = bx0; i < Nx + 1 - bx1; i++) {
        x_coords[i-bx0] = static_cast<float> (X1 + (i - 1) * dx);
    }

    for (int j = by0; j < Ny + 1 - by1; j++) {
        y_coords[j-by0] = static_cast<float> (Y1 + (j - 1) * dy);
    }

    for (int k = bz0; k < Nz + 1 - bz1; k++) {
        z_coords[k-bz0] = static_cast<float> (Z1 + (k - 1) * dz);
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
    for (int k = bz0; k < Nz - bz1; k++) {
        for (int j = by0; j < Ny - by1; j++) {
            for (int i = bx0; i < Nx - bx1; i++) {
                size_t indexData = IX(i, j, k, Nx, Ny);
                size_t indexVTK = IX(i - bx0, j - by0, k - bz0, Nx - bx0 - bx1, Ny - by0 - by1);
                x_centres[indexVTK] = x_coords[i - bx0] + static_cast<float> (0.5 * dx);
                y_centres[indexVTK] = y_coords[j - by0] + static_cast<float> (0.5 * dy);
                z_centres[indexVTK] = z_coords[k - bz0] + static_cast<float> (0.5 * dz);
                u_vel[indexVTK] = static_cast<float>(u[indexData]);
                v_vel[indexVTK] = static_cast<float>(v[indexData]);
                w_vel[indexVTK] = static_cast<float>(w[indexData]);
                pres[indexVTK] = static_cast<float>(p[indexData]);
                Temp[indexVTK] = static_cast<float>(T[indexData]);
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
