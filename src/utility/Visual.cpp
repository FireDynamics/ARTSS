/// \file       Visual.cpp
/// \brief      Writes all the data needed for post-processing to csv and to a single binary file in VTK legacy format in binary
/// \details    Prepares the arrays in a correct format and writes the structured grid and its zonal variables (\a u,\a v, \a div, \a p,\a d,\a T)
/// \date       Jul 12, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

#include "Visual.h"
#include "visit_writer.h"   //( https://wci.llnl.gov/codes/visit/ )
#include "Parameters.h"
#include "../Domain.h"

Visual::Visual() {
    auto params = Parameters::getInstance();
    auto domain = Domain::getInstance();

    m_nx = static_cast<int>(domain->Getnx());
    m_ny = static_cast<int>(domain->Getny());
    m_nz = static_cast<int>(domain->Getnz());

    m_x1 = domain->Getx1();
    m_y1 = domain->Gety1();
    m_z1 = domain->Getz1();

    m_X1 = domain->GetX1();
    m_Y1 = domain->GetY1();
    m_Z1 = domain->GetZ1();

    m_dx = domain->Getdx();
    m_dy = domain->Getdy();
    m_dz = domain->Getdz();

    m_Nx = static_cast<int>(domain->GetNx());
    m_Ny = static_cast<int>(domain->GetNy());
    m_Nz = static_cast<int>(domain->GetNz());

    m_size = static_cast<int>(domain->GetSize());

    real dt     = params->getReal("physical_parameters/dt");
    real t_end  = params->getReal("physical_parameters/t_end");
    m_Nt        = int(t_end/dt);
}

// ==================================== Visualize ====================================
// ***************************************************************************************
/// \brief  Save data in vtk file
/// \param  solver      solver
/// \param  t           time
/// \param  fname       xml-file name (via argument)
// ***************************************************************************************
void Visual::Visualize(SolverI* solver, const real t, const char *fname){

    auto params = Parameters::getInstance();

    real dt = params->getReal("physical_parameters/dt");
    int n   = static_cast<int>(std::round(t/dt));

    if (t == 0.){
        if (params->get("solver/solution/available")=="Yes") m_solution.CalcAnalyticalSolution(0.);
        vtkWriteStep(fname, n,  solver->GetU0(), solver->GetV0(), solver->GetW0(), \
                    solver->GetP0(), solver->GetRhs(), solver->GetT0(), solver->GetC0(), \
                    solver->GetSight(), solver->GetNu_t(), solver->GetS_T(),\
                    m_solution.GetU0(), m_solution.GetV0(), m_solution.GetW0(), \
                    m_solution.GetP0(), m_solution.GetT0());
    }
    else {
        if (params->get("solver/solution/available")=="Yes") m_solution.CalcAnalyticalSolution(t);
        vtkWriteStep(fname, n,  solver->GetU(), solver->GetV(), solver->GetW(), \
                    solver->GetP(), solver->GetRhs(), solver->GetT(), solver->GetC(), \
                    solver->GetSight(),solver->GetNu_t(), solver->GetS_T(),\
                    m_solution.GetU(), m_solution.GetV(), m_solution.GetW(), \
                    m_solution.GetP(), m_solution.GetT());
    }
}

//================================= Write to vtk and csv ==================================
// ***************************************************************************************
/// \brief  Writes the structured grid and its variables at time step \a n
/// \param  fname   xml-file name (via argument)
/// \param  n       time step
/// \param  u       constant input value (\a x -velocity)
/// \param  v       constant input value (\a y -velocity)
/// \param  w       constant input value (\a z -velocity)
/// \param  p       constant input value (pressure)
/// \param  div     constant input value (divergence)
/// \param  T       constant input value (temperature)
/// \param  C       constant input value (concentration)
/// \param  s       constant input value (sight)
/// \param  nu_t    constant input value (turbulent viscosity)
/// \param  S_T     constant input value (energy source)
/// \param  ua      constant input value (analytical \a x -velocity)
/// \param  va      constant input value (analytical \a y -velocity)
/// \param  wa      constant input value (analytical \a z -velocity)
/// \param  pa      constant input value (analytical pressure)
/// \param  Ta      constant input value (analytical temperature)
// ***************************************************************************************
void Visual::vtkWriteStep(const char *fname, const int n, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr div, read_ptr T, read_ptr C, read_ptr s, read_ptr nu_t, read_ptr S_T, read_ptr ua, read_ptr va, read_ptr wa, read_ptr pa, read_ptr Ta){

    auto params = Parameters::getInstance();

    int n_plots     = params->getInt("visualization/n_plots");

    bool save_csv = false;
    if (params->get("visualization/save_csv")=="Yes") save_csv = true;

    // Visualization in time step n
    if (n == 0 || fmod(n, float(m_Nt/n_plots)) == 0 || n == m_Nt)
    {
    std::string filename = RemoveExtension(fname);
        std::string csv = RemoveExtension(fname);

        // numerical solution
        filename.append("_num_");
        csv.append("_num_");
        std::ostringstream tstep;
        tstep << std::setw( 6 ) << std::setfill( '0' ) << n;
        filename.append(tstep.str());
        csv.append(tstep.str());
        filename.append(".vtk");
        csv.append(".csv");
        const char * filen = filename.c_str();
        vtkPrepareAndWrite(filen, u, v, w, p, div, T, C, s, nu_t, S_T);
        if(save_csv) csvPrepareAndWrite(csv, u, v, w, p, div, T, C, s, nu_t, S_T);


        // analytical solution
        if (params->get("solver/solution/available")=="Yes"){
      std::string filename_a = RemoveExtension(fname);
            std::string csv_a = RemoveExtension(fname);
            filename_a = filename_a.append("_ana_");
            csv_a.append("_ana_");
            std::ostringstream tstep_a;
            tstep_a << std::setw( 6 ) << std::setfill( '0' ) << n;
            filename_a.append(tstep_a.str());
            csv_a.append(tstep_a.str());
            filename_a.append(".vtk");
            csv_a.append(".csv");
            const char * filen_a = filename_a.c_str();
            vtkPrepareAndWrite(filen_a, ua, va, wa, pa, Ta, s);
            if(save_csv) csvPrepareAndWrite(csv_a, ua, va, wa, pa, Ta, s);
        }
    }
}

//================================= Remove extension ==================================
// ***************************************************************************************
/// \brief  Removes extension from filename
/// \param  filename        xml-file name (via argument)
// ***************************************************************************************
std::string Visual::RemoveExtension(const std::string& filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}

//================================= Visualization (VTK) ==================================
// ***************************************************************************************
/// \brief  Prepares the (numerical) arrays in a correct format and writes the structured grid and its variables
/// \param  fname   xml-file name (via argument)
/// \param  u       constant input value (\a x -velocity)
/// \param  v       constant input value (\a y -velocity)
/// \param  w       constant input value (\a z -velocity)
/// \param  p       constant input value (pressure)
/// \param  div     constant input value (divergence)
/// \param  T       constant input value (temperature)
/// \param  C       constant input value (concentration)
/// \param  s       constant input value (sight)
/// \param  nu_t    constant input value (turbulent viscosity)
/// \param  S_T     constant input values (energy source)
// ***************************************************************************************
void Visual::vtkPrepareAndWrite(const char *fname, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr div, read_ptr T, read_ptr C, read_ptr s, read_ptr nu_t, read_ptr S_T) {

// Initialize variables
    int nvars       = 13; // Number of variables
    int vardims[]   = {1,1,1,1,1,1,1,1,1,1,1,1,1};      // Dimensions of variables (x,y,z,u,v,w,p,div,T,C,s,nu_t)
    int centering[] = {0,0,0,0,0,0,0,0,0,0,0,0,0};  // Whether the variables are centered in a cell: 0 for zonal!
    const char *varnames[] = {  "x-coords", "y-coords", "z-coords", "x-velocity", "y-velocity", "z-velocity", "pressure", "divergence",\
                                "temperature", "concentration", "sight", "turb_visc", "source_T"};

    int dims[] = {m_Nx+1, m_Ny+1, m_Nz+1};          // Dimensions of the rectilinear array (+1 for zonal values)

    float* xcoords = new float[(m_Nx+1)];
    float* ycoords = new float[(m_Ny+1)];
    float* zcoords = new float[(m_Nz+1)];

// Initialize grid
    // faces of the grid cells
    for(int i=0;i<m_Nx+1;i++){
         xcoords[i] = (float) (m_X1 + (i-1)*m_dx);
    }

    for(int j=0;j<m_Ny+1;j++){
         ycoords[j] = (float) (m_Y1 + (j-1)*m_dy);
    }

    for(int k=0;k<m_Nz+1;k++)
         zcoords[k] = (float) (m_Z1 + (k-1)*m_dz);

    float* xcenters = new float[m_size];
    float* ycenters = new float[m_size];
    float* zcenters = new float[m_size];

    // centers of the grid cells
    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                xcenters[IX(i,j,k,m_Nx,m_Ny)] = xcoords[i] + (float) (0.5*m_dx);
          }
      }
    }

    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                ycenters[IX(i,j,k,m_Nx,m_Ny)] = ycoords[j] + (float) (0.5*m_dy);
          }
      }
    }

    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                zcenters[IX(i,j,k,m_Nx,m_Ny)] = zcoords[k] + (float) (0.5*m_dz);
          }
      }
    }

// Cast variables to floats
    // velocities
    float* u_vel = new float[m_size];
    float* v_vel = new float[m_size];
    float* w_vel = new float[m_size];

    for (int i=0;i<m_size; i++){
        u_vel[i] = float (u[i]);
        v_vel[i] = float (v[i]);
        w_vel[i] = float (w[i]);
    }

    // pressure
    float* pres = new float[m_size];

    for(int i=0;i<m_size;i++)
        pres[i] = float (p[i]);

    // divergence
    float* veldiv = new float[m_size];

    for(int i=0;i<m_size;i++)
        veldiv[i] = float (div[i]);

    // temperature
    float* Temp = new float[m_size];

    for(int i=0;i<m_size;i++)
        Temp[i] = float (T[i]);

    // smoke concentration
    float* Con = new float[m_size];

    for(int i=0;i<m_size;i++)
        Con[i] = float (C[i]);

    // boundary sight
    float* Sight = new float[m_size];

    for(int i=0;i<m_size;i++)
        Sight[i] = float (s[i]);

    // turbulent viscosity
    float* turb_visc = new float[m_size];

    for(int i=0;i<m_size;i++)
        turb_visc[i] = float (nu_t[i]);

    // energy source
    float* source_T = new float[m_size];

    for(int i=0;i<m_size;i++)
        source_T[i] = float (S_T[i]);

// Summarize pointers to variables in an array
    float *vars[] = {   (float *)xcenters, (float *)ycenters, (float *)zcenters, (float *)u_vel, (float*)v_vel, (float*)w_vel, \
                        (float*)pres, (float*) veldiv, (float*) Temp, (float*) Con, (float*) Sight, (float*) turb_visc, (float*) source_T};

// Use visit_writer to write data on mesh
    write_rectilinear_mesh(fname, 1, dims, xcoords, ycoords, zcoords, nvars, vardims, centering, varnames, vars);

// Clean up
    delete[] (xcoords);
    delete[] (ycoords);
    delete[] (zcoords);
    delete[] (xcenters);
    delete[] (ycenters);
    delete[] (zcenters);
    delete[] (u_vel);
    delete[] (v_vel);
    delete[] (w_vel);
    delete[] (pres);
    delete[] (veldiv);
    delete[] (Temp);
    delete[] (Con);
    delete[] (Sight);
    delete[] (turb_visc);
    delete[] (source_T);
}

//================================= Visualization (VTK) ==================================
// ***************************************************************************************
/// \brief  Prepares the (analytical) arrays in a correct format and writes the structured grid and its variables
/// \param  fname   xml-file name (via argument)
/// \param  u       constant input value (\a x -velocity)
/// \param  v       constant input value (\a y -velocity)
/// \param  w       constant input value (\a z -velocity)
/// \param  p       constant input value (pressure)
/// \param  T       constant input value (temperature)
/// \param  s       constant input value (sight)
// ***************************************************************************************
void Visual::vtkPrepareAndWrite(const char *fname, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr T, read_ptr s) {

// Initialize variables
    int nvars       = 9;                        // Number of variables
    int vardims[]   = {1,1,1,1,1,1,1,1,1};      // Dimensions of variables (x,y,z,u,v,w,p,div,T,s)
    int centering[] = {0,0,0,0,0,0,0,0,0};      // Whether the variables are centered in a cell: 0 for zonal!
    const char *varnames[] = {"x-coords", "y-coords", "z-coords", "x-velocity", "y-velocity", "z-velocity", "pressure", "temperature", "sight"};

    int dims[] = {m_Nx+1, m_Ny+1, m_Nz+1}; // Dimensions of the rectilinear array (+1 for zonal values)

    float* xcoords = new float[(m_Nx+1)];
    float* ycoords = new float[(m_Ny+1)];
    float* zcoords = new float[(m_Nz+1)];

// Initialize grid
    // faces of the grid cells
    for(int i=0;i<m_Nx+1;i++){
         xcoords[i] = (float) (m_X1 + (i-1)*m_dx);
    }

    for(int j=0;j<m_Ny+1;j++){
         ycoords[j] = (float) (m_Y1 + (j-1)*m_dy);
    }

    for(int k=0;k<m_Nz+1;k++)
         zcoords[k] = (float) (m_Z1 + (k-1)*m_dz);

    float* xcenters = new float[m_size];
    float* ycenters = new float[m_size];
    float* zcenters = new float[m_size];

    // centers of the grid cells
    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                xcenters[IX(i,j,k,m_Nx,m_Ny)] = xcoords[i] + (float) (0.5*m_dx);
          }
      }
    }

    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                ycenters[IX(i,j,k,m_Nx,m_Ny)] = ycoords[j] + (float) (0.5*m_dy);
          }
      }
    }

    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                zcenters[IX(i,j,k,m_Nx,m_Ny)] = zcoords[k] + (float) (0.5*m_dz);
          }
      }
    }

// Cast variables to floats
    // velocities
    float* u_vel = new float[m_size];
    float* v_vel = new float[m_size];
    float* w_vel = new float[m_size];

    for (int i=0;i<m_size; i++){
        u_vel[i] = float (u[i]);
        v_vel[i] = float (v[i]);
        w_vel[i] = float (w[i]);
    }

    // pressure
    float* pres = new float[m_size];

    for(int i=0;i<m_size;i++)
        pres[i] = float (p[i]);

    // temperature
    float* Temp = new float[m_size];

    for(int i=0;i<m_size;i++)
        Temp[i] = float (T[i]);

    // boundary sight
    float* Sight = new float[m_size];

    for(int i=0;i<m_size;i++)
        Sight[i] = float (s[i]);

// Summarize pointers to variables in an array
    float *vars[] = {   (float *)xcenters, (float *)ycenters, (float *)zcenters, (float *)u_vel, (float*)v_vel, (float*)w_vel, \
                        (float*)pres, (float*) Temp, (float*) Sight};

// Use visit_writer to write data on mesh
    write_rectilinear_mesh(fname, 1, dims, xcoords, ycoords, zcoords, nvars, vardims, centering, varnames, vars);

// Clean up
    delete[] (xcoords);
    delete[] (ycoords);
    delete[] (zcoords);
    delete[] (xcenters);
    delete[] (ycenters);
    delete[] (zcenters);
    delete[] (u_vel);
    delete[] (v_vel);
    delete[] (w_vel);
    delete[] (pres);
    delete[] (Temp);
    delete[] (Sight);
}

//================================= Saving variables (.csv) ==================================
// ***************************************************************************************
/// \brief  This function is to prepare the arrays in a correct format and writes the structured grid and its variables to csv
/// \param  fname       name of csv file
/// \param  u           constant input value (\a x -velocity)
/// \param  v           constant input value (\a y -velocity)
/// \param  w           constant input value (\a z -velocity)
/// \param  p           constant input value (pressure)
/// \param  div         constant input value (divergence)
/// \param  T           constant input value (temperature)
/// \param  C           constant input value (concentration)
/// \param  s           constant input value (sight)
/// \param  nu_t        constant input value (turbulent viscosity)
/// \param  S_T         constant input value (temperature source)
// ***************************************************************************************
void Visual::csvPrepareAndWrite(std::string fname, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr div, read_ptr T, read_ptr C, read_ptr s, read_ptr nu_t, read_ptr S_T) {

// Initialize variables
    const char *varnames[] = {  "i", "j", "k", "idx", "x-coords (m)", "y-coords (m)", "z-coords (m)", \
                                "x-velocity (m/s)", "y-velocity (m/s)", "z-velocity (m/s)", "pressure (kg/(m s^2))", \
                                "divergence (1/s)", "temperature (Celsius)", "concentration (g/m^3)", "sight", \
                                "turb viscosity (m^2/s)", "temperature source (K/s)"};

    int nvars = 17;                 // Number of variables

    float* xcoords = new float[(m_Nx+1)];
    float* ycoords = new float[(m_Ny+1)];
    float* zcoords = new float[(m_Nz+1)];

// Initialize grid
    // faces of the grid cells
    for(int i=0;i<m_Nx+1;i++){
        xcoords[i] = (float) (m_X1 + (i-1)*m_dx);
    }

    for(int j=0;j<m_Ny+1;j++){
        ycoords[j] = (float) (m_Y1 + (j-1)*m_dy);
    }

    for(int k=0;k<m_Nz+1;k++)
        zcoords[k] = (float) (m_Z1 + (k-1)*m_dz);

    float* xcenters = new float[m_size];
    float* ycenters = new float[m_size];
    float* zcenters = new float[m_size];

    // centers of the grid cells
    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                xcenters[IX(i,j,k,m_Nx,m_Ny)] = xcoords[i] + (float) (0.5*m_dx);
          }
      }
    }

    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                ycenters[IX(i,j,k,m_Nx,m_Ny)] = ycoords[j] + (float) (0.5*m_dy);
          }
      }
    }

    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                zcenters[IX(i,j,k,m_Nx,m_Ny)] = zcoords[k] + (float) (0.5*m_dz);
          }
      }
    }

// Cast variables to floats
    // velocities
    float* u_vel = new float[m_size];
    float* v_vel = new float[m_size];
    float* w_vel = new float[m_size];

    for (int i=0;i<m_size; i++){
        u_vel[i] = float (u[i]);
        v_vel[i] = float (v[i]);
        w_vel[i] = float (w[i]);
    }

    // pressure
    float* pres = new float[m_size];

    for(int i=0;i<m_size;i++)
        pres[i] = float (p[i]);

    // divergence
    float* veldiv = new float[m_size];

    for(int i=0;i<m_size;i++)
        veldiv[i] = float (div[i]);

    // temperature
    float* Temp = new float[m_size];

    for(int i=0;i<m_size;i++)
        Temp[i] = float (T[i]);

    // smoke concentration
    float* Con = new float[m_size];

    for(int i=0;i<m_size;i++)
        Con[i] = float (C[i]);

    // boundary sight
    float* Sight = new float[m_size];

    for(int i=0;i<m_size;i++)
        Sight[i] = float (s[i]);

    // turbulent viscosity
    float* turb_visc = new float[m_size];

    for(int i=0;i<m_size;i++)
        turb_visc[i] = float (nu_t[i]);

    // temperature source
    float* temp_source = new float[m_size];

    for(int i=0;i<m_size;i++)
        temp_source[i] = float (S_T[i]);

    // Summarize pointers to variables in an array
    float *vars[] = {   (float *)xcenters, (float *)ycenters, (float *)zcenters, (float *)u_vel, (float*)v_vel, (float*)w_vel, \
                        (float*)veldiv, (float*)pres, (float*)Temp, (float*)Con, (float*)Sight, (float*)turb_visc, (float*)temp_source};

// Write data to csv
    const char delimiter = ',';
  std::ofstream outputFile;
    outputFile.open(fname, std::ofstream::out);

    // varnames as column titles
    for(int i=0;i<nvars-1;i++){
        outputFile<<varnames[i]<<delimiter;
    }
    // last column
    outputFile<<varnames[nvars-1]<<"\n";

    // write variables to csv
    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                outputFile<<i<<delimiter<<j<<delimiter<<k<<delimiter<<IX(i,j,k,m_Nx,m_Ny)<<delimiter\
                       <<std::setprecision(16)<<vars[0][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[1][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[2][IX(i,j,k,m_Nx,m_Ny)]\
                       <<delimiter<<vars[3][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[4][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[5][IX(i,j,k,m_Nx,m_Ny)]\
                       <<delimiter<<vars[6][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[7][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[8][IX(i,j,k,m_Nx,m_Ny)]\
                       <<delimiter<<vars[9][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[10][IX(i,j,k,m_Nx,m_Ny)]\
                       <<delimiter<<vars[11][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[12][IX(i,j,k,m_Nx,m_Ny)]<<"\n";
          }
      }
    }

// Clean up
    delete[] (xcoords);
    delete[] (ycoords);
    delete[] (zcoords);
    delete[] (xcenters);
    delete[] (ycenters);
    delete[] (zcenters);
    delete[] (u_vel);
    delete[] (v_vel);
    delete[] (w_vel);
    delete[] (pres);
    delete[] (veldiv);
    delete[] (Temp);
    delete[] (Con);
    delete[] (Sight);
    delete[] (turb_visc);
    delete[] (temp_source);

    outputFile.close();
}

void Visual::csvPrepareAndWrite(std::string fname, read_ptr u, read_ptr v, read_ptr w, read_ptr p, read_ptr T, read_ptr s) {

// Initialize variables
    const char *varnames[] = {  "i", "j", "k", "idx", "x-coords (m)", "y-coords (m)", "z-coords (m)", \
                                "x-velocity (m/s)", "y-velocity (m/s)", "z-velocity (m/s)", "pressure (kg/(m s^2))", \
                                "temperature (Celsius)", "sight"};

    int nvars = 13; // Number of variables

    float* xcoords = new float[(m_Nx+1)];
    float* ycoords = new float[(m_Ny+1)];
    float* zcoords = new float[(m_Nz+1)];

// Initialize grid
    // faces of the grid cells
    for(int i=0;i<m_Nx+1;i++){
        xcoords[i] = (float) (m_X1 + (i-1)*m_dx);
    }

    for(int j=0;j<m_Ny+1;j++){
        ycoords[j] = (float) (m_Y1 + (j-1)*m_dy);
    }

    for(int k=0;k<m_Nz+1;k++)
        zcoords[k] = (float) (m_Z1 + (k-1)*m_dz);

    float* xcenters = new float[m_size];
    float* ycenters = new float[m_size];
    float* zcenters = new float[m_size];

    // centers of the grid cells
    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                xcenters[IX(i,j,k,m_Nx,m_Ny)] = xcoords[i] + (float) (0.5*m_dx);
          }
      }
    }

    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                ycenters[IX(i,j,k,m_Nx,m_Ny)] = ycoords[j] + (float) (0.5*m_dy);
          }
      }
    }

    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                zcenters[IX(i,j,k,m_Nx,m_Ny)] = zcoords[k] + (float) (0.5*m_dz);
          }
      }
    }

// Cast variables to floats
    // velocities
    float* u_vel = new float[m_size];
    float* v_vel = new float[m_size];
    float* w_vel = new float[m_size];

    for (int i=0;i<m_size; i++){
        u_vel[i] = float (u[i]);
        v_vel[i] = float (v[i]);
        w_vel[i] = float (w[i]);
    }

    // pressure
    float* pres = new float[m_size];

    for(int i=0;i<m_size;i++)
        pres[i] = float (p[i]);

    // temperature
    float* Temp = new float[m_size];

    for(int i=0;i<m_size;i++)
        Temp[i] = float (T[i]);

    // boundary sight
    float* Sight = new float[m_size];

    for(int i=0;i<m_size;i++)
        Sight[i] = float (s[i]);

// Summarize pointers to variables in an array
    float *vars[] = {   (float *)xcenters, (float *)ycenters, (float *)zcenters, (float *)u_vel, (float*)v_vel, (float*)w_vel, \
                        (float*)pres, (float*) Temp, (float*) Sight};

// Write data to csv
    const char delimiter = ',';
  std::ofstream outputFile;
    outputFile.open(fname, std::ofstream::out);

    // varnames as column titles
    for(int i=0;i<nvars-1;i++){
      outputFile<<varnames[i]<<delimiter;
    }
    // last column
    outputFile<<varnames[nvars-1]<<"\n";

    // write variables to csv
    for (int k=0; k<m_Nz; k++) {\
        for (int j=0; j<m_Ny; j++) {\
            for (int i=0; i<m_Nx; i++) {
                outputFile<<i<<delimiter<<j<<delimiter<<k<<delimiter<<IX(i,j,k,m_Nx,m_Ny)<<delimiter\
                       <<std::setprecision(16)<<vars[0][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[1][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[2][IX(i,j,k,m_Nx,m_Ny)]\
                       <<delimiter<<vars[3][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[4][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[5][IX(i,j,k,m_Nx,m_Ny)]\
                       <<delimiter<<vars[6][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[7][IX(i,j,k,m_Nx,m_Ny)]<<delimiter<<vars[8][IX(i,j,k,m_Nx,m_Ny)]<<"\n";
          }
      }
    }

// Clean up
    delete[] (xcoords);
    delete[] (ycoords);
    delete[] (zcoords);
    delete[] (xcenters);
    delete[] (ycenters);
    delete[] (zcenters);
    delete[] (u_vel);
    delete[] (v_vel);
    delete[] (w_vel);
    delete[] (pres);
    delete[] (Temp);
    delete[] (Sight);

    outputFile.close();
}
