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

Visual::Visual() {
    auto params = Parameters::getInstance();
    std::string fname = params->get_filename();
    m_filename = RemoveExtension(fname);

    m_save_csv = (params->get("visualisation/save_csv") == "Yes");
    m_save_csv = (params->get("visualisation/save_vtk") == "Yes");

    m_has_analytical_solution = (params->get("solver/solution/available") == "Yes");
}

void Visual::visualise(ISolver *solver, const real t) {
    // TODO 0 fields for t == 0. (change ? remove time independent solutions)

    /*              solver->GetU0(), solver->GetV0(), solver->GetW0(), \
                    solver->GetP0(), solver->GetRhs(), solver->GetT0(), solver->GetC0(), \
                    solver->GetSight(), solver->GetNu_t(), solver->GetS_T(), \
                    m_solution.GetU0(), m_solution.GetV0(), m_solution.GetW0(), \
                    m_solution.GetP0(), m_solution.GetT0());
                    */
    if (m_has_analytical_solution){
         m_solution.CalcAnalyticalSolution(t);
    }
    if (m_save_vtk) {
        VTKWriter::write_numerical(solver, createFilename(m_filename, t, false));
        if (m_has_analytical_solution){
            VTKWriter::write_analytical(&m_solution, createFilename(m_filename, t, true));
        }
    }

    if (m_save_csv) {
        CSVWriter::write_numerical(solver, createFilename(m_filename, t, false));
        if (m_has_analytical_solution){
            CSVWriter::write_analytical(&m_solution, createFilename(m_filename, t, true));
        }
    }
}
void Visual::initialiseGrid(float *x_centres, float *y_centres, float *z_centres, int Nx, int Ny, int Nz, real dx, real dy, real dz) {
    Domain *domain = Domain::getInstance();
    real X1 = domain->GetX1();
    real Y1 = domain->GetY1();
    real Z1 = domain->GetZ1();

    // Initialize grid
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                x_centres[index] = static_cast<float> (X1 + (i-0.5) * dx);
                y_centres[index] = static_cast<float> (Y1 + (j-0.5) * dy);
                z_centres[index] = static_cast<float> (Z1 + (k-0.5) * dz);
            }
        }
    }
}

void Visual::prepareFields(read_ptr *fields, float **vars, int size){
    Domain *domain = Domain::getInstance();

    int Nx = static_cast<int>(domain->GetNx());
    int Ny = static_cast<int>(domain->GetNy());
    int Nz = static_cast<int>(domain->GetNz());

    // Cast variables to floats
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                for (int v = 0; v < size; v++){
                    vars[v][index] = static_cast<float>(fields[v][index]);
                }
            }
        }
    }
}

std::string Visual::createFilename(std::string filename, const real t, bool analytical){
    std::string fname = std::move(filename);
    if (analytical){
        fname.append("_ana_");
    }else {
        fname.append("_num_");
    }
    std::ostringstream tstep;
    tstep << std::setw(6) << std::setfill('0') << t;
    fname.append(tstep.str());
    return fname;
}

//================================= Remove extension ==================================
// ***************************************************************************************
/// \brief  Removes extension from filename
/// \param  filename    xml-file name (via argument)
// ***************************************************************************************
std::string Visual::RemoveExtension(const std::string &filename) {
    size_t lastdot = filename.find_last_of('.');
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot);
}
