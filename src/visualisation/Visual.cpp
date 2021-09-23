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

Visual::Visual(const Solution &solution) : m_solution(solution) {
    auto params = Parameters::getInstance();
    m_filename = Utility::remove_extension(params->get_filename());

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
}

void Visual::visualise(const FieldController &field_controller, real t) {
    int n = static_cast<int> (std::round(t / m_dt));
    std::string filename = create_filename(m_filename, n, false);
    if (m_save_vtk) {
        if (fmod(n, m_vtk_plots) == 0 || t >= m_t_end) {
            VTKWriter::write_numerical(field_controller, filename);
            if (m_has_analytical_solution) {
                VTKWriter::write_analytical(m_solution, filename);
            }
        }
    }

    if (m_save_csv) {
        if (fmod(n, m_csv_plots) == 0 || t >= m_t_end) {
            CSVWriter::write_numerical(field_controller, filename);
            if (m_has_analytical_solution) {
                CSVWriter::write_analytical(m_solution, filename);
            }
        }
    }
}

void Visual::write_csv(const FieldController &field_controller, const std::string& filename){
    // local variables and parameters for GPU
    field_controller.field_u->update_host();
    field_controller.field_v->update_host();
    field_controller.field_w->update_host();
    field_controller.field_p->update_host();
    field_controller.field_rhs->update_host();
    field_controller.field_T->update_host();
    field_controller.field_concentration->update_host();
    field_controller.field_source_T->update_host();
    field_controller.field_source_concentration->update_host();
    field_controller.field_nu_t->update_host();
    CSVWriter::write_numerical(field_controller, filename);
}

void Visual::write_vtk(const FieldController &field_controller, const std::string& filename){
    field_controller.field_u->update_host();
    field_controller.field_v->update_host();
    field_controller.field_w->update_host();
    field_controller.field_p->update_host();
    field_controller.field_rhs->update_host();
    field_controller.field_T->update_host();
    field_controller.field_concentration->update_host();
    field_controller.field_source_T->update_host();
    field_controller.field_source_concentration->update_host();
    field_controller.field_nu_t->update_host();
    VTKWriter::write_numerical(field_controller, filename);
}

void Visual::write_vtk_debug(const FieldController &field_controller, const std::string& filename){
    field_controller.field_u->update_host();
    field_controller.field_v->update_host();
    field_controller.field_w->update_host();
    field_controller.field_p->update_host();
    field_controller.field_rhs->update_host();
    field_controller.field_T->update_host();
    field_controller.field_concentration->update_host();
    field_controller.field_source_T->update_host();
    field_controller.field_source_concentration->update_host();
    field_controller.field_nu_t->update_host();
    field_controller.field_force_x->update_host();
    field_controller.field_force_y->update_host();
    field_controller.field_force_z->update_host();
    field_controller.field_kappa_t->update_host();
    field_controller.field_gamma_t->update_host();
    VTKWriter::write_numerical_debug(field_controller, filename);
}

void Visual::initialise_grid(real *x_coords, real *y_coords, real *z_coords, int Nx, int Ny, int Nz, real dx, real dy, real dz) {
    Domain *domain = Domain::getInstance();
    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    // Initialize grid
    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                size_t index = IX(i, j, k, Nx, Ny);
                x_coords[index] = (X1 + (i - 0.5) * dx);
                y_coords[index] = (Y1 + (j - 0.5) * dy);
                z_coords[index] = (Z1 + (k - 0.5) * dz);
            }
        }
    }
}

std::string Visual::create_filename(const std::string &filename,
                                    int counter, bool analytical) {
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
