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

void Visual::write_csv(const FieldController &field_controller, std::string filename){
    // local variables and parameters for GPU
    auto u = field_controller.field_u;
    auto v = field_controller.field_v;
    auto w = field_controller.field_w;
    auto p = field_controller.field_p;
    auto rhs = field_controller.field_rhs;
    auto T = field_controller.field_T;
    auto C = field_controller.field_concentration;
    auto S_T = field_controller.field_source_T;
    auto S_C = field_controller.field_source_concentration;
    auto nu_t = field_controller.field_nu_t;

    auto d_u = u->data;
    auto d_v = v->data;
    auto d_w = w->data;
    auto d_p = p->data;
    auto d_rhs = rhs->data;
    auto d_T = T->data;
    auto d_C = C->data;
    auto d_S_T = S_T->data;
    auto d_S_C = S_C->data;
    auto d_nu_t = nu_t->data;

    auto bsize = Domain::getInstance()->get_size();
#pragma acc update host(d_u[:bsize])
#pragma acc update host(d_v[:bsize])
#pragma acc update host(d_w[:bsize])
#pragma acc update host(d_p[:bsize])
#pragma acc update host(d_rhs[:bsize])
#pragma acc update host(d_T[:bsize])
#pragma acc update host(d_C[:bsize])
#pragma acc update host(d_nu_t[:bsize])
#pragma acc update host(d_S_T[:bsize]) wait
    CSVWriter::write_numerical(field_controller, filename);
}

void Visual::write_vtk(const FieldController &field_controller, std::string filename){
    // local variables and parameters for GPU
    auto u = field_controller.field_u;
    auto v = field_controller.field_v;
    auto w = field_controller.field_w;
    auto p = field_controller.field_p;
    auto rhs = field_controller.field_rhs;
    auto T = field_controller.field_T;
    auto C = field_controller.field_concentration;
    auto S_T = field_controller.field_source_T;
    auto S_C = field_controller.field_source_concentration;
    auto nu_t = field_controller.field_nu_t;

    auto d_u = u->data;
    auto d_v = v->data;
    auto d_w = w->data;
    auto d_p = p->data;
    auto d_rhs = rhs->data;
    auto d_T = T->data;
    auto d_C = C->data;
    auto d_S_T = S_T->data;
    auto d_S_C = S_C->data;
    auto d_nu_t = nu_t->data;

    auto bsize = Domain::getInstance()->get_size();
#pragma acc update host(d_u[:bsize])
#pragma acc update host(d_v[:bsize])
#pragma acc update host(d_w[:bsize])
#pragma acc update host(d_p[:bsize])
#pragma acc update host(d_rhs[:bsize])
#pragma acc update host(d_T[:bsize])
#pragma acc update host(d_C[:bsize])
#pragma acc update host(d_nu_t[:bsize])
#pragma acc update host(d_S_T[:bsize]) wait
    VTKWriter::write_numerical(field_controller, filename);
}

void Visual::write_vtk_debug(const FieldController &field_controller, std::string filename){
    // local variables and parameters for GPU
    auto d_u = field_controller.get_field_u_data();
    auto d_v = field_controller.get_field_v_data();
    auto d_w = field_controller.get_field_w_data();
    auto d_p = field_controller.get_field_p_data();
    auto d_rhs = field_controller.get_field_rhs_data();
    auto d_T = field_controller.get_field_T_data();
    auto d_C = field_controller.get_field_concentration_data();
    auto d_S_T = field_controller.get_field_source_T_data();
    auto d_S_C = field_controller.get_field_source_concentration_data();
    auto d_nu_t = field_controller.get_field_nu_t_data();
    auto d_f_x = field_controller.get_field_force_x_data();
    auto d_f_y = field_controller.get_field_force_y_data();
    auto d_f_z = field_controller.get_field_force_z_data();
    auto d_kappa = field_controller.get_field_kappa_data();
    auto d_gamma = field_controller.get_field_gamma_data();

    auto bsize = Domain::getInstance()->get_size();
#pragma acc update host(d_u[:bsize])
#pragma acc update host(d_v[:bsize])
#pragma acc update host(d_w[:bsize])
#pragma acc update host(d_p[:bsize])
#pragma acc update host(d_rhs[:bsize])
#pragma acc update host(d_T[:bsize])
#pragma acc update host(d_C[:bsize])
#pragma acc update host(d_nu_t[:bsize])
#pragma acc update host(d_S_T[:bsize])
#pragma acc update host(d_S_C[:bsize])
#pragma acc update host(d_f_x[:bsize])
#pragma acc update host(d_f_y[:bsize])
#pragma acc update host(d_f_z[:bsize])
#pragma acc update host(d_kappa[:bsize])
#pragma acc update host(d_gamma[:bsize]) wait
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

std::string Visual::create_filename(std::string filename, int counter, bool analytical) {
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
