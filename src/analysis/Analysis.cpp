/// \file       Analysis.cpp
/// \brief      Calculates residual, compares analytical and numerical solutions, saves variables
/// \date       Jul 11, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <vector>
#include <cmath>
#include <fstream>

#include "Analysis.h"
#include "../boundary/BoundaryController.h"
#include "../DomainData.h"
#include "../utility/Utility.h"

Analysis::Analysis(Settings::Settings const &settings, Solution &solution, bool has_analytical_solution) :
        m_settings(settings),
        m_has_analytic_solution(has_analytical_solution),
        m_solution(solution) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(m_settings, typeid(this).name());
#endif
    if (m_has_analytic_solution) {
        m_tol = m_settings.get_real("solver/solution/tol");
    } else {
#ifndef BENCHMARKING
        m_logger->info("No analytical solution available!");
#endif
    }
}

// ================================ Start analysis =============================
// *****************************************************************************
/// \brief  starts analysis to compare numerical and analytical solutions
/// \param  field_controller    pointer to solver
/// \param  t     current time
// ***************************************************************************************
void Analysis::analyse(FieldController *field_controller, real t) {
    if (!m_has_analytic_solution) {
        return;
    }
    m_solution.calc_analytical_solution(t);

#ifndef BENCHMARKING
    m_logger->info("Compare to analytical solution:");
#endif

    auto used_fields = BoundaryController::getInstance()->get_used_fields();
    for (FieldType ft : used_fields) {
        switch (ft) {
            case FieldType::U:
                compare_solutions(field_controller->get_field_u_data(),
                                  m_solution.get_return_ptr_data_u(),
                                  ft, t);
                break;
            case FieldType::V:
                compare_solutions(field_controller->get_field_v_data(),
                                  m_solution.get_return_ptr_data_v(),
                                  ft, t);
                break;
            case FieldType::W:
                compare_solutions(field_controller->get_field_w_data(),
                                  m_solution.get_return_ptr_data_w(),
                                  ft, t);
                break;
            case FieldType::P:
                compare_solutions(field_controller->get_field_p_data(),
                                  m_solution.get_return_ptr_data_p(),
                                  ft, t);
                break;
            case FieldType::T:
                compare_solutions(field_controller->get_field_T_data(),
                                  m_solution.get_return_ptr_data_T(),
                                  ft, t);
                break;
            default:  // do nothing
                break;
        }
    }
}

// ================== Compare analytical and numerical solution ================
// *****************************************************************************
/// \brief  compares analytical solution and numerical solution, returns true when verification passed
/// \param  num   numerical solution
/// \param  ana   analytical solution
/// \param  type  type of variable
/// \param  t   current time
// ***************************************************************************************
bool Analysis::compare_solutions(read_ptr num, read_ptr ana, FieldType type, real t) {
    bool verification = false;

    // Choose absolute or relative based error calculation
    real res = calc_absolute_spatial_error(num, ana);
    //real res = calc_relative_spatial_error(num, ana);

    if (res <= m_tol) {
#ifndef BENCHMARKING
        m_logger->info("{} PASSED Test at time {} with error e = {}",
                       Field::get_field_type_name(type), t, res);
#endif
        verification = true;
    } else {
#ifndef BENCHMARKING
        m_logger->warn("{} FAILED Test at time {} with error e = {}",
                       Field::get_field_type_name(type), t, res);
#endif
    }
    return verification;
}

// ============================= Calculate absolute error ======================
// *****************************************************************************
/// \brief  calculates absolute spatial error based on L2-norm
/// \param  num   numerical solution
/// \param  ana   analytical solution
// ***************************************************************************************
real Analysis::calc_absolute_spatial_error(read_ptr num, read_ptr ana) {
    real sum = 0.;

    auto boundary = BoundaryController::getInstance();
    size_t *inner_list = boundary->get_domain_inner_list_level_joined();
    size_t size_inner_list = boundary->get_size_domain_inner_list_level_joined(0);

    //// weighted 2-norm
    // absolute error
    for (size_t i = 0; i < size_inner_list; i++) {
        size_t idx = inner_list[i];
        real r = std::fabs(num[idx] - ana[idx]);
        sum += r * r;
    }

    // weight
    real nr = static_cast<real>(size_inner_list);
    real eps = std::sqrt(1. / nr * sum);

    // TODO(n16h7) add. output
#ifndef BENCHMARKING
    m_logger->info("Absolute error ||e|| = {}", eps);
#endif
    //std::cout << "num =" << num[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx, m_ny)]        << std::endl;
    //std::cout << "ana =" << ana[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx, m_ny)]        << std::endl;
    //std::cout << "num =" << num[IX((m_nx-2)/2 + 1, (m_ny-2)/2, 1, m_nx, m_ny)]    << std::endl;
    //std::cout << "ana =" << ana[IX((m_nx-2)/2 + 1, (m_ny-2)/2, 1, m_nx, m_ny)]    << std::endl;
    return eps;
}

// ============================= Calculate relative error ======================
// *****************************************************************************
/// \brief  calculates relative spatial error based on L2-norm
/// \param  num   numerical solution
/// \param  ana   analytical solution
// ***************************************************************************************
real Analysis::calc_relative_spatial_error(read_ptr num, read_ptr ana) {
    real sumr = 0.;
    real rr;

    auto boundary = BoundaryController::getInstance();
    size_t *inner_list = boundary->get_domain_inner_list_level_joined();
    size_t size_inner_list = boundary->get_size_domain_inner_list_level_joined(0);

    // relative part with norm of analytical solution as denominator
    for (size_t i = 0; i < size_inner_list; i++) {
        rr = ana[inner_list[i]];
        sumr += rr * rr;
    }

    // weight
    real nr = static_cast<real>(size_inner_list);
    real adenom = sqrt(1. / nr * sumr);

    real eps;
    real zero_tol = 10e-20;
    real epsa = calc_absolute_spatial_error(num, ana);

    // zero absolute error => zero relative error
    if (epsa <= zero_tol) {
        eps = 0.0;

        // zero denominator => take 2-norm of numerical solution as denominator
    } else if (adenom <= zero_tol) {
        sumr = 0.;

        // relative part with norm of numerical solution as quotient
        for (size_t i = 0; i < size_inner_list; i++) {
            rr = num[inner_list[i]];
            sumr += rr * rr;
        }

        real ndenom = sqrt(1. / nr * sumr);

        eps = epsa / ndenom;

    } else {
        eps = epsa / adenom;
    }

    // TODO(n16h7) add. output
#ifndef BENCHMARKING
    m_logger->info("Relative error ||e|| = {}", eps);
#endif
    //std::cout << "num =" << num[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx, m_ny)]        << std::endl;
    //std::cout << "ana =" << ana[IX((m_nx-2)/2, (m_ny-2)/2, 1, m_nx, m_ny)]      << std::endl;
    //std::cout << "num =" << num[IX((m_nx-2)/2 + 1, (m_ny-2)/2, 1, m_nx, m_ny)]  << std::endl;
    //std::cout << "ana =" << ana[IX((m_nx-2)/2 + 1, (m_ny-2)/2, 1, m_nx, m_ny)]  << std::endl;

    return eps;
}

// ======== Calculate absolute error at center to be averaged over time ========
// *****************************************************************************
/// \brief  calculates absolute spatial error at time t at midpoint based on L2-norm
/// \param  field_controller    pointer to field_controller
/// \param  t     current time
/// \param  sum     pointer to sum for (u,p,T results)
// ***************************************************************************************
void Analysis::calc_L2_norm_mid_point(FieldController *field_controller, real t, real *sum) {
    auto boundary = BoundaryController::getInstance();
    size_t *inner_list = boundary->get_domain_inner_list_level_joined();

    size_t ix = inner_list[boundary->get_size_domain_inner_list_level_joined(0) / 2];
    //take median of indices in inner_list to get center point ix
    //std::nth_element(inner_list.begin(), inner_list.begin() + inner_list.size()/2, inner_list.end());
    //size_t ix = inner_list[inner_list.size()/2];

    if (m_has_analytic_solution) {
        m_solution.calc_analytical_solution(t);

        // local variables and parameters
        auto d_ua = m_solution.get_return_ptr_data_u();
        auto d_pa = m_solution.get_return_ptr_data_p();
        auto d_Ta = m_solution.get_return_ptr_data_T();

        auto d_u = field_controller->get_field_u_data();
        auto d_p = field_controller->get_field_p_data();
        auto d_T = field_controller->get_field_T_data();

        real ru = fabs((d_u[ix] - d_ua[ix]));
        real rp = fabs((d_p[ix] - d_pa[ix]));
        real rT = fabs((d_T[ix] - d_Ta[ix]));
        sum[0] += ru * ru;
        sum[1] += rp * rp;
        sum[2] += rT * rT;
    }
}

// =========================== Calculate RMS error ============================
// *****************************************************************************
/// \brief  calculates absolute spatial error at time t at midpoint based on L2-norm
/// \param  solver    pointer to solver
/// \param  t     current time
/// \param  sum     pointer to sum for (u,p,T results)
// ***************************************************************************************
void Analysis::calc_RMS_error(real sum_u, real sum_p, real sum_T) {
    if (!m_has_analytic_solution) {
        return;
    }

    if (m_has_analytic_solution) {
        // local variables and parameters
        real dt = m_settings.get_real("physical_parameters/dt");
        real t_end = m_settings.get_real("physical_parameters/t_end");
        auto Nt = static_cast<size_t>(std::round(t_end / dt));
        real rNt = 1. / static_cast<real>(Nt);

        real epsu = sqrt(rNt * sum_u);

#ifndef BENCHMARKING
        m_logger->info("RMS error of u at domain center is e_RMS = {}", epsu);
#endif

        std::vector<FieldType> v_fields = BoundaryController::getInstance()->get_used_fields();
        if (std::count(v_fields.begin(), v_fields.end(), FieldType::P)) {
            real epsp = sqrt(rNt * sum_p);
#ifndef BENCHMARKING
            m_logger->info("RMS error of p at domain center is e_RMS = {}", epsp);
#endif
        }
        if(std::count(v_fields.begin(), v_fields.end(), FieldType::T)) {
            real epsT = sqrt(rNt * sum_T);
#ifndef BENCHMARKING
            m_logger->info("RMS error of T at domain center is e_RMS = {}", epsT);
#endif
        }
    }
}

// ========================== Check Von Neumann condition ======================
// *****************************************************************************
/// \brief  checks Von Neumann condition on time step (returns true or false)
/// \param  u     x-velocity field
/// \param  dt      time step size
// ***************************************************************************************
bool Analysis::check_time_step_VN(const real dt) {
    bool VN_check;
    auto domain = DomainData::getInstance();

    // local variables and parameters
    real nu = m_settings.get_real("physical_parameters/nu");

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    real dx2sum = (dx * dx + dy * dy + dz * dz);
    real rdx2 = 1. / dx2sum;

    real VN = dt * nu * rdx2;

    VN_check = VN < 0.5;

#ifndef BENCHMARKING
    m_logger->info("VN = {}", VN);
#endif

    return VN_check;
}

// ================================ Calc CFL ===================================
// *****************************************************************************
/// \brief  Returns the max CFL overall cells
/// \param  u     x-velocity field
/// \param  v     y-velocity field
/// \param  w     z-velocity field
/// \param  dt      time step size
// *****************************************************************************
real Analysis::calc_CFL(Field const &u, Field const &v, Field const &w, real dt) const {
    real cfl_max = 0;  // highest seen C. C is always positive, so 0 is a lower bound
    real cfl_local;    // C in the local cell

    auto boundary = BoundaryController::getInstance();
    auto domain = DomainData::getInstance();

    // local variables and parameters
    size_t *inner_list = boundary->get_domain_inner_list_level_joined();
    size_t size_inner_list = boundary->get_size_domain_inner_list_level_joined(0);

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    // calc C for every cell and get the maximum
#pragma acc data present(u, v, w)
#pragma acc parallel loop reduction(max:cfl_max)
    for (size_t i = 0; i < size_inner_list; i++) {
        size_t idx = inner_list[i];
        // \frac{C}{\Delta t} = \frac{\Delta u}{\Delta x} +
        //                      \frac{\Delta v}{\Delta y} +
        //                      \frac{\Delta w}{\Delta z} +
        // https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition#The_two_and_general_n-dimensional_case
        cfl_local = std::fabs(u[idx]) / dx +
                    std::fabs(v[idx]) / dy +
                    std::fabs(w[idx]) / dz;
        cfl_max = std::max(cfl_max, cfl_local);
    }

    return dt * cfl_max;
}

// =============================== Save variables ==============================
// *****************************************************************************
/// \brief  saves variables in .dat files
/// \param  field_controller    pointer to field controller
// ***************************************************************************************
void Analysis::save_variables_in_file(FieldController *field_controller) {
    auto boundary = BoundaryController::getInstance();
    std::vector<FieldType> v_fields = boundary->get_used_fields();

    Field *fields[number_of_field_types];
    fields[FieldType::RHO] = &field_controller->get_field_concentration();
    fields[FieldType::U] = &field_controller->get_field_u();
    fields[FieldType::V] = &field_controller->get_field_v();
    fields[FieldType::W] = &field_controller->get_field_w();
    fields[FieldType::P] = &field_controller->get_field_p();
    fields[FieldType::T] = &field_controller->get_field_T();

    for (auto &v_field: v_fields) {
        write_file(*fields[v_field], Field::get_field_type_name(v_field));
        write_obstacles(*fields[v_field], Field::get_field_type_name(v_field));
    }
}

void Analysis::write_file(const Field &field, const std::string &filename) {
    std::ofstream out;
    out.open(filename + ".dat", std::ofstream::out);
    size_t size = field.get_size();
    real *data = field.data;
    for (size_t index = 0; index < size; index++) {
        out << data[index] << std::endl;
    }
    out.close();
}

void Analysis::write_obstacles(const Field &field, const std::string &filename) {
    BoundaryController *boundary = BoundaryController::getInstance();
    size_t *obstacle_list = boundary->get_obstacle_list_level_joined();
    size_t size = boundary->get_slice_size_obstacle_list_level_joined(0);
    real *data = field.data;
    if (size > 0) {  // do not create (empty) file if there are no obstacles
        size_t start = boundary->get_obstacle_list_level_joined_start(0);
        size_t end = boundary->get_obstacle_list_level_joined_end(0);

        std::ofstream out_obstacle;
        out_obstacle.open(filename + "_obstacle.dat", std::ofstream::out);
        for (size_t idx = start; idx <= end; idx++) {
            out_obstacle << obstacle_list[idx] << ";" << data[obstacle_list[idx]] << std::endl;
        }
        out_obstacle.close();
    }
}
