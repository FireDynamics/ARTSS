/// \file       TimeIntegration.cpp
/// \brief      Runs the time loop
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <chrono>
#include <vector>
#include <iostream>

#include "TimeIntegration.h"
#include "utility/Parameters.h"
#include "Domain.h"

// ==================================== Constructor ====================================
// ***************************************************************************************
/// \brief  Constructor
/// \param  solver_controller class representative for solver
// ***************************************************************************************
TimeIntegration::TimeIntegration(SolverController *sc) {
    auto params = Parameters::getInstance();

    m_dt = params->get_real("physical_parameters/dt");
    m_t_end = params->get_real("physical_parameters/t_end");
    m_t_cur = m_dt;        // since t=0 already handled in setup

    m_solver_controller = sc;
    m_field_controller = m_solver_controller->get_field_controller();

    m_adaption = new Adaption(m_field_controller);
#ifndef BENCHMARKING
    m_solution = new Solution();
    m_analysis = new Analysis(m_solution);
    m_visual = new Visual(m_solution);
#endif
}

void TimeIntegration::set_up() {
#ifndef BENCHMARKING
    m_analysis->analyse(m_field_controller, 0.);
    m_visual->visualise(m_field_controller, 0.);
#endif
}

void TimeIntegration::run() {
    Domain *domain = Domain::getInstance();
    set_up();
    std::cout << "Start calculating and timing...\n" << std::endl;
    // TODO Logger

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    {
        // local variables and parameters for GPU
        auto u = m_field_controller->field_u;
        auto v = m_field_controller->field_v;
        auto w = m_field_controller->field_w;
        auto u0 = m_field_controller->field_u0;
        auto v0 = m_field_controller->field_v0;
        auto w0 = m_field_controller->field_w0;
        auto u_tmp = m_field_controller->field_u_tmp;
        auto v_tmp = m_field_controller->field_v_tmp;
        auto w_tmp = m_field_controller->field_w_tmp;
        auto p = m_field_controller->field_p;
        auto p0 = m_field_controller->field_p0;
        auto rhs = m_field_controller->field_rhs;
        auto T = m_field_controller->field_T;
        auto T0 = m_field_controller->field_T0;
        auto T_tmp = m_field_controller->field_T_tmp;
        auto T_a = m_field_controller->field_T_ambient;
        auto C = m_field_controller->field_concentration;
        auto C0 = m_field_controller->field_concentration0;
        auto C_tmp = m_field_controller->field_concentration_tmp;
        auto f_x = m_field_controller->field_force_x;
        auto f_y = m_field_controller->field_force_y;
        auto f_z = m_field_controller->field_force_z;
        auto S_T = m_field_controller->field_source_T;
        auto S_C = m_field_controller->field_source_concentration;
        auto nu_t = m_field_controller->field_nu_t;
        auto kappa_t = m_field_controller->field_kappa_t;
        auto gamma_t = m_field_controller->field_gamma_t;

        auto d_u = u->data;
        auto d_v = v->data;
        auto d_w = w->data;
        auto d_u0 = u0->data;
        auto d_v0 = v0->data;
        auto d_w0 = w0->data;
        auto d_u_tmp = u_tmp->data;
        auto d_v_tmp = v_tmp->data;
        auto d_w_tmp = w_tmp->data;
        auto d_p = p->data;
        auto d_p0 = p0->data;
        auto d_rhs = rhs->data;
        auto d_T = T->data;
        auto d_T0 = T0->data;
        auto d_T_tmp = T_tmp->data;
        auto d_T_a = T_a->data;
        auto d_C = C->data;
        auto d_C0 = C0->data;
        auto d_C_tmp = C_tmp->data;
        auto d_f_x = f_x->data;
        auto d_f_y = f_y->data;
        auto d_f_z = f_z->data;
        auto d_S_T = S_T->data;
        auto d_S_C = S_C->data;
        auto d_nu_t = nu_t->data;
        auto d_kappa_t = kappa_t->data;
        auto d_gamma_t = gamma_t->data;

        auto bsize = domain->get_size();

        auto t_cur = m_t_cur;
        auto t_end = m_t_end;
        auto dt = m_dt;

        // preparation RMS error
#ifndef BENCHMARKING
        real Sumu = 0., Sump = 0., SumT = 0.;
        real Sum[3];
        Sum[0] = Sumu, Sum[1] = Sump, Sum[2] = SumT;

        /*
        bool CFL_check = m_analysis->check_time_step_CFL(u, v, w, dt);
        if (!CFL_check) {
            real cfl_dt = m_analysis->set_DT_with_CFL(u, v, w);
            std::cout << "CFL condition not met! A corresponding dt would be: " << cfl_dt << std::endl;
            //TODO Logger warning
        }
        bool VN_check = m_analysis->check_time_step_VN(u, dt);
        if (!VN_check) {
            std::cout << "Von Neumann condition not met!" << std::endl;
            //TODO Logger warning
        }
         */
#endif

// copyin all variables
#pragma acc enter data copyin(  d_u[:bsize], d_u0[:bsize], d_u_tmp[:bsize], \
                                d_v[:bsize], d_v0[:bsize], d_v_tmp[:bsize], \
                                d_w[:bsize], d_w0[:bsize], d_w_tmp[:bsize], \
                                d_p[:bsize], d_p0[:bsize], d_rhs[:bsize], \
                                d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], d_T_a[:bsize], \
                                d_C[:bsize], d_C0[:bsize], d_C_tmp[:bsize], \
                                d_f_x[:bsize], d_f_y[:bsize], d_f_z[:bsize], d_S_T[:bsize], d_S_C[:bsize], \
                                d_nu_t[:bsize], d_kappa_t[:bsize], d_gamma_t[:bsize])

        int iteration_step = 1;
        // std::ofstream file;
        // file.open(adaption->get_write_runtime_name(), ios::app);
        // std::chrono::time_point<std::chrono::system_clock> iter_start, iter_end;
        while (t_cur < t_end + dt / 2) {
            //iter_start = std::chrono::system_clock::now();
#ifndef BENCHMARKING
            std::cout << "\nt_cur=" << t_cur << std::endl;
            //TODO Logger
#endif

            // Calculate
            m_solver_controller->solver_do_step(t_cur, false);
#ifndef BENCHMARKING
            // Visualize
#pragma acc update host(d_u[:bsize])
#pragma acc update host(d_v[:bsize])
#pragma acc update host(d_w[:bsize])
#pragma acc update host(d_p[:bsize])
#pragma acc update host(d_rhs[:bsize])
#pragma acc update host(d_T[:bsize])
#pragma acc update host(d_C[:bsize])
#pragma acc update host(d_nu_t[:bsize])
#pragma acc update host(d_S_T[:bsize]) wait    // all in one update does not work!

            m_visual->visualise(m_field_controller, t_cur);
            if (m_adaption->is_data_extraction_before_enabled()) m_adaption->extractData(m_adaption->get_before_name(), m_adaption->get_before_height(), t_cur);
            // Error Calculation
            // RMS error at midposize_t at each time step Nt
            m_analysis->calc_L2_norm_mid_point(m_field_controller, t_cur, Sum);
#endif
            // update
            m_adaption->run(t_cur);
#ifndef BENCHMARKING
            if (m_adaption->is_data_extraction_after_enabled()) m_adaption->extractData(m_adaption->get_after_name(), m_adaption->get_after_height(), t_cur);
#endif
            m_solver_controller->update_sources(t_cur, false);
            m_field_controller->update_data(false);

            // iter_end = std::chrono::system_clock::now();
            // long ms = std::chrono::duration_cast<std::chrono::microseconds>(iter_end - iter_start).count();
            // file << "t_cur: "<<t_cur << " runtime: " << ms << " microsec\n";
            iteration_step++;
            t_cur = iteration_step * dt;

        } // end while
        // file.close();
        // Sum up RMS error
#ifndef BENCHMARKING
        m_analysis->calc_RMS_error(Sum[0], Sum[1], Sum[2]);
#endif

#pragma acc wait
// delete unnecessary variables copyout numerical solution
#pragma acc exit data delete(d_u0[:bsize], d_u_tmp[:bsize], d_v0[:bsize], d_v_tmp[:bsize], d_w0[:bsize], d_w_tmp[:bsize], d_p0[:bsize], d_T0[:bsize], d_T_tmp[:bsize], d_T_a[:bsize], d_C0[:bsize], d_C_tmp[:bsize], d_f_x[:bsize], d_f_y[:bsize], d_f_z[:bsize], d_S_T[:bsize], d_S_C[:bsize], d_nu_t[:bsize], d_kappa_t[:bsize], d_gamma_t[:bsize])
#pragma acc exit data copyout(d_u[:bsize], d_v[:bsize], d_w[:bsize], d_p[:bsize], d_rhs[:bsize], d_T[:bsize], d_C[:bsize])

    } // end RANGE

    std::cout << "\nDone calculating and timing ...\n" << std::endl;
    //TODO Logger

    // stop timer
    end = std::chrono::system_clock::now();
    long ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Global Time: " << ms << "ms" << std::endl;
    //TODO Logger
#ifndef BENCHMARKING
    if (m_adaption->is_data_extraction_endresult_enabled()) {
        m_adaption->extractData(m_adaption->get_endresult_name());
    }
    // testing correct output (when changing implementation/ calculating on GPU)
    m_analysis->save_variables_in_file(m_field_controller);
    m_analysis->analyse(m_field_controller, m_t_end);
#endif
    delete m_adaption;
    delete m_analysis;
    delete m_solution;
    delete m_visual;
}
