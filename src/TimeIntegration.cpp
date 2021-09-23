/// \file       TimeIntegration.cpp
/// \brief      Runs the time loop
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "TimeIntegration.h"
#include "utility/Parameters.h"
#include "Domain.h"

// ==================================== Constructor ====================================
// ***************************************************************************************
/// \brief  Constructor
/// \param  solver_controller class representative for solver
// ***************************************************************************************
TimeIntegration::TimeIntegration(SolverController *sc) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    Parameters *params = Parameters::getInstance();

    m_dt = params->get_real("physical_parameters/dt");
    m_t_end = params->get_real("physical_parameters/t_end");
    m_t_cur = m_dt;        // since t=0 already handled in setup

    m_solver_controller = sc;
    m_field_controller = m_solver_controller->get_field_controller();

    m_adaption = new Adaption(m_field_controller);
#ifndef BENCHMARKING
    std::string initial_condition = params->get("initial_conditions/usr_fct");
    bool has_analytical_solution = (params->get("solver/solution/available") == "Yes");
    m_solution = new Solution(initial_condition, has_analytical_solution);
    m_analysis = new Analysis(m_solution);
    m_visual = new Visual(*m_solution);
#endif
}

void TimeIntegration::run() {
    Domain *domain = Domain::getInstance();

    // local variables and parameters for GPU
    Field &u = *m_field_controller->field_u;
    Field &v = *m_field_controller->field_v;
    Field &w = *m_field_controller->field_w;
    Field &p = *m_field_controller->field_p;
    Field &rhs = *m_field_controller->field_rhs;
    Field &T = *m_field_controller->field_T;
    Field &C = *m_field_controller->field_concentration;
    Field &S_T = *m_field_controller->field_source_T;
    Field &S_C = *m_field_controller->field_source_concentration;
    Field &nu_t = *m_field_controller->field_nu_t;

#ifndef BENCHMARKING
    u.update_host();
    v.update_host();
    w.update_host();
    p.update_host();
    rhs.update_host();
    T.update_host();
    C.update_host();
    nu_t.update_host();
    S_T.update_host();
    S_C.update_host();
#pragma acc wait

    m_analysis->analyse(m_field_controller, 0.);
    m_visual->visualise(*m_field_controller, 0.);
    m_logger->info("Start calculating and timing...");
#else
    std::cout << "Start calculating and timing...\n" << std::endl;
#endif
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    {
        auto t_cur = m_t_cur;
        auto t_end = m_t_end;
        auto dt = m_dt;

        // preparation RMS error
#ifndef BENCHMARKING
        real Sumu = 0., Sump = 0., SumT = 0.;
        real Sum[3];
        Sum[0] = Sumu, Sum[1] = Sump, Sum[2] = SumT;
#endif

        int iteration_step = 1;
        // std::ofstream file;
        // file.open(adaption->get_write_runtime_name(), ios::app);
        // std::chrono::time_point<std::chrono::system_clock> iter_start, iter_end;
        while (t_cur < t_end + dt / 2) {
            //iter_start = std::chrono::system_clock::now();
#ifndef BENCHMARKING
            m_logger->info("t_cur = {:.5f}", t_cur);
#endif

            // Calculate
            m_solver_controller->solver_do_step(t_cur, false);
#ifndef BENCHMARKING
            // Visualize
            u.update_host();
            v.update_host();
            w.update_host();
            p.update_host();
            rhs.update_host();
            T.update_host();
            C.update_host();
            nu_t.update_host();
            S_T.update_host();
            S_C.update_host();
#pragma acc wait

            m_solution->calc_analytical_solution(t_cur);

            m_visual->visualise(*m_field_controller, t_cur);
            if (m_adaption->is_data_extraction_before_enabled()) {
                m_adaption->extractData(m_adaption->get_before_name(),
                                        m_adaption->get_before_height(),
                                        t_cur);
            }

            // Error Calculation
            // RMS error at midposize_t at each time step Nt
            m_analysis->calc_L2_norm_mid_point(m_field_controller, t_cur, Sum);

            // check CFL
            real cfl = m_analysis->calc_CFL(&u, &v, &w, dt);

            // CFL condition not met
            if (cfl > 1) {
                m_logger->warn("CFL condition not met. CFL={}, dt={}", cfl, dt);
                m_logger->warn("To lower th CFL value a smaller dt must be selected. Proposed CFL value of <= 0.8 yields to dt <= {}", dt*0.8/cfl);
            } else {
                m_logger->info("CFL = {}", cfl);
            }
            // bool VN_check = ana.check_time_step_VN(dt);
            // if(!VN_check)
            //     std::cout<<"Von Neumann condition not met!"<<std::endl;
#endif
            // update
            m_adaption->run(t_cur);
#ifndef BENCHMARKING
            if (m_adaption->is_data_extraction_after_enabled()) {
                m_adaption->extractData(m_adaption->get_after_name(), m_adaption->get_after_height(), t_cur);
            }
#endif
            m_solver_controller->update_sources(t_cur, false);
            m_field_controller->update_data(false);

            // iter_end = std::chrono::system_clock::now();
            // long ms = std::chrono::duration_cast<std::chrono::microseconds>(iter_end - iter_start).count();
            // file << "t_cur: "<<t_cur << " runtime: " << ms << " microsec\n";
            iteration_step++;
            t_cur = iteration_step * dt;
        }
        // file.close();
        // Sum up RMS error
#ifndef BENCHMARKING
        m_analysis->calc_RMS_error(Sum[0], Sum[1], Sum[2]);
#endif

#pragma acc wait
    u.update_host();
    v.update_host();
    w.update_host();
    p.update_host();
    rhs.update_host();
    T.update_host();
    C.update_host();
    nu_t.update_host();
    S_T.update_host();
    S_C.update_host();
#pragma acc wait
    }

#ifndef BENCHMARKING
    m_logger->info("Done calculating and timing ...");
#else
    std::cout << "Done calculating and timing ..." << std::endl;
#endif

    // stop timer
    end = std::chrono::system_clock::now();
    long ms = std::chrono::duration_cast < std::chrono::milliseconds > (end - start).count();

#ifndef BENCHMARKING
    m_logger->info("Global Time: {}ms", ms);
    if (m_adaption->is_data_extraction_endresult_enabled()) {
        m_adaption->extractData(m_adaption->get_endresult_name());
    }
    // testing correct output (when changing implementation/ calculating on GPU)
    m_analysis->save_variables_in_file(m_field_controller);
    m_analysis->analyse(m_field_controller, m_t_end);
    delete m_analysis;
    delete m_solution;
    delete m_visual;
#else
    std::cout << "Global Time: " << ms << "ms" << std::endl;
#endif
    delete m_adaption;
}
