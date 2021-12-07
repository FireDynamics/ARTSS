/// \file       TimeIntegration.cpp
/// \brief      Runs the time loop
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "TimeIntegration.h"
#include "DomainData.h"

// ==================================== Constructor ====================================
// ***************************************************************************************
/// \brief  Constructor
/// \param  solver_controller class representative for solver
// ***************************************************************************************
TimeIntegration::TimeIntegration(Settings::Settings const &settings, SolverController *sc) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    m_dt = settings.get_real("physical_parameters/dt");
    m_t_end = settings.get_real("physical_parameters/t_end");
    m_t_cur = m_dt;        // since t=0 already handled in setup

    m_solver_controller = sc;
    m_field_controller = m_solver_controller->get_field_controller();

    m_adaption = new Adaption(settings, m_field_controller);
#ifndef BENCHMARKING
    std::string initial_condition = settings.get("initial_conditions/usr_fct");
    bool has_analytical_solution = settings.get_bool("solver/solution/available");
    m_solution = new Solution(settings, initial_condition, has_analytical_solution);
    m_analysis = new Analysis(settings, *m_solution, has_analytical_solution);
    m_visual = new Visual(settings, *m_solution, has_analytical_solution);
#endif
}

void TimeIntegration::run() {
    Field &u = m_field_controller->get_field_u();
    Field &v = m_field_controller->get_field_v();
    Field &w = m_field_controller->get_field_w();

#ifndef BENCHMARKING
    m_field_controller->update_host();
    m_analysis->analyse(m_field_controller, 0.);
    m_visual->visualise(*m_field_controller, 0.);
    Visual::write_vtk_debug(*m_field_controller, "initial_steckler");
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
            m_field_controller->update_host();
            // Visualize
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
            real cfl = m_analysis->calc_CFL(u, v, w, dt);

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

    }
    // stop timer
    end = std::chrono::system_clock::now();
    long ms = std::chrono::duration_cast < std::chrono::milliseconds > (end - start).count();

#ifndef BENCHMARKING
    m_logger->info("Done calculating and timing ...");
    m_logger->info("Global Time: {}ms", ms);
    m_field_controller->update_host();
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
    std::cout << "Done calculating and timing ..." << std::endl;
    std::cout << "Global Time: " << ms << "ms" << std::endl;
#endif
    delete m_adaption;
}
