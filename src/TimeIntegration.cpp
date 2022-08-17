/// \file       TimeIntegration.cpp
/// \brief      Runs the time loop
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <unistd.h>
#include "TimeIntegration.h"
#include "domain/DomainData.h"
#ifdef ASSIMILATION
#include <mpi.h>
#endif

// ==================================== Constructor ====================================
// ***************************************************************************************
/// \brief  Constructor
/// \param  solver_controller class representative for solver
// ***************************************************************************************
TimeIntegration::TimeIntegration(const Settings::Settings &settings, SolverController *sc) :
        m_settings(settings), m_solver_controller(sc) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    auto domain_data = DomainData::getInstance();
    m_dt = domain_data->get_physical_parameters().dt;
    m_t_end = domain_data->get_physical_parameters().t_end;
    m_t_cur = m_dt;        // since t=0 already handled in setup

    m_field_controller = m_solver_controller->get_field_controller();

    m_adaption = new Adaption(m_settings.adaption_parameters, m_field_controller, m_settings.filename);
#ifdef ASSIMILATION
    m_data_assimilation = new DataAssimilation(*m_solver_controller, m_field_controller, m_settings);
#endif
#ifndef BENCHMARKING
    m_solution = new Solution(m_settings.initial_conditions_parameters, m_settings.solver_parameters.solution);
    m_analysis = new Analysis(m_settings.solver_parameters.solution, *m_solution);
    m_visual = new Visual(m_settings.visualisation_parameters, *m_solution, m_settings.filename);
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
#ifdef ASSIMILATION
            if (m_data_assimilation->requires_rollback(t_cur)) {
                t_cur = m_data_assimilation->get_new_time_value();
                iteration_step = static_cast<int>(t_cur / dt);
                m_logger->info("ROLLBACK to time step {} (step: {})", t_cur, iteration_step);
                m_data_assimilation->initiate_rollback(t_cur);
            } else {
                m_data_assimilation->save_data(t_cur);
            }
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

#ifdef ASSIMILATION
    bool simulation_is_running = false;
    MPI_Request request;
    MPI_Isend(&simulation_is_running, 1, MPI_LOGICAL, 1, 77, MPI_COMM_WORLD, &request);
#endif

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
