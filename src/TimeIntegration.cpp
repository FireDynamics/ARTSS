/// \file       TimeIntegration.cpp
/// \brief      Runs the time loop
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "TimeIntegration.h"
#include "utility/Parameters.h"
#include "Domain.h"
#include "visualisation/VTKWriter.h"

// ==================================== Constructor ====================================
// ***************************************************************************************
/// \brief  Constructor
/// \param  solver_controller class representative for solver
// ***************************************************************************************
TimeIntegration::TimeIntegration(SolverController *sc) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    auto params = Parameters::getInstance();

    m_dt = params->get_real("physical_parameters/dt");
    m_t_end = params->get_real("physical_parameters/t_end");
    m_t_cur = m_dt;        // since t=0 already handled in setup

    m_solver_controller = sc;
    m_field_controller = m_solver_controller->get_field_controller();
    int counter = 0;
    VTKWriter::write_numerical(*m_field_controller, "debug_ti_" + std::to_string(counter));
    counter++;

    m_adaption = new Adaption(m_field_controller);
    VTKWriter::write_numerical(*m_field_controller, "debug_ti_" + std::to_string(counter));
    counter++;
#ifndef BENCHMARKING
    m_solution = new Solution();
    VTKWriter::write_numerical(*m_field_controller, "debug_ti_" + std::to_string(counter));
    counter++;
    m_analysis = new Analysis(m_solution);
    VTKWriter::write_numerical(*m_field_controller, "debug_ti_" + std::to_string(counter));
    counter++;
    m_visual = new Visual(*m_solution);
    VTKWriter::write_numerical(*m_field_controller, "debug_ti_" + std::to_string(counter));
    counter++;
#endif
}

void TimeIntegration::run() {
    Domain *domain = Domain::getInstance();

    // local variables and parameters for GPU
    auto u = m_field_controller->field_u;
    auto v = m_field_controller->field_v;
    auto w = m_field_controller->field_w;
    auto p = m_field_controller->field_p;
    auto rhs = m_field_controller->field_rhs;
    auto T = m_field_controller->field_T;
    auto C = m_field_controller->field_concentration;
    auto S_T = m_field_controller->field_source_T;
    auto S_C = m_field_controller->field_source_concentration;
    auto nu_t = m_field_controller->field_nu_t;

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

    auto bsize = domain->get_size();

#ifndef BENCHMARKING
#pragma acc update host(d_u[:bsize])
#pragma acc update host(d_v[:bsize])
#pragma acc update host(d_w[:bsize])
#pragma acc update host(d_p[:bsize])
#pragma acc update host(d_rhs[:bsize])
#pragma acc update host(d_T[:bsize])
#pragma acc update host(d_C[:bsize])
#pragma acc update host(d_nu_t[:bsize])
#pragma acc update host(d_S_T[:bsize]) wait    // all in one update does not work!
    int counter = 0;
    VTKWriter::write_numerical(*m_field_controller, "debug_ti2_" + std::to_string(counter));
    counter++;
    m_analysis->analyse(m_field_controller, 0.);
    VTKWriter::write_numerical(*m_field_controller, "debug_ti2_" + std::to_string(counter));
    counter++;
    m_visual->visualise(*m_field_controller, 0.);
    VTKWriter::write_numerical(*m_field_controller, "debug_ti2_" + std::to_string(counter));
    counter++;
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
#pragma acc update host(d_u[:bsize])
#pragma acc update host(d_v[:bsize])
#pragma acc update host(d_w[:bsize])
#pragma acc update host(d_p[:bsize])
#pragma acc update host(d_rhs[:bsize])
#pragma acc update host(d_T[:bsize])
#pragma acc update host(d_C[:bsize])
#pragma acc update host(d_nu_t[:bsize])
#pragma acc update host(d_S_T[:bsize]) wait    // all in one update does not work!

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
            // bool VN_check = ana.check_time_step_VN(u, dt);
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

#pragma acc update host(d_u[:bsize])
#pragma acc update host(d_v[:bsize])
#pragma acc update host(d_w[:bsize])
#pragma acc update host(d_p[:bsize])
#pragma acc update host(d_rhs[:bsize])
#pragma acc update host(d_T[:bsize])
#pragma acc update host(d_C[:bsize]) wait

    } // end RANGE

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
