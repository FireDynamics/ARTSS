/// \file       TimeIntegration.cpp
/// \brief      Runs the time loop
/// \date       May 20, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <chrono>
#include <vector>

#include "TimeIntegration.h"
#include "utility/Parameters.h"
#include "analysis/Analysis.h"
#include "adaption/Adaption.h"
#include "Domain.h"

#ifndef BENCHMARKING

#include "visualisation/Visual.h"

#endif

// =============================== Constructor ===============================
// *****************************************************************************
/// \brief  Constructor
/// \param  isolv pointer to solver
/// \param  fname filename of xml-input (via argument)
// ***************************************************************************************
TimeIntegration::TimeIntegration(ISolver *isolv) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    auto params = Parameters::getInstance();
    auto domain = Domain::getInstance();

    m_dt = params->get_real("physical_parameters/dt");
    m_t_end = params->get_real("physical_parameters/t_end");
    m_t_cur = m_dt;        // since t=0 already handled in setup

    m_size = domain->get_size();

    this->m_solver = isolv;
}

void TimeIntegration::run() {
    Adaption *adaption;
    adaption = new Adaption(m_solver);
#ifndef BENCHMARKING
    Solution *solution = new Solution();
    Analysis *analysis = new Analysis(solution);
    analysis->analyse(m_solver, 0.);

    // Visualise
    Visual *visual = new Visual(solution);
    visual->visualise(m_solver, 0.);
    m_logger->info("Start calculating and timing...");
#endif

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    {
        // local variables and parameters for GPU
        auto u = m_solver->u;
        auto v = m_solver->v;
        auto w = m_solver->w;
        auto u0 = m_solver->u0;
        auto v0 = m_solver->v0;
        auto w0 = m_solver->w0;
        auto u_tmp = m_solver->u_tmp;
        auto v_tmp = m_solver->v_tmp;
        auto w_tmp = m_solver->w_tmp;
        auto p = m_solver->p;
        auto p0 = m_solver->p0;
        auto rhs = m_solver->rhs;
        auto T = m_solver->T;
        auto T0 = m_solver->T0;
        auto T_tmp = m_solver->T_tmp;
        auto T_a = m_solver->T_ambient;
        auto C = m_solver->concentration;
        auto C0 = m_solver->concentration0;
        auto C_tmp = m_solver->concentration_tmp;
        auto f_x = m_solver->f_x;
        auto f_y = m_solver->f_y;
        auto f_z = m_solver->f_z;
        auto S_T = m_solver->S_T;
        auto S_C = m_solver->S_concentration;
        auto nu_t = m_solver->nu_t;
        auto kappa_t = m_solver->kappa_t;
        auto gamma_t = m_solver->gamma_t;

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

        auto bsize = m_size;

        auto t_cur = m_t_cur;
        auto t_end = m_t_end;
        auto dt = m_dt;

        // preparation RMS error
#ifndef BENCHMARKING
        real Sumu = 0., Sump = 0., SumT = 0.;
        real Sum[3];
        Sum[0] = Sumu, Sum[1] = Sump, Sum[2] = SumT;
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

        // initialize boundary cells
        m_solver->set_up_boundary(false);
        int iteration_step = 1;
        // std::ofstream file;
        // file.open(adaption->get_write_runtime_name(), ios::app);
        // std::chrono::time_point<std::chrono::system_clock> iter_start, iter_end;
        while (t_cur < t_end + dt/2) {

            //iter_start = std::chrono::system_clock::now();
#ifndef BENCHMARKING
            m_logger->info("t_cur = {:.5f}", t_cur);
#endif

            // Calculate
            m_solver->do_step(t_cur, false);

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

            visual->visualise(m_solver, t_cur);
            if (adaption->is_data_extraction_before_enabled()) adaption->extractData(adaption->get_before_name(), adaption->get_before_height(), t_cur);
            // Error Calculation
            // RMS error at midposize_t at each time step Nt
            analysis->calc_L2_norm_mid_point(m_solver, t_cur, Sum);

            // To check CFL and VN, comment out
            /*bool CFL_check = ana.check_time_step_CFL(u, v, w, dt);
            if(!CFL_check){
                std::cout<<"CFL condition not met!\n";
                // To change dt, comment out
                //dt = ana.set_DT_with_CFL(u, v, w);
                //std::cout<<" Setting dt = "<<dt<<std::endl;
            }
            bool VN_check = ana.check_time_step_VN(u, dt);
            if(!VN_check)
                std::cout<<"Von Neumann condition not met!"<<std::endl;
            */
#endif
            // update
            adaption->run(t_cur);
#ifndef BENCHMARKING
            if (adaption->is_data_extraction_after_enabled()) adaption->extractData(adaption->get_after_name(), adaption->get_after_height(), t_cur);
#endif
        m_solver->update_sources(t_cur, false);
        m_solver->update_data(false);

            // iter_end = std::chrono::system_clock::now();
            // long ms = std::chrono::duration_cast<std::chrono::microseconds>(iter_end - iter_start).count();
            // file << "t_cur: "<<t_cur << " runtime: " << ms << " microsec\n";
            iteration_step++;
            t_cur = iteration_step * dt;

        } // end while
        // file.close();
        // Sum up RMS error
#ifndef BENCHMARKING
        analysis->calc_RMS_error(Sum[0], Sum[1], Sum[2]);
#endif

#pragma acc wait
// delete unnecessary variables copyout numerical solution
#pragma acc exit data delete(d_u0[:bsize], d_u_tmp[:bsize], d_v0[:bsize], d_v_tmp[:bsize], d_w0[:bsize], d_w_tmp[:bsize], d_p0[:bsize], d_T0[:bsize], d_T_tmp[:bsize], d_T_a[:bsize], d_C0[:bsize], d_C_tmp[:bsize], d_f_x[:bsize], d_f_y[:bsize], d_f_z[:bsize], d_S_T[:bsize], d_S_C[:bsize], d_nu_t[:bsize], d_kappa_t[:bsize], d_gamma_t[:bsize])
#pragma acc exit data copyout(d_u[:bsize], d_v[:bsize], d_w[:bsize], d_p[:bsize], d_rhs[:bsize], d_T[:bsize], d_C[:bsize])

    } // end RANGE

#ifndef BENCHMARKING
    m_logger->info("Done calculating and timing ...");
#endif

    // stop timer
    end = std::chrono::system_clock::now();
    long ms = std::chrono::duration_cast < std::chrono::milliseconds > (end - start).count();

#ifndef BENCHMARKING
    m_logger->info("Global Time: {}ms", ms);
    if (adaption->is_data_extraction_endresult_enabled()) {
        adaption->extractData(adaption->get_endresult_name());
    }
    // testing correct output (when changing implementation/ calculating on GPU)
    analysis->save_variables_in_file(m_solver);
    analysis->analyse(m_solver, m_t_end);
#endif

    delete (adaption);
}
