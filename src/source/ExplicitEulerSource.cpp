/// \file       ExplicitSource.cpp
/// \brief      Adding source via Explicit Euler
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "ExplicitEulerSource.h"
#include "../boundary/BoundaryController.h"

ExplicitEulerSource::ExplicitEulerSource(Settings const &settings) {
    m_dt = settings.get_real("physical_parameters/dt");
    m_dir_vel = settings.get("solver/source/dir");

    if (m_dir_vel.find('x') == std::string::npos &&
        m_dir_vel.find('y') == std::string::npos &&
        m_dir_vel.find('z') == std::string::npos) {
#ifndef BENCHMARKING
        m_logger = Utility::create_logger(settings, typeid(ExplicitEulerSource).name());
        m_logger->error("unknown direction -> exit");
#endif
        std::exit(1);
        // TODO(issue 6) Error handling
    }
}

//==================================== Add Source ======================================
// ***************************************************************************************
/// \brief  adds source \f$ \partial_t \phi_1 = S_{\phi} \f$ via Explicit Euler
/// \param  out_x  output pointer in x-direction
/// \param  out_y  output pointer in y-direction
/// \param  out_z  output pointer in z-direction
/// \param  S_x    Source pointer in x-direction
/// \param  S_y    Source pointer in y-direction
/// \param  S_z    Source pointer in z-direction
/// \param  sync  synchronous kernel launching (true, default: false)
// ***************************************************************************************
void ExplicitEulerSource::add_source(
        Field &out_x, Field &out_y, Field &out_z,
        const Field &s_x, const Field &s_y, const Field &s_z,
        bool sync) {
    auto dt = m_dt;
    auto dir = m_dir_vel;

    auto boundary = BoundaryController::getInstance();
    size_t *d_inner_list = boundary->get_inner_list_level_joined();
    auto bsize_i = boundary->get_size_inner_list();

#pragma acc data present(out_x, out_y, out_z, s_x, s_y, s_z)
    {
        // check directions of source
        // x - direction
        if (dir.find('x') != std::string::npos) {
#pragma acc parallel loop independent present(out_x, s_x, d_inner_list[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_inner_list[j];
                out_x[i] += dt * s_x[i];
            }

            boundary->apply_boundary(out_x, sync);
        } // end x- direction

        // y - direction
        if (dir.find('y') != std::string::npos) {
#pragma acc parallel loop independent present(out_y, s_y, d_inner_list[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_inner_list[j];
                out_y[i] += dt * s_y[i];
            }

            boundary->apply_boundary(out_y, sync);
        } // end y- direction

        // z - direction
        if (dir.find('z') != std::string::npos) {
#pragma acc parallel loop independent present(out_z, s_z, d_inner_list[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_inner_list[j];
                out_z[i] += dt * s_z[i];
            }

            boundary->apply_boundary(out_z, sync);
        } // end z- direction

        if (sync) {
#pragma acc wait
        }
    }
}

// ***************************************************************************************
/// \brief  adds source \f$ \partial_t \phi_1 = S_{\phi} \f$ via Explicit Euler
/// \param  out   output pointer
/// \param  S   Source pointer
/// \param  sync  synchronous kernel launching (true, default: false)
// ***************************************************************************************
void ExplicitEulerSource::add_source(Field &out, Field const &s, bool sync) {
    auto boundary = BoundaryController::getInstance();
    size_t *d_inner_list = boundary->get_inner_list_level_joined();
    auto bsize_i = boundary->get_size_inner_list();

#pragma acc data present(out, s)
    {
#pragma acc parallel loop independent present(out, s, d_inner_list[:bsize_i]) async
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_inner_list[j];
            out[i] += m_dt * s[i];
        }

        boundary->apply_boundary(out, sync);

        if (sync) {
#pragma acc wait
        }
    }
}
