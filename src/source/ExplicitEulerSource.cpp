/// \file       ExplicitSource.cpp
/// \brief      Adding source via Explicit Euler
/// \date       Dec 2, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "ExplicitEulerSource.h"
#include "../utility/Parameters.h"
#include "../boundary/BoundaryController.h"

ExplicitEulerSource::ExplicitEulerSource() {
    auto params = Parameters::getInstance();

    m_dt = params->get_real("physical_parameters/dt");
    m_dir_vel = params->get("solver/source/dir");

    if (m_dir_vel.find('x') == std::string::npos && m_dir_vel.find('y') == std::string::npos && m_dir_vel.find('z') == std::string::npos) {
#ifndef BENCHMARKING
        m_logger = Utility::create_logger(typeid(ExplicitEulerSource).name());
        m_logger->error("unknown direction -> exit");
#endif
        std::exit(1);
        //TODO Error handling
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
void ExplicitEulerSource::add_source(Field *out_x, Field *out_y, Field *out_z, Field *S_x, Field *S_y, Field *S_z, bool sync) {
    // local variables and parameters for GPU
    size_t level = out_x->get_level();
    FieldType type = out_x->get_type();

    auto d_outx = out_x->data;
    auto d_outy = out_y->data;
    auto d_outz = out_z->data;
    auto d_Sx = S_x->data;
    auto d_Sy = S_y->data;
    auto d_Sz = S_z->data;

    auto dt = m_dt;
    auto dir = m_dir_vel;

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

#pragma acc data present(outx, outy, outz, Sx, Sy, Sz)
    {
        //check directions of source
        //x - direction
        if (dir.find('x') != std::string::npos) {
#pragma acc parallel loop independent present(d_outx, d_Sx, d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
                d_outx[i] += dt * d_Sx[i];
            }

            boundary->applyBoundary(d_outx, level, type, sync);
        }

        // y - direction
        if (dir.find('y') != std::string::npos) {
#pragma acc parallel loop independent present(d_outy, d_Sy, d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
                d_outy[i] += dt * d_Sy[i];
            }

            boundary->applyBoundary(d_outy, level, type, sync);
        }

        // z - direction
        if (dir.find('z') != std::string::npos) {
#pragma acc parallel loop independent present(d_outz, d_Sz, d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
                d_outz[i] += dt * d_Sz[i];
            }

            boundary->applyBoundary(d_outz, level, type, sync);
        }

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
void ExplicitEulerSource::add_source(Field *out, Field *S, bool sync) {
    // local variables and parameters for GPU
    size_t level = out->get_level();
    FieldType type = out->get_type();

    auto d_out = out->data;
    auto d_S = S->data;

    auto dt = m_dt;

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

#pragma acc data present(out, S)
    {
#pragma acc parallel loop independent present(out, S, d_iList[:bsize_i]) async
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            d_out[i] += dt * d_S[i];
        }

        boundary->applyBoundary(d_out, level, type, sync);

        if (sync) {
#pragma acc wait
        }
    }
}
