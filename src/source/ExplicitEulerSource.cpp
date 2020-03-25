/// \file 		ExplicitSource.h
/// \brief 		Adding source via Explicit Euler
/// \date 		Dec 2, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <iostream>
#include "ExplicitEulerSource.h"
#include "../utility/Parameters.h"
#include "../boundary/BoundaryController.h"
#include "../Domain.h"

ExplicitEulerSource::ExplicitEulerSource() {

    auto params = Parameters::getInstance();

    m_dt = params->getReal("physical_parameters/dt");
    m_dir_vel = params->get("solver/source/dir");

    if (m_dir_vel.find('x') == std::string::npos && m_dir_vel.find('y') == std::string::npos && m_dir_vel.find('z') == std::string::npos) {
        std::cout << "unknown direction -> exit" << std::endl;
        std::flush(std::cout);
        std::exit(1);
        //TODO Error handling + Logger
    }
}

//==================================== Add Source ======================================
// ***************************************************************************************
/// \brief  adds source \f$ \partial_t \phi_1 = S_{\phi} \f$ via Explicit Euler
/// \param  outx	output pointer in x-direction
/// \param  outy	output pointer in y-direction
/// \param  outz	output pointer in z-direction
/// \param  Sx		Source pointer in x-direction
/// \param  Sy		Source pointer in y-direction
/// \param  Sz		Source pointer in z-direction
/// \param	sync	synchronous kernel launching (true, default: false)
// ***************************************************************************************
void ExplicitEulerSource::addSource(Field *outx, Field *outy, Field *outz, Field *Sx, Field *Sy, Field *Sz, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    size_t level = outx->GetLevel();
    size_t bsize = domain->GetSize(level);
    FieldType type = outx->GetType();

    auto d_outx = outx->data;
    auto d_outy = outy->data;
    auto d_outz = outz->data;
    auto d_Sx = Sx->data;
    auto d_Sy = Sy->data;
    auto d_Sz = Sz->data;

    auto dt = m_dt;
    auto dir = m_dir_vel;

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

#pragma acc data present(d_outx[:bsize], d_outy[:bsize], d_outz[:bsize], d_Sx[:bsize], d_Sy[:bsize], d_Sz[:bsize])
    {
        //check directions of source
        //x - direction
        if (dir.find('x') != std::string::npos) {
#pragma acc parallel loop independent present(d_outx[:bsize], d_Sx[:bsize], d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
                d_outx[i] += dt * d_Sx[i];
            }

            boundary->applyBoundary(d_outx, level, type, sync);
        } // end x- direction

        //y - direction
        if (dir.find('y') != std::string::npos) {
#pragma acc parallel loop independent present(d_outy[:bsize], d_Sy[:bsize], d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
                d_outy[i] += dt * d_Sy[i];
            }

            boundary->applyBoundary(d_outy, level, type, sync);
        } // end y- direction

        //z - direction
        if (dir.find('z') != std::string::npos) {
#pragma acc parallel loop independent present(d_outz[:bsize], d_Sz[:bsize], d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
                d_outz[i] += dt * d_Sz[i];
            }

            boundary->applyBoundary(d_outz, level, type, sync);
        } // end z- direction

        if (sync) {
#pragma acc wait
        }
    }//end acc data
}

// ***************************************************************************************
/// \brief  adds source \f$ \partial_t \phi_1 = S_{\phi} \f$ via Explicit Euler
/// \param  out		output pointer
/// \param  S		Source pointer
/// \param	sync	synchronous kernel launching (true, default: false)
// ***************************************************************************************
void ExplicitEulerSource::addSource(Field *out, Field *S, bool sync) {

    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    size_t level = out->GetLevel();
    auto bsize = domain->GetSize(level);
    FieldType type = out->GetType();

    auto d_out = out->data;
    auto d_S = S->data;

    auto dt = m_dt;

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

#pragma acc data present(d_out[:bsize], d_S[:bsize])
    {
#pragma acc parallel loop independent present(d_out[:bsize], d_S[:bsize], d_iList[:bsize_i]) async
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            d_out[i] += dt * d_S[i];
        }

        boundary->applyBoundary(d_out, level, type, sync);

        if (sync) {
#pragma acc wait
        }
    }//end acc data
}
