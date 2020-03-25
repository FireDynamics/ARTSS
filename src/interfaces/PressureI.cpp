/// \file 		PressureI.cpp
/// \brief 		Interface for pressure method
/// \date 		Sep 14, 2016
/// \author 	Severt
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "PressureI.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"

//======================================== Divergence ====================================
// ***************************************************************************************
/// \brief  calculates divergence \f$ \nabla \cdot u\f$  via central finite differences
/// \param  out		output pointer (\f$ \nabla \cdot u\f$)
/// \param  in_x	input pointer (x -velocity)
/// \param  in_y	input pointer (y -velocity)
/// \param  in_z	input pointer (z -velocity)
/// \param  sync	synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void PressureI::Divergence(Field *out, const Field *inx, const Field *iny, const Field *inz, bool sync) {

    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    auto d_out = out->data;
    auto d_inx = inx->data;
    auto d_iny = iny->data;
    auto d_inz = inz->data;

    auto Nx = domain->GetNx(out->GetLevel());
    auto Ny = domain->GetNy(out->GetLevel());
    auto dx = domain->Getdx(out->GetLevel());
    auto dy = domain->Getdy(out->GetLevel());
    auto dz = domain->Getdz(out->GetLevel());
    auto rdx = 1. / dx;
    auto rdy = 1. / dy;
    auto rdz = 1. / dz;

    auto size = domain->GetSize(out->GetLevel());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t *d_bList = boundary->get_boundaryList_level_joined();

    auto bsize_i = boundary->getSize_innerList();
    auto bsize_b = boundary->getSize_boundaryList();

#pragma acc data present(d_out[:size], d_inx[:size], d_iny[:size], d_inz[:size], d_iList[:bsize_i], d_bList[:bsize_b])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            d_out[i] = 0.5 * rdx * (d_inx[i + 1] - d_inx[i - 1]) \
 + 0.5 * rdy * (d_iny[i + Nx] - d_iny[i - Nx]) \
 + 0.5 * rdz * (d_inz[i + Nx * Ny] - d_inz[i - Nx * Ny]);
        }

//boundaries
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t i = d_bList[j];
            d_out[i] = 0.;
        }

        if (sync) {
#pragma acc wait
        }
    }//end data region
}

//======================================== Projection ====================================
// ***************************************************************************************
/// \brief  subtracts pressure gradient (\f$ \phi_4 = \phi_3 - \nabla p\f$) via finite differences
/// 		to make \a \f$ \phi_4 \f$ divergence-free
/// \param  out_u	output pointer (x -velocity)
/// \param  out_v	output pointer (y -velocity)
/// \param  out_w	output pointer (z -velocity)
/// \param  in_u	input pointer (x -velocity)
/// \param  in_v	input pointer (y -velocity)
/// \param  in_w	input pointer (z -velocity)
/// \param  in_p	input pointer (pressure)
/// \param  sync	synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void PressureI::Project(Field *outu, Field *outv, Field *outw, const Field *inu, const Field *inv, const Field *inw, const Field *inp, bool sync) {

    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    auto d_outu = outu->data;
    auto d_outv = outv->data;
    auto d_outw = outw->data;

    auto d_inu = inu->data;
    auto d_inv = inv->data;
    auto d_inw = inw->data;
    auto d_inp = inp->data;

    auto Nx = domain->GetNx(outu->GetLevel());
    auto Ny = domain->GetNy(outu->GetLevel());

    auto dx = domain->Getdx(outu->GetLevel());
    auto dy = domain->Getdy(outu->GetLevel());
    auto dz = domain->Getdz(outu->GetLevel());

    auto rdx = 1. / dx;
    auto rdy = 1. / dy;
    auto rdz = 1. / dz;

    auto size = domain->GetSize(outu->GetLevel());

    auto typeu = outu->GetType();
    auto typev = outv->GetType();
    auto typew = outw->GetType();

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();

    auto bsize_i = boundary->getSize_innerList();

#pragma acc data present(d_outu[:size], d_outv[:size], d_outw[:size], d_inu[:size], d_inv[:size], d_inw[:size], d_inp[:size], d_iList[:bsize_i])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            d_outu[i] = d_inu[i] - 0.5 * rdx * (d_inp[i + 1] - d_inp[i - 1]);
            d_outv[i] = d_inv[i] - 0.5 * rdy * (d_inp[i + Nx] - d_inp[i - Nx]);
            d_outw[i] = d_inw[i] - 0.5 * rdz * (d_inp[i + Nx * Ny] - d_inp[i - Nx * Ny]);
        }

        //boundaries
        boundary->applyBoundary(d_outu, typeu, false);
        boundary->applyBoundary(d_outv, typev, false);
        boundary->applyBoundary(d_outw, typew, false);

        if (sync) {
#pragma acc wait
        }
    }//end data region
}
