/// \file       IPressure.cpp
/// \brief      Interface for pressure method
/// \date       Sep 14, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "IPressure.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"

//======================================== Divergence ====================================
// ***************************************************************************************
/// \brief  calculates divergence \f$ \nabla \cdot u\f$  via central finite differences
/// \param  out   output pointer (\f$ \nabla \cdot u\f$)
/// \param  in_x  input pointer (x -velocity)
/// \param  in_y  input pointer (y -velocity)
/// \param  in_z  input pointer (z -velocity)
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void IPressure::divergence(Field *out, const Field *in_x, const Field *in_y, const Field *in_z, bool sync) {

    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    auto d_out = out->data;
    auto d_inx = in_x->data;
    auto d_iny = in_y->data;
    auto d_inz = in_z->data;

    auto Nx = domain->get_Nx(out->GetLevel());
    auto Ny = domain->get_Ny(out->GetLevel());
    auto dx = domain->get_dx(out->GetLevel());
    auto dy = domain->get_dy(out->GetLevel());
    auto dz = domain->get_dz(out->GetLevel());
    auto rdx = 1. / dx;
    auto rdy = 1. / dy;
    auto rdz = 1. / dz;

    auto size = domain->get_size(out->GetLevel());

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
///     to make \a \f$ \phi_4 \f$ divergence-free
/// \param  out_u output pointer (x -velocity)
/// \param  out_v output pointer (y -velocity)
/// \param  out_w output pointer (z -velocity)
/// \param  in_u  input pointer (x -velocity)
/// \param  in_v  input pointer (y -velocity)
/// \param  in_w  input pointer (z -velocity)
/// \param  in_p  input pointer (pressure)
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void IPressure::projection(Field *out_u, Field *out_v, Field *out_w, const Field *in_u, const Field *in_v, const Field *in_w, const Field *in_p, bool sync) {

    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    auto d_outu = out_u->data;
    auto d_outv = out_v->data;
    auto d_outw = out_w->data;

    auto d_inu = in_u->data;
    auto d_inv = in_v->data;
    auto d_inw = in_w->data;
    auto d_inp = in_p->data;

    auto Nx = domain->get_Nx(out_u->GetLevel());
    auto Ny = domain->get_Ny(out_u->GetLevel());

    auto dx = domain->get_dx(out_u->GetLevel());
    auto dy = domain->get_dy(out_u->GetLevel());
    auto dz = domain->get_dz(out_u->GetLevel());

    auto rdx = 1. / dx;
    auto rdy = 1. / dy;
    auto rdz = 1. / dz;

    auto size = domain->get_size(out_u->GetLevel());

    auto typeu = out_u->GetType();
    auto typev = out_v->GetType();
    auto typew = out_w->GetType();

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
