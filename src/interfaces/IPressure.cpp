/// \file       IPressure.cpp
/// \brief      Interface for pressure method
/// \date       Sep 14, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include "IPressure.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"
#include "../visualisation/Visual.h"

//======================================== Divergence ====================================
// ***************************************************************************************
/// \brief  calculates divergence \f$ \nabla \cdot u\f$  via central finite differences
/// \param  out   output pointer (\f$ \nabla \cdot u\f$)
/// \param  in_x  input pointer (x -velocity)
/// \param  in_y  input pointer (y -velocity)
/// \param  in_z  input pointer (z -velocity)
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void IPressure::divergence(Field &out,
        Field const &in_x, Field const &in_y, Field const &in_z, bool sync) {
    auto domain = Domain::getInstance();

    auto Nx = domain->get_Nx(out.get_level());
    auto Ny = domain->get_Ny(out.get_level());
    auto dx = domain->get_dx(out.get_level());
    auto dy = domain->get_dy(out.get_level());
    auto dz = domain->get_dz(out.get_level());
    auto rdx = 1. / dx;
    auto rdy = 1. / dy;
    auto rdz = 1. / dz;

    auto size = domain->get_size(out.get_level());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t *d_bList = boundary->get_boundaryList_level_joined();

    auto bsize_i = boundary->getSize_innerList();
    auto bsize_b = boundary->getSize_boundaryList();

#pragma acc data present(d_iList[:bsize_i], d_bList[:bsize_b])
#pragma acc data present(out, inx, iny, inz)
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            out[i] = 0.5 * rdx * (in_x[i + 1] - in_x[i - 1]) \
                     + 0.5 * rdy * (in_y[i + Nx] - in_y[i - Nx]) \
                     + 0.5 * rdz * (in_z[i + Nx * Ny] - in_z[i - Nx * Ny]);
        }

// boundaries
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t i = d_bList[j];
            out[i] = 0.;
        }

        if (sync) {
#pragma acc wait
        }
    }
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
void IPressure::projection(Field &out_u, Field &out_v, Field &out_w,
        Field const &in_u, Field const &in_v, Field const &in_w,
        Field const &in_p, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    auto Nx = domain->get_Nx(out_u.get_level());
    auto Ny = domain->get_Ny(out_u.get_level());

    auto dx = domain->get_dx(out_u.get_level());
    auto dy = domain->get_dy(out_u.get_level());
    auto dz = domain->get_dz(out_u.get_level());

    auto rdx = 1. / dx;
    auto rdy = 1. / dy;
    auto rdz = 1. / dz;

    auto typeu = out_u.get_type();
    auto typev = out_v.get_type();
    auto typew = out_w.get_type();

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();

    auto bsize_i = boundary->getSize_innerList();

#pragma acc data present(d_iList[:bsize_i])
#pragma acc data present(outu, outv, outw, inu, inv, inw, inp)
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            out_u[i] = in_u[i] - 0.5 * rdx * (in_p[i + 1] - in_p[i - 1]);
            out_v[i] = in_v[i] - 0.5 * rdy * (in_p[i + Nx] - in_p[i - Nx]);
            out_w[i] = in_w[i] - 0.5 * rdz * (in_p[i + Nx * Ny] - in_p[i - Nx * Ny]);
        }

        // boundaries
        boundary->applyBoundary(out_u.data, typeu, false);
        boundary->applyBoundary(out_v.data, typev, false);
        boundary->applyBoundary(out_w.data, typew, false);

        if (sync) {
#pragma acc wait
        }
    }
}
