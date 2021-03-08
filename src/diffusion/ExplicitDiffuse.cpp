/// \file       ExplicitDiffuse.cpp
/// \brief      Solves diffusion equation with an explicit method
/// \date       December 12, 2019
/// \author     Max Boehler
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>

#ifdef _OPENACC
#include <accelmath.h>
#endif

#include "ExplicitDiffuse.h"


ExplicitDiffuse::ExplicitDiffuse() :
    m_dt(Parameters::getInstance()->get_real("physical_parameters/dt")),
    m_cs(Parameters::getInstance()->get_real("solver/turbulence/Cs")),
    m_domain(*Domain::getInstance()),
    m_boundary(BoundaryController::getInstance()),
    m_logger(Utility::create_logger("ExplicitDiffuse")) {
}

// ExplicitDiffuse::ExplicitDiffuse(const Domain &domain, const BoundaryController &boundary, real dt, real cs) :
//     m_dt(dt),
//     m_cs(cs),
//     m_domain(domain),
//     m_boundary(&boundary) {
// }

//====================================== Diffuse ===============================================
// ***************************************************************************************
/// \brief  solves diffusion equation \f$ \partial_t \phi_2 = \nu \ nabla^2 \phi_2 \f$
/// \param  out     output pointer
/// \param  in      input pointer
/// \param  b       source pointer
/// \param  D       diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ExplicitDiffuse::diffuse(Field *out, const Field &in, const Field &b,
    const Field &u, const Field &v, const Field &w,
    real D, bool sync) {
    // local variables and parameters for GPU
    FieldType type = out->get_type();

    auto d_out = out->data;

    auto boundary = BoundaryController::getInstance();

    ExplicitStep(out, in, D, sync);
    boundary->applyBoundary(d_out, type, sync);
}

//====================================== Turbulent Diffuse ===============================================
// ***************************************************************************************
/// \brief  solves diffusion equation \f$ \partial_t \phi_2 = \nu \ nabla^2 \phi_2 \f$
/// \param  out     output pointer
/// \param  in      input pointer
/// \param  b       source pointer
/// \param  D       diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ExplicitDiffuse::diffuse(Field *out, const Field &in, const Field &b,
    const Field &u, const Field &v, const Field &w,
    real D, const Field &EV, bool sync) {

    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    auto bsize = domain->get_size(out->get_level());
    FieldType type = out->get_type();

    auto d_out = out->data;

    auto boundary = BoundaryController::getInstance();

    ExplicitStep(out, in, u, v, w, EV, D, m_cs, sync);
    boundary->applyBoundary(d_out, type, sync);
}


void ExplicitDiffuse::ExplicitStep(Field *out, const Field &in, real D, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(out->get_level()); //due to unnecessary parameter passing of *this
    const size_t Ny = domain->get_Ny(out->get_level());

    auto bsize = domain->get_size(out->get_level());

    auto d_out = out->data;
    auto d_in = in.data;

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

    real rdx = D / (domain->get_dx(out->get_level()) * domain->get_dx(out->get_level()));
    real rdy = D / (domain->get_dy(out->get_level()) * domain->get_dy(out->get_level()));
    real rdz = D / (domain->get_dz(out->get_level()) * domain->get_dz(out->get_level()));

#pragma acc parallel loop independent present(d_out[:bsize], d_in[:bsize], d_iList[:bsize_i]) async
    for (size_t ii = 0; ii < bsize_i; ++ii) {
        const size_t idx = d_iList[ii];

        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);

        d_out[idx] = d_in[idx] +
                     m_dt * ((d_in[IX(i + 1, j, k, Nx, Ny)] - 2 * d_in[idx] + d_in[IX(i - 1, j, k, Nx, Ny)]) * rdx
                           + (d_in[IX(i, j + 1, k, Nx, Ny)] - 2 * d_in[idx] + d_in[IX(i, j - 1, k, Nx, Ny)]) * rdy
                           + (d_in[IX(i, j, k + 1, Nx, Ny)] - 2 * d_in[idx] + d_in[IX(i, j, k - 1, Nx, Ny)]) * rdz
                     );
    }

    if (sync) {
#pragma acc wait
    }
}


// Turbulent Diffuse
void ExplicitDiffuse::ExplicitStep(Field *out, const Field &in,
        const Field &u, const Field &v, const Field &w,
        const Field &EV, real D, real C, bool sync) {
    // local variables and parameters for GPU
    const auto level = out->get_level();
    const size_t Nx = m_domain.get_Nx(level);
    const size_t Ny = m_domain.get_Ny(level);
    const auto dx = m_domain.get_dx(level);
    const auto dy = m_domain.get_dy(level);
    const auto dz = m_domain.get_dz(level);

    auto bsize = m_domain.get_size(level);

    auto d_out = out->data;
    auto d_in = in.data;

    size_t *d_iList = m_boundary->get_innerList_level_joined();
    auto bsize_i = m_boundary->getSize_innerList();

    const real rdx = 1 / dx;
    const real rdy = 1 / dy;
    const real rdz = 1 / dz;

    const real rdxx = 1 / (dx * dx);
    const real rdyy = 1 / (dy * dy);
    const real rdzz = 1 / (dz * dz);

    // p. 3 (9)
    const auto delta = cbrt(m_domain.get_dx(level) * m_domain.get_dy(level) * m_domain.get_dz(level));
    const auto pre_factor = D * C * C * delta * delta;  // D * (C_S*\Delta)^2

#pragma acc parallel loop independent present(d_out[:bsize], d_in[:bsize], d_iList[:bsize_i], d_ev[:bsize]) async
    for (size_t ii = 0; ii < bsize_i; ++ii) {
        const size_t idx = d_iList[ii];

        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);

        const size_t ix_im = IX(i-1, j, k, Nx, Ny);
        const size_t ix_jm = IX(i, j-1, k, Nx, Ny);
        const size_t ix_km = IX(i, j, k-1, Nx, Ny);

        const size_t ix_i1 = IX(i+1, j, k, Nx, Ny);
        const size_t ix_j1 = IX(i, j+1, k, Nx, Ny);
        const size_t ix_k1 = IX(i, j, k+1, Nx, Ny);

        const auto u_c = u.data[idx];
        const auto v_c = v.data[idx];
        const auto w_c = w.data[idx];
        const auto rho = EV.data[idx];

        const auto u_x1 = u.data[ix_i1];
        const auto v_x1 = v.data[ix_i1];
        const auto w_x1 = w.data[ix_i1];

        const auto u_y1 = u.data[ix_j1];
        const auto v_y1 = v.data[ix_j1];
        const auto w_y1 = w.data[ix_j1];

        const auto u_z1 = u.data[ix_k1];
        const auto v_z1 = v.data[ix_k1];
        const auto w_z1 = w.data[ix_k1];

        const auto d_u_x = (u_x1 - u_c) * rdx;  // p u / p x (4.34 FDS_TR)
        const auto d_v_x = (v_x1 - v_c) * rdx;
        const auto d_w_x = (w_x1 - w_c) * rdx;

        const auto d_u_y = (u_y1 - u_c) * rdy;
        const auto d_v_y = (v_y1 - v_c) * rdy;
        const auto d_w_y = (w_y1 - w_c) * rdy;

        const auto d_u_z = (u_z1 - u_c) * rdz;
        const auto d_v_z = (v_z1 - v_c) * rdz;
        const auto d_w_z = (w_z1 - w_c) * rdz;

        // (grad(u))^2 (4.33 FDS_TR)
        const auto grad_u_sqr = d_u_x * d_u_x + d_v_y * d_v_y + d_w_z * d_w_z;

        // |S|^2 = (4.33 FDS_TR)
        const auto magSSquare = 2.0 * d_u_x * d_u_x
            + 2.0 * d_v_y * d_v_y
            + 2.0 * d_w_z * d_w_z
            + (d_u_y + d_v_x) * (d_u_y + d_v_x)
            + (d_u_z + d_w_x) * (d_u_z + d_w_x)
            + (d_v_z + d_w_y) * (d_v_y + d_w_y)
            - 2.0 / 3.0 * grad_u_sqr;

        if (magSSquare <= 0.0)
            m_logger->warn("AHHHHH {}", magSSquare);

        const auto magS = sqrt(std::max(magSSquare, 0.0));
        const auto nu = rho * pre_factor * magS;

        // (3.26 FDS_TR)
        d_out[idx] = m_dt * nu *
                        ((d_in[ix_i1] - 2 * d_in[idx] + d_in[ix_im]) * rdxx
                       + (d_in[ix_j1] - 2 * d_in[idx] + d_in[ix_jm]) * rdyy
                       + (d_in[ix_k1] - 2 * d_in[idx] + d_in[ix_km]) * rdzz);
    }

    if (sync) {
#pragma acc wait
    }
}
