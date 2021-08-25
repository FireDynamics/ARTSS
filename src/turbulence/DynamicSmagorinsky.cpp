/// \file         DynamicSmagorinsky.cpp
/// \brief        calculates eddy viscosity based on Dynamic Smagorinsky-Lilly LES model
/// \date         September 8, 2016
/// \author       Suryanarayana Maddu
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>

#ifdef _OPENACC
#include <accelmath.h>
#endif

#include "DynamicSmagorinsky.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"

DynamicSmagorinsky::DynamicSmagorinsky() :
    u_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    v_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    w_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),

    uu(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    vv(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    ww(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    uv(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    uw(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    vw(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),

    uu_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    vv_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    ww_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    uv_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    uw_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    vw_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),

    uf_uf(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    vf_vf(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    wf_wf(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    uf_vf(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    uf_wf(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    vf_wf(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),

    L11(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    L22(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    L33(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    L12(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    L13(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    L23(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),

    S11(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    S22(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    S33(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    S12(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    S13(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    S23(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),

    S11_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    S22_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    S33_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    S12_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    S13_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    S23_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),

    P11(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    P22(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    P33(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    P12(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    P13(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    P23(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),

    P11_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    P22_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    P33_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    P12_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    P13_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    P23_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),

    M11(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    M22(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    M33(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    M12(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    M13(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    M23(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),

    S_bar(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    S_bar_f(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()),
    Cs(FieldType::U, 0.0, 0, Domain::getInstance()->get_size()) {
    auto params = Parameters::getInstance();
    m_nu = params->get_real("physical_parameters/nu");

    // Variables related to Dynamic Smagorinsky
    u_f.copyin();
    v_f.copyin();
    w_f.copyin();

    uu.copyin();
    vv.copyin();
    ww.copyin();
    uv.copyin();
    uw.copyin();
    vw.copyin();

    uf_uf.copyin();
    vf_vf.copyin();
    wf_wf.copyin();
    uf_vf.copyin();
    uf_wf.copyin();
    vf_wf.copyin();

    uu_f.copyin();
    vv_f.copyin();
    ww_f.copyin();
    uv_f.copyin();
    uw_f.copyin();
    vw_f.copyin();

    L11.copyin();
    L22.copyin();
    L33.copyin();
    L12.copyin();
    L13.copyin();
    L23.copyin();

    S11.copyin();
    S22.copyin();
    S33.copyin();
    S12.copyin();
    S13.copyin();
    S23.copyin();

    P11.copyin();
    P22.copyin();
    P33.copyin();
    P12.copyin();
    P13.copyin();
    P23.copyin();

    S11_f.copyin();
    S22_f.copyin();
    S33_f.copyin();
    S12_f.copyin();
    S13_f.copyin();
    S23_f.copyin();

    P11_f.copyin();
    P22_f.copyin();
    P33_f.copyin();
    P12_f.copyin();
    P13_f.copyin();
    P23_f.copyin();

    S_bar.copyin();
    S_bar_f.copyin();

    M11.copyin();
    M22.copyin();
    M33.copyin();
    M12.copyin();
    M13.copyin();
    M23.copyin();

    Cs.copyin();
}

//============================ Calculate turbulent viscosity =============================
// ***************************************************************************************
/// \brief  calculates turbulent viscosity
/// \param  ev            output pointer
/// \param  in_u          input pointer of x-velocity
/// \param  in_v          input pointer of y-velocity
/// \param  in_w          input pointer of z-velocity
/// \param  sync          synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void DynamicSmagorinsky::CalcTurbViscosity(
        Field &ev,
        Field const &in_u, Field const &in_v, Field const &in_w, bool sync) {
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(in_u.get_level());
    const size_t Ny = domain->get_Ny(in_v.get_level());

    const real dx = domain->get_dx(in_u.get_level());
    const real dy = domain->get_dy(in_v.get_level());
    const real dz = domain->get_dz(in_w.get_level());

    const real rdx = 1. / dx;
    const real rdy = 1. / dy;
    const real rdz = 1. / dz;

    const real alpha = 2.0;

    const real delta_s = cbrt(dx * dy * dz);

    real num = 0;
    real den = 0;
    real sum = 0;

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_inner_list_level_joined();
    auto bsize_i = boundary->get_size_inner_list();

// Velocity filter
    ExplicitFiltering(u_f, in_u, sync);
    ExplicitFiltering(v_f, in_v, sync);
    ExplicitFiltering(w_f, in_w, sync);

#pragma acc parallel loop independent present(u_f, v_f, w_f, \
                        uu, vv, ww, \
                        uv, uw, vw, \
                        u_f, v_f, w_f, \
                        uf_uf, vf_vf, wf_wf, \
                        uf_vf, uf_wf, vf_wf, \
                        d_iList[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
        // product of velocities
        uu[i] = u_f[i] * u_f[i];
        vv[i] = v_f[i] * v_f[i];
        ww[i] = w_f[i] * w_f[i];
        uv[i] = u_f[i] * v_f[i];
        uw[i] = u_f[i] * w_f[i];
        vw[i] = v_f[i] * w_f[i];

        // product of filtered velocities
        uf_uf[i] = u_f[i] * u_f[i];
        vf_vf[i] = v_f[i] * v_f[i];
        wf_wf[i] = w_f[i] * w_f[i];
        uf_vf[i] = u_f[i] * v_f[i];
        uf_wf[i] = u_f[i] * w_f[i];
        vf_wf[i] = v_f[i] * w_f[i];
    }

    if (sync) {
#pragma acc wait
    }

// calculation  of the filter of velocity products
    ExplicitFiltering(uu_f, uu, sync);
    ExplicitFiltering(vv_f, vv, sync);
    ExplicitFiltering(ww_f, ww, sync);
    ExplicitFiltering(uv_f, uv, sync);
    ExplicitFiltering(vw_f, vw, sync);
    ExplicitFiltering(ww_f, ww, sync);

#pragma acc parallel loop independent present(  L11, L22, L33, L12, L13, L23, \
                        S11, S22, S33, S12, S13, S23, \
                        S_bar, \
                        P11, P22, P33, P12, P13, P23, \
                        uu_f, vv_f, ww_f, \
                        uv_f, uw_f, vw_f, \
                        uf_uf, vf_vf, wf_wf, \
                        uf_vf, uf_wf, vf_wf, \
                        u_f, v_f, w_f, \
                        d_iList[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];

        // Leonard stress
        L11[i] = uu_f[i] - uf_uf[i];
        L22[i] = vv_f[i] - vf_vf[i];
        L33[i] = ww_f[i] - wf_wf[i];
        L12[i] = uv_f[i] - uf_vf[i];
        L13[i] = uw_f[i] - uf_wf[i];
        L23[i] = vw_f[i] - vf_wf[i];

        // strain tensor
        S11[i] = (u_f[i + 1] - u_f[i - 1]) * 0.5 * rdx;
        S22[i] = (v_f[i + Nx] - v_f[i - Nx]) * 0.5 * rdy;
        S33[i] = (w_f[i + Nx * Ny] - w_f[i - Nx * Ny]) * 0.5 * rdz;

        S12[i] = 0.5 * ((u_f[i + Nx] - u_f[i - Nx]) * 0.5 * rdy \
 + (v_f[i + 1] - v_f[i - 1]) * 0.5 * rdx);

        S13[i] = 0.5 * ((u_f[i + Nx * Ny] - u_f[i - Nx * Ny]) * 0.5 * rdz  \
 + (w_f[i + 1] - w_f[i - 1]) * 0.5 * rdx);

        S23[i] = 0.5 * ((v_f[i + Nx * Ny] - v_f[i - Nx * Ny]) * 0.5 * rdz  \
 + (w_f[i + Nx] - w_f[i - Nx]) * 0.5 * rdy);

        // modulus of strain tensor
        S_bar[i] = sqrt(2. * (S11[i] * S11[i] \
 + S22[i] * S22[i] \
 + S33[i] * S33[i] \
 + 2. * (S12[i] * S12[i]) \
 + 2. * (S13[i] * S13[i]) \
 + 2. * (S23[i] * S23[i])));

        // product of strain modulus and strain tensor
        P11[i] = S_bar[i] * S11[i];
        P22[i] = S_bar[i] * S22[i];
        P33[i] = S_bar[i] * S33[i];
        P12[i] = S_bar[i] * S12[i];
        P13[i] = S_bar[i] * S13[i];
        P23[i] = S_bar[i] * S23[i];
    }

    if (sync) {
#pragma acc wait
    }

    // filtering the strain tensor
    ExplicitFiltering(S11_f, S11, sync);
    ExplicitFiltering(S22_f, S22, sync);
    ExplicitFiltering(S33_f, S33, sync);
    ExplicitFiltering(S12_f, S12, sync);
    ExplicitFiltering(S13_f, S13, sync);
    ExplicitFiltering(S23_f, S23, sync);

    // filtering the product of strain tensor modulus and strain tensor
    ExplicitFiltering(P11_f, P11, sync);
    ExplicitFiltering(P22_f, P22, sync);
    ExplicitFiltering(P33_f, P33, sync);
    ExplicitFiltering(P12_f, P12, sync);
    ExplicitFiltering(P13_f, P13, sync);
    ExplicitFiltering(P23_f, P23, sync);

#pragma acc parallel loop independent present(  S_bar_f, S11_f, S22_f, S33_f, \
                        S12_f, S13_f, S23_f, \
                        M11, M22, M33, \
                        M12, M13, M23, \
                        P11_f, P22_f, P33_f, \
                        P12_f, P13_f, P23_f, \
                        L11, L22, L33, \
                        L12, L13, L23, \
                        Cs, \
                        d_iList[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
        // modulus of filtered strain tensor
        S_bar_f[i] = sqrt(2. * (S11_f[i] * S11_f[i] \
 + S22_f[i] * S22_f[i] \
 + S33_f[i] * S33_f[i] \
 + 2. * (S12_f[i] * S12_f[i]) \
 + 2. * (S13_f[i] * S13_f[i]) \
 + 2. * (S23_f[i] * S23_f[i])));

        // High frequency resolved terms
        M11[i] = 2.0 * delta_s * delta_s * (P11_f[i] - alpha * alpha * S_bar_f[i] * S11_f[i]);
        M22[i] = 2.0 * delta_s * delta_s * (P22_f[i] - alpha * alpha * S_bar_f[i] * S22_f[i]);
        M33[i] = 2.0 * delta_s * delta_s * (P33_f[i] - alpha * alpha * S_bar_f[i] * S33_f[i]);
        M12[i] = 2.0 * delta_s * delta_s * (P12_f[i] - alpha * alpha * S_bar_f[i] * S12_f[i]);
        M13[i] = 2.0 * delta_s * delta_s * (P13_f[i] - alpha * alpha * S_bar_f[i] * S13_f[i]);
        M23[i] = 2.0 * delta_s * delta_s * (P23_f[i] - alpha * alpha * S_bar_f[i] * S23_f[i]);

        num = L11[i] * M11[i] + L22[i] * M22[i] * L33[i] * M33[i]
              + 2.0 * L12[i] * M12[i] + 2.0 * L13[i] * M13[i] + 2.0 * L23[i] * M23[i];
        den = M11[i] * M11[i] + M22[i] * M22[i] * M33[i] * M33[i]
              + 2.0 * M12[i] * M12[i] + 2.0 * M13[i] * M13[i] + 2.0 * M23[i] * M23[i];

        // dynamic constant
        Cs[i] = num / den;
    }

    if (sync) {
#pragma acc wait
    }

    // local averaging of the coefficients
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
        sum = 0;
        for (size_t l = 0; l < 3; l++) {
            for (size_t m = 0; m < 3; m++) {
                for (size_t n = 0; n < 3; n++) {
                    sum = sum + Cs[i - 1 + l + Nx * (m - 1) + Nx * Ny * (n - 1)];
                }
            }
        }
        Cs[i] = sum / 27.0;
        ev[i] = 2.0 * Cs[i] * delta_s * delta_s * S_bar[i];
    }

    // negative coefficients are allowed unless they don't make the effective viscosity zero. In our case d_ev is very small to do that
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
        if ((ev[i] + m_nu) < 0) {
            ev[i] = 0;
        }
    }
}

//============================ Explicit filtering =============================
// ***************************************************************************************
/// \brief  explicitly filters variables
/// \param  out           output pointer
/// \param  in            input pointer
/// \param  sync          synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void DynamicSmagorinsky::ExplicitFiltering(Field &out, Field const &in, bool sync) {
    auto domain = Domain::getInstance();

    const size_t Nx = domain->get_Nx(out.get_level());
    const size_t Ny = domain->get_Ny(out.get_level());
    real sum = 0;

    //Implement a discrete filter by trapezoidal or simpsons rule.
    real a[3] = {1. / 4., 1. / 2., 1. / 4.};  //trapezoidal weights
    //real a[3] = {1./6.,1./3.,1./6.};  //simpsons weights

    //Construction by product combination for trapezoidal
    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_inner_list_level_joined();
    auto bsize_i = boundary->get_size_inner_list();

#pragma acc parallel loop independent present(out, in, a[:3], d_iList[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
        sum = 0;
#pragma acc loop independent collapse(3)
        for (size_t l = 0; l < 3; l++) {
            for (size_t m = 0; m < 3; m++) {
                for (size_t n = 0; n < 3; n++) {
                    sum = sum + a[l] * a[m] * a[n] * in[i - 1 + l + Nx * (m - 1) + Nx * Ny * (n - 1)];
                }
            }
        }
        out[i] = sum;
    }

    if (sync) {
#pragma acc wait
    }
}
