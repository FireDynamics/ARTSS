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

DynamicSmagorinsky::DynamicSmagorinsky() {

    auto params = Parameters::getInstance();

    m_nu = params->get_real("physical_parameters/nu");

    u_f = new Field(FieldType::U, 0.0);
    v_f = new Field(FieldType::U, 0.0);
    w_f = new Field(FieldType::U, 0.0);

    uu = new Field(FieldType::U, 0.0);
    vv = new Field(FieldType::U, 0.0);
    ww = new Field(FieldType::U, 0.0);
    uv = new Field(FieldType::U, 0.0);
    uw = new Field(FieldType::U, 0.0);
    vw = new Field(FieldType::U, 0.0);

    uu_f = new Field(FieldType::U, 0.0);
    vv_f = new Field(FieldType::U, 0.0);
    ww_f = new Field(FieldType::U, 0.0);
    uv_f = new Field(FieldType::U, 0.0);
    uw_f = new Field(FieldType::U, 0.0);
    vw_f = new Field(FieldType::U, 0.0);

    uf_uf = new Field(FieldType::U, 0.0);
    vf_vf = new Field(FieldType::U, 0.0);
    wf_wf = new Field(FieldType::U, 0.0);
    uf_vf = new Field(FieldType::U, 0.0);
    uf_wf = new Field(FieldType::U, 0.0);
    vf_wf = new Field(FieldType::U, 0.0);

    L11 = new Field(FieldType::U, 0.0);
    L22 = new Field(FieldType::U, 0.0);
    L33 = new Field(FieldType::U, 0.0);
    L12 = new Field(FieldType::U, 0.0);
    L13 = new Field(FieldType::U, 0.0);
    L23 = new Field(FieldType::U, 0.0);

    S11 = new Field(FieldType::U, 0.0);
    S22 = new Field(FieldType::U, 0.0);
    S33 = new Field(FieldType::U, 0.0);
    S12 = new Field(FieldType::U, 0.0);
    S13 = new Field(FieldType::U, 0.0);
    S23 = new Field(FieldType::U, 0.0);

    S11_f = new Field(FieldType::U, 0.0);
    S22_f = new Field(FieldType::U, 0.0);
    S33_f = new Field(FieldType::U, 0.0);
    S12_f = new Field(FieldType::U, 0.0);
    S13_f = new Field(FieldType::U, 0.0);
    S23_f = new Field(FieldType::U, 0.0);

    P11 = new Field(FieldType::U, 0.0);
    P22 = new Field(FieldType::U, 0.0);
    P33 = new Field(FieldType::U, 0.0);
    P12 = new Field(FieldType::U, 0.0);
    P13 = new Field(FieldType::U, 0.0);
    P23 = new Field(FieldType::U, 0.0);

    P11_f = new Field(FieldType::U, 0.0);
    P22_f = new Field(FieldType::U, 0.0);
    P33_f = new Field(FieldType::U, 0.0);
    P12_f = new Field(FieldType::U, 0.0);
    P13_f = new Field(FieldType::U, 0.0);
    P23_f = new Field(FieldType::U, 0.0);

    M11 = new Field(FieldType::U, 0.0);
    M22 = new Field(FieldType::U, 0.0);
    M33 = new Field(FieldType::U, 0.0);
    M12 = new Field(FieldType::U, 0.0);
    M13 = new Field(FieldType::U, 0.0);
    M23 = new Field(FieldType::U, 0.0);

    S_bar = new Field(FieldType::U, 0.0);
    S_bar_f = new Field(FieldType::U, 0.0);

    Cs = new Field(FieldType::U, 0.0);

    // Variables related to Dynamic Smagorinsky
    auto d_u_f = u_f->data;
    auto d_v_f = v_f->data;
    auto d_w_f = w_f->data;

    auto d_uu = uu->data;
    auto d_vv = vv->data;
    auto d_ww = ww->data;
    auto d_uv = uv->data;
    auto d_uw = uw->data;
    auto d_vw = vw->data;

    auto d_uf_uf = uf_uf->data;
    auto d_vf_vf = vf_vf->data;
    auto d_wf_wf = wf_wf->data;
    auto d_uf_vf = uf_vf->data;
    auto d_uf_wf = uf_wf->data;
    auto d_vf_wf = vf_wf->data;

    auto d_uu_f = uu_f->data;
    auto d_vv_f = vv_f->data;
    auto d_ww_f = ww_f->data;
    auto d_uv_f = uv_f->data;
    auto d_uw_f = uw_f->data;
    auto d_vw_f = vw_f->data;

    auto d_L11 = L11->data;
    auto d_L22 = L22->data;
    auto d_L33 = L33->data;
    auto d_L12 = L12->data;
    auto d_L13 = L13->data;
    auto d_L23 = L23->data;

    auto d_S11 = S11->data;
    auto d_S22 = S22->data;
    auto d_S33 = S33->data;
    auto d_S12 = S12->data;
    auto d_S13 = S13->data;
    auto d_S23 = S23->data;

    auto d_P11 = P11->data;
    auto d_P22 = P22->data;
    auto d_P33 = P33->data;
    auto d_P12 = P12->data;
    auto d_P13 = P13->data;
    auto d_P23 = P23->data;

    auto d_S11_f = S11_f->data;
    auto d_S22_f = S22_f->data;
    auto d_S33_f = S33_f->data;
    auto d_S12_f = S12_f->data;
    auto d_S13_f = S13_f->data;
    auto d_S23_f = S23_f->data;

    auto d_P11_f = P11_f->data;
    auto d_P22_f = P22_f->data;
    auto d_P33_f = P33_f->data;
    auto d_P12_f = P12_f->data;
    auto d_P13_f = P13_f->data;
    auto d_P23_f = P23_f->data;

    auto d_S_bar = S_bar->data;
    auto d_S_bar_f = S_bar_f->data;

    auto d_M11 = M11->data;
    auto d_M22 = M22->data;
    auto d_M33 = M33->data;
    auto d_M12 = M12->data;
    auto d_M13 = M13->data;
    auto d_M23 = M23->data;

    auto d_Cs = Cs->data;

    size_t bsize = Domain::getInstance()->get_size(u_f->GetLevel());

#pragma acc enter data copyin(d_u_f[:bsize], d_v_f[:bsize], d_w_f[:bsize])
#pragma acc enter data copyin(d_uu[:bsize], d_vv[:bsize], d_ww[:bsize], d_uv[:bsize], d_uw[:bsize], d_vw[:bsize])
#pragma acc enter data copyin(d_uu_f[:bsize], d_vv_f[:bsize], d_ww_f[:bsize], d_uv_f[:bsize], d_uw_f[:bsize], d_vw_f[:bsize])
#pragma acc enter data copyin(d_uf_uf[:bsize], d_vf_vf[:bsize], d_wf_wf[:bsize], d_uf_vf[:bsize], d_uf_wf[:bsize], d_vf_wf[:bsize])

#pragma acc enter data copyin(d_L11[:bsize], d_L22[:bsize], d_L33[:bsize], d_L12[:bsize], d_L13[:bsize], d_L23[:bsize])
#pragma acc enter data copyin(d_S11[:bsize], d_S22[:bsize], d_S33[:bsize], d_S12[:bsize], d_S13[:bsize], d_S23[:bsize])
#pragma acc enter data copyin(d_S11_f[:bsize], d_S22_f[:bsize], d_S33_f[:bsize], d_S12_f[:bsize], d_S13_f[:bsize], d_S23_f[:bsize])
#pragma acc enter data copyin(d_P11[:bsize], d_P22[:bsize], d_P33[:bsize], d_P12[:bsize], d_P13[:bsize], d_P23[:bsize])
#pragma acc enter data copyin(d_P11_f[:bsize], d_P22_f[:bsize], d_P33_f[:bsize], d_P12_f[:bsize], d_P13_f[:bsize], d_P23_f[:bsize])
#pragma acc enter data copyin(d_M11[:bsize], d_M22[:bsize], d_M33[:bsize], d_M12[:bsize], d_M13[:bsize], d_M23[:bsize])

#pragma acc enter data copyin(d_S_bar[:bsize], d_S_bar_f[:bsize])

#pragma acc enter data copyin(d_Cs[:bsize])
}

DynamicSmagorinsky::~DynamicSmagorinsky() {

    auto d_u_f = u_f->data;
    auto d_v_f = v_f->data;
    auto d_w_f = w_f->data;

    auto d_uu = uu->data;
    auto d_vv = vv->data;
    auto d_ww = ww->data;
    auto d_uv = uv->data;
    auto d_uw = uw->data;
    auto d_vw = vw->data;

    auto d_uf_uf = uf_uf->data;
    auto d_vf_vf = vf_vf->data;
    auto d_wf_wf = wf_wf->data;
    auto d_uf_vf = uf_vf->data;
    auto d_uf_wf = uf_wf->data;
    auto d_vf_wf = vf_wf->data;

    auto d_uu_f = uu_f->data;
    auto d_vv_f = vv_f->data;
    auto d_ww_f = ww_f->data;
    auto d_uv_f = uv_f->data;
    auto d_uw_f = uw_f->data;
    auto d_vw_f = vw_f->data;

    auto d_L11 = L11->data;
    auto d_L22 = L22->data;
    auto d_L33 = L33->data;
    auto d_L12 = L12->data;
    auto d_L13 = L13->data;
    auto d_L23 = L23->data;

    auto d_S11 = S11->data;
    auto d_S22 = S22->data;
    auto d_S33 = S33->data;
    auto d_S12 = S12->data;
    auto d_S13 = S13->data;
    auto d_S23 = S23->data;

    auto d_P11 = P11->data;
    auto d_P22 = P22->data;
    auto d_P33 = P33->data;
    auto d_P12 = P12->data;
    auto d_P13 = P13->data;
    auto d_P23 = P23->data;

    auto d_S11_f = S11_f->data;
    auto d_S22_f = S22_f->data;
    auto d_S33_f = S33_f->data;
    auto d_S12_f = S12_f->data;
    auto d_S13_f = S13_f->data;
    auto d_S23_f = S23_f->data;

    auto d_P11_f = P11_f->data;
    auto d_P22_f = P22_f->data;
    auto d_P33_f = P33_f->data;
    auto d_P12_f = P12_f->data;
    auto d_P13_f = P13_f->data;
    auto d_P23_f = P23_f->data;

    auto d_S_bar = S_bar->data;
    auto d_S_bar_f = S_bar_f->data;

    auto d_M11 = M11->data;
    auto d_M22 = M22->data;
    auto d_M33 = M33->data;
    auto d_M12 = M12->data;
    auto d_M13 = M13->data;
    auto d_M23 = M23->data;

    auto d_Cs = Cs->data;

    size_t bsize = Domain::getInstance()->get_size(u_f->GetLevel());

#pragma acc exit data delete(d_u_f[:bsize], d_v_f[:bsize], d_w_f[:bsize])
    delete u_f;
    delete v_f;
    delete w_f;

#pragma acc exit data delete(d_uu[:bsize], d_vv[:bsize], d_ww[:bsize], d_uv[:bsize], d_uw[:bsize], d_vw[:bsize])
    delete uu;
    delete vv;
    delete ww;
    delete uv;
    delete uw;
    delete vw;

#pragma acc exit data delete(d_uu_f[:bsize], d_vv_f[:bsize], d_ww_f[:bsize], d_uv_f[:bsize], d_uw_f[:bsize], d_vw_f[:bsize])
    delete uu_f;
    delete vv_f;
    delete ww_f;
    delete uv_f;
    delete uw_f;
    delete vw_f;

#pragma acc exit data delete(d_uf_uf[:bsize], d_vf_vf[:bsize], d_wf_wf[:bsize], d_uf_vf[:bsize], d_uf_wf[:bsize], d_vf_wf[:bsize])
    delete uf_uf;
    delete vf_vf;
    delete wf_wf;
    delete uf_vf;
    delete uf_wf;
    delete vf_wf;

#pragma acc exit data delete(d_L11[:bsize], d_L22[:bsize], d_L33[:bsize], d_L12[:bsize], d_L13[:bsize], d_L23[:bsize])
    delete L11;
    delete L22;
    delete L33;
    delete L12;
    delete L13;
    delete L23;

#pragma acc exit data delete(d_S11[:bsize], d_S22[:bsize], d_S33[:bsize], d_S12[:bsize], d_S13[:bsize], d_S23[:bsize])
    delete S11;
    delete S22;
    delete S33;
    delete S12;
    delete S13;
    delete S23;

#pragma acc exit data delete(d_S11_f[:bsize], d_S22_f[:bsize], d_S33_f[:bsize], d_S12_f[:bsize], d_S13_f[:bsize], d_S23_f[:bsize])
    delete S11_f;
    delete S22_f;
    delete S33_f;
    delete S12_f;
    delete S13_f;
    delete S23_f;

#pragma acc exit data delete(d_P11[:bsize], d_P22[:bsize], d_P33[:bsize], d_P12[:bsize], d_P13[:bsize], d_P23[:bsize])
    delete P11;
    delete P22;
    delete P33;
    delete P12;
    delete P13;
    delete P23;

#pragma acc exit data delete(d_P11_f[:bsize], d_P22_f[:bsize], d_P33_f[:bsize], d_P12_f[:bsize], d_P13_f[:bsize], d_P23_f[:bsize])
    delete P11_f;
    delete P22_f;
    delete P33_f;
    delete P12_f;
    delete P13_f;
    delete P23_f;

#pragma acc exit data delete(d_M11[:bsize], d_M22[:bsize], d_M33[:bsize], d_M12[:bsize], d_M13[:bsize], d_M23[:bsize])
    delete M11;
    delete M22;
    delete M33;
    delete M12;
    delete M13;
    delete M23;

#pragma acc exit data delete(d_S_bar[:bsize], d_S_bar_f[:bsize])
    delete S_bar;
    delete S_bar_f;

#pragma acc exit data delete(d_Cs[:bsize])
    delete Cs;
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
void DynamicSmagorinsky::CalcTurbViscosity(Field *ev, Field *in_u, Field *in_v, Field *in_w, bool sync) {

    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    auto d_u = in_u->data;
    auto d_v = in_v->data;
    auto d_w = in_w->data;
    auto d_ev = ev->data;

    // Variables related to Dynamic Smagorinsky
    // filtered velocities
    auto d_u_f = u_f->data;
    auto d_v_f = v_f->data;
    auto d_w_f = w_f->data;

    // velocity products
    auto d_uu = uu->data;
    auto d_vv = vv->data;
    auto d_ww = ww->data;
    auto d_uv = uv->data;
    auto d_uw = uw->data;
    auto d_vw = vw->data;

    // product of the filtered velocities
    auto d_uu_f = uu_f->data;
    auto d_vv_f = vv_f->data;
    auto d_ww_f = ww_f->data;
    auto d_uv_f = uv_f->data;
    auto d_uw_f = uw_f->data;
    auto d_vw_f = vw_f->data;

    // filters of the velocity products
    auto d_uf_uf = uf_uf->data;
    auto d_vf_vf = vf_vf->data;
    auto d_wf_wf = wf_wf->data;
    auto d_uf_vf = uf_vf->data;
    auto d_uf_wf = uf_wf->data;
    auto d_vf_wf = vf_wf->data;

    // Leonard stress
    auto d_L11 = L11->data;
    auto d_L22 = L22->data;
    auto d_L33 = L33->data;
    auto d_L12 = L12->data;
    auto d_L13 = L13->data;
    auto d_L23 = L23->data;

    // strain tensor
    auto d_S11 = S11->data;
    auto d_S22 = S22->data;
    auto d_S33 = S33->data;
    auto d_S12 = S12->data;
    auto d_S13 = S13->data;
    auto d_S23 = S23->data;

    // second filtered strain tensor
    auto d_S11_f = S11_f->data;
    auto d_S22_f = S22_f->data;
    auto d_S33_f = S33_f->data;
    auto d_S12_f = S12_f->data;
    auto d_S13_f = S13_f->data;
    auto d_S23_f = S23_f->data;

    // Product of strain modulus and strain tensor
    auto d_P11 = P11->data;
    auto d_P22 = P22->data;
    auto d_P33 = P33->data;
    auto d_P12 = P12->data;
    auto d_P13 = P13->data;
    auto d_P23 = P23->data;

    // second filtered product of strain modulus and strain tensor
    auto d_P11_f = P11_f->data;
    auto d_P22_f = P22_f->data;
    auto d_P33_f = P33_f->data;
    auto d_P12_f = P12_f->data;
    auto d_P13_f = P13_f->data;
    auto d_P23_f = P23_f->data;

    // High frequency resolved terms
    auto d_M11 = M11->data;
    auto d_M22 = M22->data;
    auto d_M33 = M33->data;
    auto d_M12 = M12->data;
    auto d_M13 = M13->data;
    auto d_M23 = M23->data;

    // modulus of strain tensor and filtered strain tensor
    auto d_S_bar = S_bar->data;
    auto d_S_bar_f = S_bar_f->data;

    // dynamic constant
    auto d_Cs = Cs->data;

    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(in_u->GetLevel());
    const size_t Ny = domain->get_Ny(in_v->GetLevel());
    const size_t Nz = domain->get_Nz(in_w->GetLevel());

    const real dx = domain->get_dx(in_u->GetLevel());
    const real dy = domain->get_dy(in_v->GetLevel());
    const real dz = domain->get_dz(in_w->GetLevel());

    size_t bsize = domain->get_size(in_u->GetLevel());

    const real rdx = 1. / dx;
    const real rdy = 1. / dy;
    const real rdz = 1. / dz;

    const real alpha = 2.0;

    const real delta_s = cbrt(dx * dy * dz);

    real num = 0;
    real den = 0;
    real sum = 0;

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

// Velocity filter
    ExplicitFiltering(u_f, in_u, sync);
    ExplicitFiltering(v_f, in_v, sync);
    ExplicitFiltering(w_f, in_w, sync);

#pragma acc parallel loop independent present(  d_u[:bsize], d_v[:bsize], d_w[:bsize], \
                        d_uu[:bsize], d_vv[:bsize], d_ww[:bsize], \
                        d_uv[:bsize], d_uw[:bsize], d_vw[:bsize], \
                        d_u_f[:bsize], d_v_f[:bsize], d_w_f[:bsize], \
                        d_uf_uf[:bsize], d_vf_vf[:bsize], d_wf_wf[:bsize], \
                        d_uf_vf[:bsize], d_uf_wf[:bsize], d_vf_wf[:bsize], \
                        d_iList[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
        // product of velocities
        d_uu[i] = d_u[i] * d_u[i];
        d_vv[i] = d_v[i] * d_v[i];
        d_ww[i] = d_w[i] * d_w[i];
        d_uv[i] = d_u[i] * d_v[i];
        d_uw[i] = d_u[i] * d_w[i];
        d_vw[i] = d_v[i] * d_w[i];

        // product of filtered velocities
        d_uf_uf[i] = d_u_f[i] * d_u_f[i];
        d_vf_vf[i] = d_v_f[i] * d_v_f[i];
        d_wf_wf[i] = d_w_f[i] * d_w_f[i];
        d_uf_vf[i] = d_u_f[i] * d_v_f[i];
        d_uf_wf[i] = d_u_f[i] * d_w_f[i];
        d_vf_wf[i] = d_v_f[i] * d_w_f[i];
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

#pragma acc parallel loop independent present(  d_L11[:bsize], d_L22[:bsize], d_L33[:bsize], d_L12[:bsize], d_L13[:bsize], d_L23[:bsize], \
                        d_S11[:bsize], d_S22[:bsize], d_S33[:bsize], d_S12[:bsize], d_S13[:bsize], d_S23[:bsize], \
                        d_S_bar[:bsize], \
                        d_P11[:bsize], d_P22[:bsize], d_P33[:bsize], d_P12[:bsize], d_P13[:bsize], d_P23[:bsize], \
                        d_uu_f[:bsize], d_vv_f[:bsize], d_ww_f[:bsize], \
                        d_uv_f[:bsize], d_uw_f[:bsize], d_vw_f[:bsize], \
                        d_uf_uf[:bsize], d_vf_vf[:bsize], d_wf_wf[:bsize], \
                        d_uf_vf[:bsize], d_uf_wf[:bsize], d_vf_wf[:bsize], \
                        d_u[:bsize], d_v[:bsize], d_w[:bsize], d_iList[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];

        // Leonard stress
        d_L11[i] = d_uu_f[i] - d_uf_uf[i];
        d_L22[i] = d_vv_f[i] - d_vf_vf[i];
        d_L33[i] = d_ww_f[i] - d_wf_wf[i];
        d_L12[i] = d_uv_f[i] - d_uf_vf[i];
        d_L13[i] = d_uw_f[i] - d_uf_wf[i];
        d_L23[i] = d_vw_f[i] - d_vf_wf[i];

        // strain tensor
        d_S11[i] = (d_u[i + 1] - d_u[i - 1]) * 0.5 * rdx;
        d_S22[i] = (d_v[i + Nx] - d_v[i - Nx]) * 0.5 * rdy;
        d_S33[i] = (d_w[i + Nx * Ny] - d_w[i - Nx * Ny]) * 0.5 * rdz;

        d_S12[i] = 0.5 * ((d_u[i + Nx] - d_u[i - Nx]) * 0.5 * rdy \
 + (d_v[i + 1] - d_v[i - 1]) * 0.5 * rdx);

        d_S13[i] = 0.5 * ((d_u[i + Nx * Ny] - d_u[i - Nx * Ny]) * 0.5 * rdz  \
 + (d_w[i + 1] - d_w[i - 1]) * 0.5 * rdx);

        d_S23[i] = 0.5 * ((d_v[i + Nx * Ny] - d_v[i - Nx * Ny]) * 0.5 * rdz  \
 + (d_w[i + Nx] - d_w[i - Nx]) * 0.5 * rdy);

        // modulus of strain tensor
        d_S_bar[i] = sqrt(2. * (d_S11[i] * d_S11[i] \
 + d_S22[i] * d_S22[i] \
 + d_S33[i] * d_S33[i] \
 + 2. * (d_S12[i] * d_S12[i]) \
 + 2. * (d_S13[i] * d_S13[i]) \
 + 2. * (d_S23[i] * d_S23[i])));

        // product of strain modulus and strain tensor
        d_P11[i] = d_S_bar[i] * d_S11[i];
        d_P22[i] = d_S_bar[i] * d_S22[i];
        d_P33[i] = d_S_bar[i] * d_S33[i];
        d_P12[i] = d_S_bar[i] * d_S12[i];
        d_P13[i] = d_S_bar[i] * d_S13[i];
        d_P23[i] = d_S_bar[i] * d_S23[i];
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

#pragma acc parallel loop independent present(  d_S_bar_f[:bsize], d_S11_f[:bsize], d_S22_f[:bsize], d_S33_f[:bsize], \
                        d_S12_f[:bsize], d_S13_f[:bsize], d_S23_f[:bsize], \
                        d_M11[:bsize], d_M22[:bsize], d_M33[:bsize], \
                        d_M12[:bsize], d_M13[:bsize], d_M23[:bsize], \
                        d_P11_f[:bsize], d_P22_f[:bsize], d_P33_f[:bsize], \
                        d_P12_f[:bsize], d_P13_f[:bsize], d_P23_f[:bsize], \
                        d_L11[:bsize], d_L22[:bsize], d_L33[:bsize], \
                        d_L12[:bsize], d_L13[:bsize], d_L23[:bsize], \
                        d_Cs[:bsize], d_iList[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
        // modulus of filtered strain tensor
        d_S_bar_f[i] = sqrt(2. * (d_S11_f[i] * d_S11_f[i] \
 + d_S22_f[i] * d_S22_f[i] \
 + d_S33_f[i] * d_S33_f[i] \
 + 2. * (d_S12_f[i] * d_S12_f[i]) \
 + 2. * (d_S13_f[i] * d_S13_f[i]) \
 + 2. * (d_S23_f[i] * d_S23_f[i])));

        // High frequency resolved terms
        d_M11[i] = 2.0 * delta_s * delta_s * (d_P11_f[i] - alpha * alpha * d_S_bar_f[i] * d_S11_f[i]);
        d_M22[i] = 2.0 * delta_s * delta_s * (d_P22_f[i] - alpha * alpha * d_S_bar_f[i] * d_S22_f[i]);
        d_M33[i] = 2.0 * delta_s * delta_s * (d_P33_f[i] - alpha * alpha * d_S_bar_f[i] * d_S33_f[i]);
        d_M12[i] = 2.0 * delta_s * delta_s * (d_P12_f[i] - alpha * alpha * d_S_bar_f[i] * d_S12_f[i]);
        d_M13[i] = 2.0 * delta_s * delta_s * (d_P13_f[i] - alpha * alpha * d_S_bar_f[i] * d_S13_f[i]);
        d_M23[i] = 2.0 * delta_s * delta_s * (d_P23_f[i] - alpha * alpha * d_S_bar_f[i] * d_S23_f[i]);

        num = d_L11[i] * d_M11[i] + d_L22[i] * d_M22[i] * d_L33[i] * d_M33[i]
              + 2.0 * d_L12[i] * d_M12[i] + 2.0 * d_L13[i] * d_M13[i] + 2.0 * d_L23[i] * d_M23[i];
        den = d_M11[i] * d_M11[i] + d_M22[i] * d_M22[i] * d_M33[i] * d_M33[i]
              + 2.0 * d_M12[i] * d_M12[i] + 2.0 * d_M13[i] * d_M13[i] + 2.0 * d_M23[i] * d_M23[i];

        // dynamic constant
        d_Cs[i] = num / den;
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
                    sum = sum + d_Cs[i - 1 + l + Nx * (m - 1) + Nx * Ny * (n - 1)];
                }
            }
        }
        d_Cs[i] = sum / 27.0;
        d_ev[i] = 2.0 * d_Cs[i] * delta_s * delta_s * d_S_bar[i];
    }

    // negative coefficients are allowed unless they don't make the effective viscosity zero. In our case d_ev is very small to do that
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
        if ((d_ev[i] + m_nu) < 0) {
            d_ev[i] = 0;
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
void DynamicSmagorinsky::ExplicitFiltering(Field *out, const Field *in, bool sync) {

    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    auto d_out = out->data;
    auto d_in = in->data;

    const size_t Nx = domain->get_Nx(out->GetLevel());
    const size_t Ny = domain->get_Ny(out->GetLevel());
    const size_t bsize = domain->get_size(out->GetLevel());
    real sum = 0;

    //Implement a discrete filter by trapezoidal or simpsons rule.
    real a[3] = {1. / 4., 1. / 2., 1. / 4.};  //trapezoidal weights
    //real a[3] = {1./6.,1./3.,1./6.};  //simpsons weights

    //Construction by product combination for trapezoidal
    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

#pragma acc parallel loop independent present(d_out[:bsize], d_in[:bsize], a[:3], d_iList[:bsize_i]) async
    for (size_t j = 0; j < bsize_i; ++j) {
        const size_t i = d_iList[j];
        sum = 0;
#pragma acc loop independent collapse(3)
        for (size_t l = 0; l < 3; l++) {
            for (size_t m = 0; m < 3; m++) {
                for (size_t n = 0; n < 3; n++) {
                    sum = sum + a[l] * a[m] * a[n] * d_in[i - 1 + l + Nx * (m - 1) + Nx * Ny * (n - 1)];
                }
            }
        }
        d_out[i] = sum;
    }

    if (sync) {
#pragma acc wait
    }
}
