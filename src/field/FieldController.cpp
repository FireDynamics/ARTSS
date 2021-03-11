/// \file       FieldController.cpp
/// \brief      
/// \date       Aug 12, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "FieldController.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"

FieldController::FieldController(Domain const &domain):
    // Variables
    // Velocities
    field_u(Field(FieldType::U, 0.0, 0, domain.get_size())),
    field_v(Field(FieldType::V, 0.0, 0, domain.get_size())),
    field_w(Field(FieldType::W, 0.0, 0, domain.get_size())),

    field_u0(Field(FieldType::U, 0.0, 0, domain.get_size())),
    field_v0(Field(FieldType::V, 0.0, 0, domain.get_size())),
    field_w0(Field(FieldType::W, 0.0, 0, domain.get_size())),

    field_u_tmp(Field(FieldType::U, 0.0, 0, domain.get_size())),
    field_v_tmp(Field(FieldType::V, 0.0, 0, domain.get_size())),
    field_w_tmp(Field(FieldType::W, 0.0, 0, domain.get_size())),

    // Turbulent diffusivity
    field_nu_t(Field(FieldType::U, 0.0, 0, domain.get_size())),
    field_kappa_t(Field(FieldType::T, 0.0, 0, domain.get_size())),
    field_gamma_t(Field(FieldType::RHO, 0.0, 0, domain.get_size())),

    // Pressure
    field_p(Field(FieldType::P, 0.0, 0, domain.get_size())),
    field_p0(Field(FieldType::P, 0.0, 0, domain.get_size())),
    field_rhs(Field(FieldType::P, 0.0, 0, domain.get_size())),

    // Temperature
    field_T(Field(FieldType::T, 0.0, 0, domain.get_size())),
    field_T0(Field(FieldType::T, 0.0, 0, domain.get_size())),
    field_T_tmp(Field(FieldType::T, 0.0, 0, domain.get_size())),
    field_T_ambient(Field(FieldType::T, 300, 0, domain.get_size())),

    // Concentration
    field_concentration(Field(FieldType::RHO, 0.0, 0, domain.get_size())),
    field_concentration0(Field(FieldType::RHO, 0.0, 0, domain.get_size())),
    field_concentration_tmp(Field(FieldType::RHO, 0.0, 0, domain.get_size())),

    // Forces
    field_force_x(Field(FieldType::U, 0.0, 0, domain.get_size())),
    field_force_y(Field(FieldType::V, 0.0, 0, domain.get_size())),
    field_force_z(Field(FieldType::W, 0.0, 0, domain.get_size())),

    // Sources
    field_source_T(Field(FieldType::T, 0.0, 0, domain.get_size())),
    field_source_concentration(Field(FieldType::RHO, 0.0, 0, domain.get_size())),

    // Fields for sight of boundaries
    sight(Field(FieldType::RHO, 1.0, 0, domain.get_size())) {
    field_u.copyin();
    field_v.copyin();
    field_w.copyin();
    field_p.copyin();
    field_rhs.copyin();
    field_T.copyin();
    field_T_ambient.copyin();
    field_concentration.copyin();
    field_force_x.copyin();
    field_force_y.copyin();
    field_force_z.copyin();
    field_source_T.copyin();
    field_source_concentration.copyin();
    field_nu_t.copyin();
    field_kappa_t.copyin();
    field_gamma_t.copyin();
}

// ========================================== Set up boundary =======================================
// ***************************************************************************************
/// \brief  initializes boundary cells
// ***************************************************************************************
void FieldController::set_up_boundary() {
    auto boundary = BoundaryController::getInstance();
    boundary->applyBoundary(field_u.data, field_u.get_type());
    boundary->applyBoundary(field_v.data, field_v.get_type());
    boundary->applyBoundary(field_w.data, field_w.get_type());
    boundary->applyBoundary(field_p.data, field_p.get_type());
    boundary->applyBoundary(field_T.data, field_T.get_type());
    boundary->applyBoundary(field_concentration.data, field_concentration.get_type());

    // TODO necessary?
    boundary->applyBoundary(field_T_ambient.data, field_T_ambient.get_type());
}

void FieldController::set_up_temporary_fields() {
    field_u0.copyin();
    field_v0.copyin();
    field_w0.copyin();
    field_u_tmp.copyin();
    field_v_tmp.copyin();
    field_w_tmp.copyin();
    field_p0.copyin();
    field_T0.copyin();
    field_T_tmp.copyin();
    field_concentration0.copyin();
    field_concentration_tmp.copyin();
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates variables for the next iteration step or time dependent parameters such as temperature source function
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void FieldController::update_data(bool sync) {
    // local variables and parameters for GPU
    auto field_C = field_concentration.data;
    auto field_C0 = field_concentration0.data;
    auto field_C_tmp = field_concentration_tmp.data;

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t bsize_i = boundary->getSize_innerList();
    size_t *d_bList = boundary->get_boundaryList_level_joined();
    size_t bsize_b = boundary->getSize_boundaryList();
    size_t *d_oList = boundary->get_obstacleList();
    size_t bsize_o = boundary->getSize_obstacleList();

#pragma acc data present(d_iList[:bsize_i], d_bList[:bsize_b], d_oList[:bsize_o])
#pragma acc data present(field_u, field_v, field_w, field_u0, field_v0, field_w0, field_u_tmp, field_v_tmp, field_w_tmp, \
                         field_p, field_p0, field_T, field_T0, field_T_tmp, field_C, field_C0, field_C_tmp)
    {
        // inner
#pragma acc parallel loop independent present(d_iList[:bsize_i], \
                                              field_u, field_v, field_w, \
                                              field_u0, field_v0, field_w0, field_u_tmp, field_v_tmp, field_d_w_tmp, \
                                              field_p, field_p0, field_T, field_T0, field_T_tmp, \
                                              field_C, field_C0, field_C_tmp) async
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t idx = d_iList[j];
            field_u0[idx] = field_u[idx];
            field_v0[idx] = field_v[idx];
            field_w0[idx] = field_w[idx];
            field_u_tmp[idx] = field_u[idx];
            field_v_tmp[idx] = field_v[idx];
            field_w_tmp[idx] = field_w[idx];
            field_p0[idx] = field_p[idx];
            field_T0[idx] = field_T[idx];
            field_T_tmp[idx] = field_T[idx];
            field_C0[idx] = field_C[idx];
            field_C_tmp[idx] = field_C[idx];
        }
        // boundary
#pragma acc parallel loop independent present(d_bList[:bsize_b], \
                                              field_u, field_v, field_w, field_u0, field_v0, field_w0, field_u_tmp, field_v_tmp, field_w_tmp, \
                                              field_p, field_p0, field_T, field_T0, field_T_tmp, field_C, field_C0, field_C_tmp) async
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t idx = d_bList[j];
            field_u0[idx] = field_u[idx];
            field_v0[idx] = field_v[idx];
            field_w0[idx] = field_w[idx];
            field_u_tmp[idx] = field_u[idx];
            field_v_tmp[idx] = field_v[idx];
            field_w_tmp[idx] = field_w[idx];
            field_p0[idx] = field_p[idx];
            field_T0[idx] = field_T[idx];
            field_T_tmp[idx] = field_T[idx];
            field_C0[idx] = field_C[idx];
            field_C_tmp[idx] = field_C[idx];
        }
        // obstacles
#pragma acc parallel loop independent present(d_bList[:bsize_b], \
                                              field_u, field_v, field_w, field_u0, field_v0, field_w0, field_u_tmp, field_v_tmp, field_w_tmp, \
                                              field_p, field_p0, field_T, field_T0, field_T_tmp, field_C, field_C0, field_C_tmp) async
        for (size_t j = 0; j < bsize_o; ++j) {
            const size_t idx = d_oList[j];
            field_u0[idx] = field_u[idx];
            field_v0[idx] = field_v[idx];
            field_w0[idx] = field_w[idx];
            field_u_tmp[idx] = field_u[idx];
            field_v_tmp[idx] = field_v[idx];
            field_w_tmp[idx] = field_w[idx];
            field_p0[idx] = field_p[idx];
            field_T0[idx] = field_T[idx];
            field_T_tmp[idx] = field_T[idx];
            field_C0[idx] = field_C[idx];
            field_C_tmp[idx] = field_C[idx];
        }

        if (sync) {
#pragma acc wait
        }
    } //end data region
}

//======================================== Couple velocity ====================================
// ***************************************************************************************
/// \brief  couples vector (sets tmp and zero-th variables to current numerical solution)
/// \param  a   current field in x- direction (const)
/// \param  a0    zero-th field in x- direction
/// \param  a_tmp temporal field in x- direction
/// \param  b   current field in y- direction (const)
/// \param  b0    zero-th field in y- direction
/// \param  b_tmp temporal field in y- direction
/// \param  c   current field in z- direction (const)
/// \param  c0    zero-th field in z- direction
/// \param  c_tmp temporal field in z- direction
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void FieldController::couple_vector(Field const &a, Field &a0, Field &a_tmp, Field const &b, Field &b0, Field &b_tmp, Field const &c, Field &c0, Field &c_tmp, bool sync) {
    // local variables and parameters for GPU
    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t bsize_i = boundary->getSize_innerList();
    size_t *d_bList = boundary->get_boundaryList_level_joined();
    size_t bsize_b = boundary->getSize_boundaryList();
    size_t *d_oList = boundary->get_obstacleList();
    size_t bsize_o = boundary->getSize_obstacleList();

#pragma acc data present(a, a0, a_tmp, b, b0, b_tmp,c, c0, c_tmp, \
                            iList[:bsize_i], bList[:bsize_b], oList[:bsize_o])
    {
        // inner
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            a0[i] = a[i];
            b0[i] = b[i];
            c0[i] = c[i];
            a_tmp[i] = a[i];
            b_tmp[i] = b[i];
            c_tmp[i] = c[i];
        }
        // boundary
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t i = d_bList[j];
            a0[i] = a[i];
            b0[i] = b[i];
            c0[i] = c[i];
            a_tmp[i] = a[i];
            b_tmp[i] = b[i];
            c_tmp[i] = c[i];
        }
        // obstacles
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_o; ++j) {
            const size_t i = d_oList[j];
            a0[i] = a[i];
            b0[i] = b[i];
            c0[i] = c[i];
            a_tmp[i] = a[i];
            b_tmp[i] = b[i];
            c_tmp[i] = c[i];
        }

        if (sync) {
#pragma acc wait
        }
    }
}

//======================================= Couple scalar ==================================
// ***************************************************************************************
/// \brief  couples vector (sets tmp and zero-th variables to current numerical solution)
/// \param  a   current field in x- direction (const)
/// \param  a0    zero-th field in x- direction
/// \param  a_tmp temporal field in x- direction
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void FieldController::couple_scalar(Field const &a, Field &a0, Field &a_tmp, bool sync) {
    // local variables and parameters for GPU
    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t bsize_i = boundary->getSize_innerList();
    size_t *d_bList = boundary->get_boundaryList_level_joined();
    size_t bsize_b = boundary->getSize_boundaryList();
    size_t *d_oList = boundary->get_obstacleList();
    size_t bsize_o = boundary->getSize_obstacleList();

#pragma acc data present(a, a0, a_tmp, iList[:bsize_i], bList[:bsize_b], oList[:bsize_o])
    {
        // inner
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            a0[i] = a[i];
            a_tmp[i] = a[i];
        }
        // boundary
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t i = d_bList[j];
            a0[i] = a[i];
            a_tmp[i] = a[i];
        }
        // obstacles
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_o; ++j) {
            const size_t i = d_oList[j];
            a0[i] = a[i];
            a_tmp[i] = a[i];
        }
        if (sync) {
#pragma acc wait
        }
    }
}

void FieldController::update_device() {
    field_u.update_dev();
    field_v.update_dev();
    field_w.update_dev();
    field_p.update_dev();
    field_rhs.update_dev();
    field_T.update_dev();
    field_T_ambient.update_dev();
    field_concentration.update_dev();
    field_force_x.update_dev();
    field_force_y.update_dev();
    field_force_z.update_dev();
    field_source_T.update_dev();
    field_source_concentration.update_dev();
    field_nu_t.update_dev();
    field_kappa_t.update_dev();
    field_gamma_t.update_dev();
}

void FieldController::update_host(){
    field_u.update_host();
    field_v.update_host();
    field_w.update_host();
    field_p.update_host();
    field_rhs.update_host();
    field_T.update_host();
    field_T_ambient.update_host();
    field_concentration.update_host();
    field_force_x.update_host();
    field_force_y.update_host();
    field_force_z.update_host();
    field_source_T.update_host();
    field_source_concentration.update_host();
    field_nu_t.update_host();
    field_kappa_t.update_host();
    field_gamma_t.update_host();
#pragma acc wait
}
