/// \file       FieldController.cpp
/// \brief      manages everything reagarding any field object
/// \date       Aug 12, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "FieldController.h"
#include <string>
#include "../Domain.h"
#include "../boundary/BoundaryController.h"
#include "../solver/SolverSelection.h"
#include "../utility/Parameters.h"

FieldController::FieldController(Domain const &domain):
    // Variables
    // Velocities
    field_u(FieldType::U),
    field_v(FieldType::V),
    field_w(FieldType::W),

    field_u0(FieldType::U),
    field_v0(FieldType::V),
    field_w0(FieldType::W),

    field_u_tmp(FieldType::U),
    field_v_tmp(FieldType::V),
    field_w_tmp(FieldType::W),

    // Turbulent diffusivity
    field_nu_t(FieldType::NU),
    field_kappa_t(FieldType::T),
    field_gamma_t(FieldType::RHO),

    // Pressure
    field_p(FieldType::P),
    field_p0(FieldType::P),
    field_rhs(FieldType::P),

    // Temperature
    field_T(FieldType::T),
    field_T0(FieldType::T),
    field_T_tmp(FieldType::T),
    field_T_ambient(FieldType::T),

    // Concentration
    field_concentration(FieldType::RHO),
    field_concentration0(FieldType::RHO),
    field_concentration_tmp(FieldType::RHO),

    // Forces
    field_force_x(FieldType::U),
    field_force_y(FieldType::V),
    field_force_z(FieldType::W),

    // Sources
    field_source_T(FieldType::T),
    field_source_concentration(FieldType::RHO),

    // Fields for sight of boundaries
    sight(FieldType::RHO, 1.0, 0, domain.get_size()) {
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
    boundary->apply_boundary(field_u.data, field_u.get_type());
    boundary->apply_boundary(field_v.data, field_v.get_type());
    boundary->apply_boundary(field_w.data, field_w.get_type());
    boundary->apply_boundary(field_p.data, field_p.get_type());
    boundary->apply_boundary(field_T.data, field_T.get_type());
    boundary->apply_boundary(field_concentration.data, field_concentration.get_type());

    // TODO necessary?
    boundary->apply_boundary(field_T_ambient.data, field_T_ambient.get_type());
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
    auto bsize = Domain::getInstance()->get_size();

    const auto d_u = field_u.data;                        //due to const correctness
    const auto d_v = field_v.data;
    const auto d_w = field_w.data;
    const auto d_u0 = field_u0.data;
    const auto d_v0 = field_v0.data;
    const auto d_w0 = field_w0.data;
    const auto d_u_tmp = field_u_tmp.data;
    const auto d_v_tmp = field_v_tmp.data;
    const auto d_w_tmp = field_w_tmp.data;
    const auto d_p = field_p.data;
    const auto d_p0 = field_p0.data;
    const auto d_T = field_T.data;
    const auto d_T0 = field_T0.data;
    const auto d_T_tmp = field_T_tmp.data;
    const auto d_C = field_concentration.data;
    const auto d_C0 = field_concentration0.data;
    const auto d_C_tmp = field_concentration_tmp.data;

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_inner_list_level_joined();
    size_t bsize_i = boundary->get_size_inner_list();
    size_t *d_bList = boundary->get_boundary_list_level_joined();
    size_t bsize_b = boundary->get_size_boundary_list();
    size_t *d_oList = boundary->get_obstacle_list();
    size_t bsize_o = boundary->get_size_obstacle_list();

#pragma acc data present(d_iList[:bsize_i], d_bList[:bsize_b], d_oList[:bsize_o], d_u[:bsize], d_v[:bsize], d_w[:bsize], \
                         d_u0[:bsize], d_v0[:bsize], d_w0[:bsize], d_u_tmp[:bsize], d_v_tmp[:bsize], d_w_tmp[:bsize], \
                         d_p[:bsize], d_p0[:bsize], d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], d_C[:bsize], d_C0[:bsize], d_C_tmp[:bsize])
    {
        // inner
#pragma acc parallel loop independent present(d_iList[:bsize_i], \
                                              d_u[:bsize], d_v[:bsize], d_w[:bsize], \
                                              d_u0[:bsize], d_v0[:bsize], d_w0[:bsize], \
                                              d_u_tmp[:bsize], d_v_tmp[:bsize], d_w_tmp[:bsize], \
                                              d_p[:bsize], d_p0[:bsize], \
                                              d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], \
                                              d_C[:bsize], d_C0[:bsize], d_C_tmp[:bsize]) async
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t idx = d_iList[j];
            d_u0[idx] = d_u[idx];
            d_v0[idx] = d_v[idx];
            d_w0[idx] = d_w[idx];
            d_u_tmp[idx] = d_u[idx];
            d_v_tmp[idx] = d_v[idx];
            d_w_tmp[idx] = d_w[idx];
            d_p0[idx] = d_p[idx];
            d_T0[idx] = d_T[idx];
            d_T_tmp[idx] = d_T[idx];
            d_C0[idx] = d_C[idx];
            d_C_tmp[idx] = d_C[idx];
        }
        // boundary
#pragma acc parallel loop independent present(d_bList[:bsize_b], d_u[:bsize], d_v[:bsize], d_w[:bsize], d_u0[:bsize], d_v0[:bsize], d_w0[:bsize], d_u_tmp[:bsize], d_v_tmp[:bsize], d_w_tmp[:bsize], d_p[:bsize], d_p0[:bsize], d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], d_C[:bsize], d_C0[:bsize], d_C_tmp[:bsize]) async
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t idx = d_bList[j];
            d_u0[idx] = d_u[idx];
            d_v0[idx] = d_v[idx];
            d_w0[idx] = d_w[idx];
            d_u_tmp[idx] = d_u[idx];
            d_v_tmp[idx] = d_v[idx];
            d_w_tmp[idx] = d_w[idx];
            d_p0[idx] = d_p[idx];
            d_T0[idx] = d_T[idx];
            d_T_tmp[idx] = d_T[idx];
            d_C0[idx] = d_C[idx];
            d_C_tmp[idx] = d_C[idx];
        }
        // obstacles
#pragma acc parallel loop independent present(d_oList[:bsize_o], d_u[:bsize], d_v[:bsize], d_w[:bsize], d_u0[:bsize], d_v0[:bsize], d_w0[:bsize], d_u_tmp[:bsize], d_v_tmp[:bsize], d_w_tmp[:bsize], d_p[:bsize], d_p0[:bsize], d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], d_C[:bsize], d_C0[:bsize], d_C_tmp[:bsize]) async
        for (size_t j = 0; j < bsize_o; ++j) {
            const size_t idx = d_oList[j];
            d_u0[idx] = d_u[idx];
            d_v0[idx] = d_v[idx];
            d_w0[idx] = d_w[idx];
            d_u_tmp[idx] = d_u[idx];
            d_v_tmp[idx] = d_v[idx];
            d_w_tmp[idx] = d_w[idx];
            d_p0[idx] = d_p[idx];
            d_T0[idx] = d_T[idx];
            d_T_tmp[idx] = d_T[idx];
            d_C0[idx] = d_C[idx];
            d_C_tmp[idx] = d_C[idx];
        }

        if (sync) {
#pragma acc wait
        }
    }  // end data region
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
void FieldController::couple_vector(
        const Field &a, Field &a0, Field &a_tmp,
        const Field &b, Field &b0, Field &b_tmp,
        const Field &c, Field &c0, Field &c_tmp, bool sync) {
    // local variables and parameters for GPU
    auto d_a = a.data;
    auto d_a0 = a0.data;
    auto d_a_tmp = a_tmp.data;
    auto d_b = b.data;
    auto d_b0 = b0.data;
    auto d_b_tmp = b_tmp.data;
    auto d_c = c.data;
    auto d_c0 = c0.data;
    auto d_c_tmp = c_tmp.data;

    auto size = Domain::getInstance()->get_size(a0.get_level());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_inner_list_level_joined();
    size_t bsize_i = boundary->get_size_inner_list();
    size_t *d_bList = boundary->get_boundary_list_level_joined();
    size_t bsize_b = boundary->get_size_boundary_list();
    size_t *d_oList = boundary->get_obstacle_list();
    size_t bsize_o = boundary->get_size_obstacle_list();

#pragma acc data present(d_a[:size], d_a0[:size], d_a_tmp[:size], d_b[:size], d_b0[:size], d_b_tmp[:size], \
                         d_c[:size], d_c0[:size], d_c_tmp[:size], d_iList[:bsize_i], d_bList[:bsize_b], d_oList[:bsize_o])
    {
        // inner
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            d_a0[i] = d_a[i];
            d_b0[i] = d_b[i];
            d_c0[i] = d_c[i];
            d_a_tmp[i] = d_a[i];
            d_b_tmp[i] = d_b[i];
            d_c_tmp[i] = d_c[i];
        }
        // boundary
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t i = d_bList[j];
            d_a0[i] = d_a[i];
            d_b0[i] = d_b[i];
            d_c0[i] = d_c[i];
            d_a_tmp[i] = d_a[i];
            d_b_tmp[i] = d_b[i];
            d_c_tmp[i] = d_c[i];
        }
        // obstacles
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_o; ++j) {
            const size_t i = d_oList[j];
            d_a0[i] = d_a[i];
            d_b0[i] = d_b[i];
            d_c0[i] = d_c[i];
            d_a_tmp[i] = d_a[i];
            d_b_tmp[i] = d_b[i];
            d_c_tmp[i] = d_c[i];
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
    a0.copy_data(a);
    a_tmp.copy_data(a);
    if (sync) {
#pragma acc wait
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

