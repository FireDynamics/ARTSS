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
    field_u(FieldType::U, 0.0, 0, domain.get_size()),
    field_v(FieldType::V, 0.0, 0, domain.get_size()),
    field_w(FieldType::W, 0.0, 0, domain.get_size()),

    field_u0(FieldType::U, 0.0, 0, domain.get_size()),
    field_v0(FieldType::V, 0.0, 0, domain.get_size()),
    field_w0(FieldType::W, 0.0, 0, domain.get_size()),

    field_u_tmp(FieldType::U, 0.0, 0, domain.get_size()),
    field_v_tmp(FieldType::V, 0.0, 0, domain.get_size()),
    field_w_tmp(FieldType::W, 0.0, 0, domain.get_size()),

    // Turbulent diffusivity
    field_nu_t(FieldType::U, 0.0, 0, domain.get_size()),
    field_kappa_t(FieldType::T, 0.0, 0, domain.get_size()),
    field_gamma_t(FieldType::RHO, 0.0, 0, domain.get_size()),

    // Pressure
    field_p(FieldType::P, 0.0, 0, domain.get_size()),
    field_p0(FieldType::P, 0.0, 0, domain.get_size()),
    field_rhs(FieldType::P, 0.0, 0, domain.get_size()),

    // Temperature
    field_T(FieldType::T, 0.0, 0, domain.get_size()),
    field_T0(FieldType::T, 0.0, 0, domain.get_size()),
    field_T_tmp(FieldType::T, 0.0, 0, domain.get_size()),
    field_T_ambient(FieldType::T, 0.0, 0, domain.get_size()),

    // Concentration
    field_concentration(FieldType::RHO, 0.0, 0, domain.get_size()),
    field_concentration0(FieldType::RHO, 0.0, 0, domain.get_size()),
    field_concentration_tmp(FieldType::RHO, 0.0, 0, domain.get_size()),

    // Forces
    field_force_x(FieldType::U, 0.0, 0, domain.get_size()),
    field_force_y(FieldType::V, 0.0, 0, domain.get_size()),
    field_force_z(FieldType::W, 0.0, 0, domain.get_size()),

    // Sources
    field_source_T(FieldType::T, 0.0, 0, domain.get_size()),
    field_source_concentration(FieldType::RHO, 0.0, 0, domain.get_size()),

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
    field_v0.copy_data(field_v);
    field_w0.copy_data(field_w);
    field_u_tmp.copy_data(field_u);
    field_v_tmp.copy_data(field_v);
    field_w_tmp.copy_data(field_w);
    field_p0.copy_data(field_p);
    field_T0.copy_data(field_T);
    field_T_tmp.copy_data(field_T);
    field_concentration0.copy_data(field_concentration);
    field_concentration_tmp.copy_data(field_concentration);

    if (sync) {
#pragma acc wait
    }
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
void FieldController::couple_vector(Field const &a, Field &a0, Field &a_tmp,
        Field const &b, Field &b0, Field &b_tmp,
        Field const &c, Field &c0, Field &c_tmp, bool sync) {
    FieldController::couple_scalar(a, a0, a_tmp, false);
    FieldController::couple_scalar(b, b0, b_tmp, false);
    FieldController::couple_scalar(c, c0, c_tmp, false);

    if (sync) {
#pragma acc wait
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

