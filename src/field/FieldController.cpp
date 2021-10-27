/// \file       FieldController.cpp
/// \brief      manages everything regarding any field object
/// \date       Aug 12, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "FieldController.h"
#include <string>
#include "../DomainData.h"
#include "../boundary/BoundaryController.h"

FieldController::FieldController():
        // Variables
        // Velocities
        field_u(FieldType::U, 0),  // initialise with 0 in case they won't be used (e.g. pressure test) to prevent uninitialised warnings in valgrind
        field_v(FieldType::V, 0),
        field_w(FieldType::W, 0),

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
        field_p(FieldType::P, 0),
        field_p0(FieldType::P),
        field_rhs(FieldType::P),

        // Temperature
        field_T(FieldType::T, 0),
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
        sight(FieldType::RHO, 1.0) {

    field_u.copyin();
    field_v.copyin();
    field_w.copyin();

    field_u0.copyin();
    field_v0.copyin();
    field_w0.copyin();

    field_u_tmp.copyin();
    field_v_tmp.copyin();
    field_w_tmp.copyin();

    field_p.copyin();
    field_p0.copyin();
    field_rhs.copyin();

    field_T.copyin();
    field_T_ambient.copyin();

    field_concentration.copyin();

    field_force_x.copyin();
    field_force_y.copyin();
    field_force_z.copyin();

    field_source_T.copyin();
    field_T0.copyin();
    field_T_tmp.copyin();

    field_source_concentration.copyin();
    field_concentration0.copyin();
    field_concentration_tmp.copyin();

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
    boundary->apply_boundary(field_u);
    boundary->apply_boundary(field_v);
    boundary->apply_boundary(field_w);
    boundary->apply_boundary(field_p);
    boundary->apply_boundary(field_T);
    boundary->apply_boundary(field_concentration);

    // TODO necessary?
    boundary->apply_boundary(field_T_ambient);
}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates variables for the next iteration step or time dependent parameters such as temperature source function
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void FieldController::update_data(bool sync) {
    // TODO parallelisable ?
    field_u0.copy_data(field_u);
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
void FieldController::couple_vector(const Field &a, Field &a0, Field &a_tmp,
                                    const Field &b, Field &b0, Field &b_tmp,
                                    const Field &c, Field &c0, Field &c_tmp,
                                    bool sync) {
    // TODO parallelisable ?
    a0.copy_data(a);
    a_tmp.copy_data(a);

    b0.copy_data(b);
    b_tmp.copy_data(b);

    c0.copy_data(c);
    c_tmp.copy_data(c);

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
void FieldController::couple_scalar(const Field &a, Field &a0, Field &a_tmp, bool sync) {
    // TODO parallelisable ?
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
#pragma acc wait
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

