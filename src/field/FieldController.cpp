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

FieldController::FieldController() {
    // Variables
    // Velocities
    field_u = new Field(FieldType::U, 0.0);
    field_v = new Field(FieldType::V, 0.0);
    field_w = new Field(FieldType::W, 0.0);

    // Turbulent diffusivity
    field_nu_t = new Field(FieldType::U, 0.0);
    field_kappa_t = new Field(FieldType::T, 0.0);
    field_gamma_t = new Field(FieldType::RHO, 0.0);

    // Pressure
    field_p = new Field(FieldType::P, 0.0);
    field_rhs = new Field(FieldType::P, 0.0);

    // Temperature
    field_T = new Field(FieldType::T, 0.0);
    field_T_ambient = new Field(FieldType::T, 0);

    // Concentration
    field_concentration = new Field(FieldType::RHO, 0.0);

    // Forces
    field_force_x = new Field(FieldType::U, 0.0);
    field_force_y = new Field(FieldType::V, 0.0);
    field_force_z = new Field(FieldType::W, 0.0);

    // Sources
    field_source_T = new Field(FieldType::T, 0.0);
    field_source_concentration = new Field(FieldType::RHO, 0.0);

    // Fields for sight of boundaries
    sight = new Field(FieldType::RHO, 1.0);

    auto d_u = field_u->data;
    auto d_v = field_v->data;
    auto d_w = field_w->data;
    auto d_p = field_p->data;
    auto d_rhs = field_rhs->data;
    auto d_T = field_T->data;
    auto d_T_a = field_T_ambient->data;
    auto d_C = field_concentration->data;
    auto d_f_x = field_force_x->data;
    auto d_f_y = field_force_y->data;
    auto d_f_z = field_force_z->data;
    auto d_S_T = field_source_T->data;
    auto d_S_C = field_source_concentration->data;
    auto d_nu_t = field_nu_t->data;
    auto d_kappa_t = field_kappa_t->data;
    auto d_gamma_t = field_gamma_t->data;

    Domain *domain = Domain::getInstance();
    auto bsize = domain->get_size();
    // copyin fields
#pragma acc enter data copyin(  d_u[:bsize], \
                                d_v[:bsize], \
                                d_w[:bsize], \
                                d_p[:bsize], d_rhs[:bsize], \
                                d_T[:bsize], d_T_a[:bsize], \
                                d_C[:bsize], \
                                d_f_x[:bsize], d_f_y[:bsize], d_f_z[:bsize], d_S_T[:bsize], d_S_C[:bsize], \
                                d_nu_t[:bsize], d_kappa_t[:bsize], d_gamma_t[:bsize])

}

// ==================================== Destructor ====================================
// ***************************************************************************************
FieldController::~FieldController() {
    Domain *domain = Domain::getInstance();
    auto bsize = domain->get_size();

    auto d_u = field_u->data;
    auto d_v = field_v->data;
    auto d_w = field_w->data;
    auto d_p = field_p->data;
    auto d_rhs = field_rhs->data;
    auto d_T = field_T->data;
    auto d_T_a = field_T_ambient->data;
    auto d_C = field_concentration->data;
    auto d_f_x = field_force_x->data;
    auto d_f_y = field_force_y->data;
    auto d_f_z = field_force_z->data;
    auto d_S_T = field_source_T->data;
    auto d_S_C = field_source_concentration->data;
    auto d_nu_t = field_nu_t->data;
    auto d_kappa_t = field_kappa_t->data;
    auto d_gamma_t = field_gamma_t->data;
#pragma acc exit data delete(d_u[:bsize], d_v[:bsize], d_w[:bsize], d_p[:bsize], d_rhs[:bsize], d_T[:bsize], d_T_a[:bsize], d_C[:bsize], d_f_x[:bsize], d_f_y[:bsize], d_f_z[:bsize], d_S_T[:bsize], d_S_C[:bsize], d_nu_t[:bsize], d_kappa_t[:bsize], d_gamma_t[:bsize])

    auto d_u0 = field_u0->data;
    auto d_v0 = field_v0->data;
    auto d_w0 = field_w0->data;
    auto d_u_tmp = field_u_tmp->data;
    auto d_v_tmp = field_v_tmp->data;
    auto d_w_tmp = field_w_tmp->data;
    auto d_p0 = field_p0->data;
    auto d_T0 = field_T0->data;
    auto d_T_tmp = field_T_tmp->data;
    auto d_C0 = field_concentration0->data;
    auto d_C_tmp = field_concentration_tmp->data;
#pragma acc exit data delete(d_u0[:bsize], d_u_tmp[:bsize], d_v0[:bsize], d_v_tmp[:bsize], d_w0[:bsize], d_w_tmp[:bsize], d_p0[:bsize], d_T0[:bsize], d_T_tmp[:bsize], d_C0[:bsize], d_C_tmp[:bsize])

    delete field_u;
    delete field_v;
    delete field_w;
    delete field_u0;
    delete field_v0;
    delete field_w0;
    delete field_u_tmp;
    delete field_v_tmp;
    delete field_w_tmp;

    delete field_nu_t;
    delete field_kappa_t;
    delete field_gamma_t;

    delete field_p;
    delete field_p0;
    delete field_rhs;

    delete field_T;
    delete field_T0;
    delete field_T_tmp;
    delete field_T_ambient;

    delete field_concentration;
    delete field_concentration0;
    delete field_concentration_tmp;

    delete field_force_x;
    delete field_force_y;
    delete field_force_z;
    delete field_source_T;
    delete field_source_concentration;

    delete sight;
}

// ========================================== Set up boundary =======================================
// ***************************************************************************************
/// \brief  initializes boundary cells
// ***************************************************************************************
void FieldController::set_up_boundary() {
    auto boundary = BoundaryController::getInstance();
    boundary->applyBoundary(field_u->data, field_u->get_type());
    boundary->applyBoundary(field_v->data, field_v->get_type());
    boundary->applyBoundary(field_w->data, field_w->get_type());
    boundary->applyBoundary(field_p->data, field_p->get_type());
    boundary->applyBoundary(field_T->data, field_T->get_type());
    boundary->applyBoundary(field_concentration->data, field_concentration->get_type());

    // TODO necessary?
    boundary->applyBoundary(field_T_ambient->data, field_T_ambient->get_type());
}

void FieldController::set_up_temporary_fields() {
    // copy constructor
    field_u0 = new Field(*field_u);
    field_u_tmp = new Field(*field_u);

    field_v0 = new Field(*field_v);
    field_v_tmp = new Field(*field_v);

    field_w0 = new Field(*field_w);
    field_w_tmp = new Field(*field_w);

    field_T0 = new Field(*field_T);
    field_T_tmp = new Field(*field_T);

    field_concentration0 = new Field(*field_concentration);
    field_concentration_tmp = new Field(*field_concentration);

    field_p0 = new Field(*field_p);

    auto d_u0 = field_u0->data;
    auto d_v0 = field_v0->data;
    auto d_w0 = field_w0->data;
    auto d_u_tmp = field_u_tmp->data;
    auto d_v_tmp = field_v_tmp->data;
    auto d_w_tmp = field_w_tmp->data;
    auto d_p0 = field_p0->data;
    auto d_T0 = field_T0->data;
    auto d_T_tmp = field_T_tmp->data;
    auto d_C0 = field_concentration0->data;
    auto d_C_tmp = field_concentration_tmp->data;

    Domain *domain = Domain::getInstance();
    auto bsize = domain->get_size();
    // copyin fields
#pragma acc enter data copyin(  d_u0[:bsize], d_u_tmp[:bsize], \
                                d_v0[:bsize], d_v_tmp[:bsize], \
                                d_w0[:bsize], d_w_tmp[:bsize], \
                                d_p0[:bsize], \
                                d_T0[:bsize], d_T_tmp[:bsize], \
                                d_C0[:bsize], d_C_tmp[:bsize])

}

//======================================= Update data ==================================
// ***************************************************************************************
/// \brief  Updates variables for the next iteration step or time dependent parameters such as temperature source function
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void FieldController::update_data(bool sync) {
    // local variables and parameters for GPU
    auto bsize = Domain::getInstance()->get_size();

    const auto d_u = field_u->data;                        //due to const correctness
    const auto d_v = field_v->data;
    const auto d_w = field_w->data;
    auto d_u0 = field_u0->data;
    auto d_v0 = field_v0->data;
    auto d_w0 = field_w0->data;
    auto d_u_tmp = field_u_tmp->data;
    auto d_v_tmp = field_v_tmp->data;
    auto d_w_tmp = field_w_tmp->data;
    const auto d_p = field_p->data;
    auto d_p0 = field_p0->data;
    const auto d_T = field_T->data;
    auto d_T0 = field_T0->data;
    auto d_T_tmp = field_T_tmp->data;
    const auto d_C = field_concentration->data;
    auto d_C0 = field_concentration0->data;
    auto d_C_tmp = field_concentration_tmp->data;

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t bsize_i = boundary->getSize_innerList();
    size_t *d_bList = boundary->get_boundaryList_level_joined();
    size_t bsize_b = boundary->getSize_boundaryList();
    size_t *d_oList = boundary->get_obstacleList();
    size_t bsize_o = boundary->getSize_obstacleList();

#pragma acc data present(    d_iList[:bsize_i], d_bList[:bsize_b], d_oList[:bsize_o], d_u[:bsize], d_v[:bsize], d_w[:bsize], \
                            d_u0[:bsize], d_v0[:bsize], d_w0[:bsize], d_u_tmp[:bsize], d_v_tmp[:bsize], d_w_tmp[:bsize], \
                            d_p[:bsize], d_p0[:bsize], d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], d_C[:bsize], d_C0[:bsize], d_C_tmp[:bsize])
    {
        // inner
#pragma acc parallel loop independent present(    d_iList[:bsize_i], d_u[:bsize], d_v[:bsize], d_w[:bsize], \
                                                d_u0[:bsize], d_v0[:bsize], d_w0[:bsize], d_u_tmp[:bsize], d_v_tmp[:bsize], d_w_tmp[:bsize], \
                                                d_p[:bsize], d_p0[:bsize], d_T[:bsize], d_T0[:bsize], d_T_tmp[:bsize], \
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
void FieldController::couple_vector(const Field *a, Field *a0, Field *a_tmp, const Field *b, Field *b0, Field *b_tmp, const Field *c, Field *c0, Field *c_tmp, bool sync) {

    // local variables and parameters for GPU
    auto d_a = a->data;
    auto d_a0 = a0->data;
    auto d_a_tmp = a_tmp->data;
    auto d_b = b->data;
    auto d_b0 = b0->data;
    auto d_b_tmp = b_tmp->data;
    auto d_c = c->data;
    auto d_c0 = c0->data;
    auto d_c_tmp = c_tmp->data;

    auto size = Domain::getInstance()->get_size(a0->get_level());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t bsize_i = boundary->getSize_innerList();
    size_t *d_bList = boundary->get_boundaryList_level_joined();
    size_t bsize_b = boundary->getSize_boundaryList();
    size_t *d_oList = boundary->get_obstacleList();
    size_t bsize_o = boundary->getSize_obstacleList();

#pragma acc data present(    d_a[:size], d_a0[:size], d_a_tmp[:size], d_b[:size], d_b0[:size], d_b_tmp[:size], \
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
    }//end data region
}

//======================================= Couple scalar ==================================
// ***************************************************************************************
/// \brief  couples vector (sets tmp and zero-th variables to current numerical solution)
/// \param  a   current field in x- direction (const)
/// \param  a0    zero-th field in x- direction
/// \param  a_tmp temporal field in x- direction
/// \param  sync  synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void FieldController::couple_scalar(const Field *a, Field *a0, Field *a_tmp, bool sync) {

    // local variables and parameters for GPU
    auto d_a = a->data;
    auto d_a0 = a0->data;
    auto d_a_tmp = a_tmp->data;

    auto size = Domain::getInstance()->get_size(a0->get_level());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t bsize_i = boundary->getSize_innerList();
    size_t *d_bList = boundary->get_boundaryList_level_joined();
    size_t bsize_b = boundary->getSize_boundaryList();
    size_t *d_oList = boundary->get_obstacleList();
    size_t bsize_o = boundary->getSize_obstacleList();

#pragma acc data present(d_a[:size], d_a0[:size], d_a_tmp[:size], d_iList[:bsize_i], d_bList[:bsize_b], d_oList[:bsize_o])
    {
        // inner
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            d_a0[i] = d_a[i];
            d_a_tmp[i] = d_a[i];
        }
        // boundary
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_b; ++j) {
            const size_t i = d_bList[j];
            d_a0[i] = d_a[i];
            d_a_tmp[i] = d_a[i];
        }
        // obstacles
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = 0; j < bsize_o; ++j) {
            const size_t i = d_oList[j];
            d_a0[i] = d_a[i];
            d_a_tmp[i] = d_a[i];
        }
        if (sync) {
#pragma acc wait
        }
    }//end data region
}
