/// \file       VCycleMG.cpp
/// \brief      Defines V-cycle of geometric multigrid method
/// \date       Sep 14, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>

#include "VCycleMG.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"
#include "../diffusion/ColoredGaussSeidelDiffuse.h"
#include "../diffusion/JacobiDiffuse.h"
#include "../solver/SolverSelection.h"
#include "../utility/Parameters.h"
#include "../utility/Utility.h"


// =============================== Constructor ===============================
// *****************************************************************************
/// \brief  Constructor
/// \param  out     pressure
/// \param  b       rhs
// *****************************************************************************
VCycleMG::VCycleMG(Field *out, Field *b) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    Parameters *params = Parameters::getInstance();
    Domain *domain = Domain::getInstance();

    m_levels = domain->get_levels();
    m_cycles = params->get_int("solver/pressure/n_cycle");
    m_relaxs = params->get_int("solver/pressure/diffusion/n_relax");

    m_dsign = -1.;
    m_w = params->get_real("solver/pressure/diffusion/w");

    // copies of out and b to prevent aliasing

// residuum
    auto *b_res1 = new Field(*b);
    m_residuum1.push_back(b_res1);
    // TODO(issue 124)
    real *data_b_res1 = b_res1->data;
    size_t bsize_b_res1 = b_res1->get_size();
#pragma acc enter data copyin(data_b_res1[:bsize_b_res1])

// error
    auto *out_err1 = new Field(*out);
    m_error1.push_back(out_err1);
    // TODO(issue 124)
    real *data_err1 = out_err1->data;
    size_t bsize_err1 = out_err1->get_size();
#pragma acc enter data copyin(data_err1[:bsize_err1])

// temporal solution
    auto *out_tmp = new Field(*out);
    m_mg_temporal_solution.push_back(out_tmp);
    // TODO(issue 124)
    real *data_mg_temporal_solution = out_tmp->data;
    size_t bsize_mg_temporal_solution = out_tmp->get_size();
#pragma acc enter data copyin(data_mg_temporal_solution[:bsize_mg_temporal_solution])


    // building fields for level + sending to GPU
    // level going up
    for (size_t i = 0; i < m_levels; ++i) {
        // build m_residuum0
        auto *r0 = new Field(FieldType::P, 0.0, i);
        m_residuum0.push_back(r0);
        // TODO(issue 124)
        real *data_residuum0 = r0->data;
        size_t bsize_residuum0 = r0->get_size();
#pragma acc enter data copyin(data_residuum0[:bsize_residuum0])

        // build m_residuum1
        auto *r1 = new Field(FieldType::P, 0.0, i + 1);
        m_residuum1.push_back(r1);
        // TODO(issue 124)
        real *data_residuum1 = r1->data;
        size_t bsize_residuum1 = r1->get_size();
#pragma acc enter data copyin(data_residuum1[:bsize_residuum1])

        //  build m_error1
        auto *e1 = new Field(FieldType::P, 0.0, i + 1);
        m_error1.push_back(e1);
        // TODO(issue 124)
        real *data_e1 = e1->data;
        size_t bsize_e1 = e1->get_size();
#pragma acc enter data copyin(data_e1[:bsize_e1])

        // build m_mg_temporal_solution
        auto *mg = new Field(FieldType::P, 0.0, i + 1);  // new field to prevent aliasing
        m_mg_temporal_solution.push_back(mg);
        // TODO(issue 124)
        real *data_mg = mg->data;
        size_t bsize_mg = mg->get_size();
#pragma acc enter data copyin(data_mg[:bsize_mg])
    } // end level loop

//  build m_err0
    m_err0.resize(m_levels + 1);

    auto *e00 = new Field(FieldType::P, 0.0, m_error1[0]->get_level());
    m_err0[0] = e00;
    // TODO(issue 124)
    real *data_err00 = e00->data;
    size_t bsize_err00 = e00->get_size();
#pragma acc enter data copyin(data_err00[:bsize_err00])

    // building fields for level
    // level going down
    for (size_t i = m_levels; i > 0; --i) {
        // build m_err0
        auto *e0 = new Field(FieldType::P, 0.0, m_error1[i - 1]->get_level());
        m_err0[i] = e0;
        // TODO(issue 124)
        real *data_err0 = m_err0[i]->data;
        size_t bsize_err0 = m_err0[i]->get_size();
#pragma acc enter data copyin(data_err0[:bsize_err0])
    }
}

VCycleMG::~VCycleMG() {
    auto domain = Domain::getInstance();

    while (!m_residuum0.empty()) {
        Field *field = m_residuum0.back();
        // TODO(issue 124)
        real *data = field->data;
        size_t bsize = field->get_size();
#pragma acc exit data delete(data[:bsize])
        delete field;
        m_residuum0.pop_back();
    }

    while (!m_residuum1.empty()) {
        Field *field = m_residuum1.back();
        real *data = field->data;
        // TODO(issue 124)
        size_t bsize = field->get_size();
#pragma acc exit data delete(data[:bsize])
        delete field;
        m_residuum1.pop_back();
    }

    while (!m_err0.empty()) {
        Field *field = m_err0.back();
        real *data = field->data;
        // TODO(issue 124)
        size_t bsize = field->get_size();
#pragma acc exit data delete(data[:bsize])
        delete field;
        m_err0.pop_back();
    }

    while (!m_error1.empty()) {
        Field *field = m_error1.back();
        real *data = field->data;
        // TODO(issue 124)
        size_t bsize = field->get_size();
#pragma acc exit data delete(data[:bsize])
        delete field;
        m_error1.pop_back();
    }

    while (!m_mg_temporal_solution.empty()) {
        Field *field = m_mg_temporal_solution.back();
        real *data = field->data;
        // TODO(issue 124)
        size_t bsize = field->get_size();
#pragma acc exit data delete(data[:bsize])
        delete field;
        m_mg_temporal_solution.pop_back();
    }
}

// =============================== Update ===============================
// ************************************************************************
/// \brief  Update input
/// \param  out     pressure
/// \param  b       rhs
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ************************************************************************
void VCycleMG::UpdateInput(Field *out, Field *b, bool sync) {
    m_error1[0]->copy_data(*out);
    m_mg_temporal_solution[0]->copy_data(*out);
    m_residuum1[0]->copy_data(*b);
}

//==================================== Pressure =================================
// *****************************************************************************
/// \brief  solves Poisson equation \f$ \nabla^2 p = rhs\f$ via geometric multigrid (VCycle)
/// \param  out         pressure
/// \param  b           rhs
/// \param  t           current time
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::pressure(Field *out, Field *b, real t, bool sync) {
    // Update first
    UpdateInput(out, b);

    // solve more accurately, in first time step
    Parameters *params = Parameters::getInstance();
    Domain *domain = Domain::getInstance();
    const real dt = params->get_real("physical_parameters/dt");
    const auto Nt = static_cast<size_t>(std::round(t / dt));

    const int set_relaxs = params->get_int("solver/pressure/diffusion/n_relax");
    const int set_cycles = params->get_int("solver/pressure/n_cycle");

    int act_cycles = 0;

    if (Nt == 1) {
        const int max_cycles = params->get_int("solver/pressure/max_cycle");
        const int max_relaxs = params->get_int("solver/pressure/diffusion/max_solve");

        real r = 10000.;
        real sum = 0;
        const real tol_res = params->get_real("solver/pressure/tol_res");

        const size_t Nx = domain->get_Nx();
        const size_t Ny = domain->get_Ny();
        size_t bsize = domain->get_size();

        const real dx = domain->get_dx();
        const real dy = domain->get_dy();
        const real dz = domain->get_dz();

        const real rdx2 = 1. / (dx * dx);
        const real rdy2 = 1. / (dy * dy);
        const real rdz2 = 1. / (dz * dz);

        auto data_out = out->data;
        auto data_b = b->data;

        auto boundary = BoundaryController::getInstance();
        auto bsize_i = boundary->get_size_inner_list();
        size_t *data_inner_list = boundary->get_inner_list_level_joined();

        while (r > tol_res && act_cycles < max_cycles && m_relaxs < max_relaxs) {
            for (size_t i = 0; i < m_cycles; i++) {
                VCycleMultigrid(out, sync);
                act_cycles++;
            }
            sum = 0.;

            // calculate residuum in inner cells
#pragma acc parallel loop independent present(data_out[:bsize], data_b[:bsize], data_inner_list[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = data_inner_list[j];
                r = data_b[i] - (rdx2 * (data_out[i - 1]       - 2 * data_out[i] + data_out[i + 1])
                              +  rdy2 * (data_out[i - Nx]      - 2 * data_out[i] + data_out[i + Nx])
                              +  rdz2 * (data_out[i - Nx * Ny] - 2 * data_out[i] + data_out[i + Nx * Ny]));
                sum += r * r;
            }
#pragma acc wait
            r = sqrt(sum);
            m_relaxs += set_relaxs;
        }
    } else {  // Nt > 1
        m_cycles = set_cycles;
        m_relaxs = set_relaxs;

        for (size_t i = 0; i < m_cycles; i++) {
            VCycleMultigrid(out, sync);
        }
    }

    if (sync) {
#pragma acc wait
    }
}

//==================================== VCycle =================================
// *****************************************************************************
/// \brief  Conducts the V-cycle Multigrid method
/// \param  field_out         pressure
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::VCycleMultigrid(Field *field_out, bool sync) {
    size_t max_level = m_levels;

    auto domain = Domain::getInstance();
    auto boundary = BoundaryController::getInstance();
    size_t *data_inner_list = boundary->get_inner_list_level_joined();
    size_t *data_boundary_list = boundary->get_boundary_list_level_joined();

//===================== No refinement, when levels = 0 =========//
    if (max_level == 0) {
        Field *field_mg_tmpi = m_mg_temporal_solution[0];
        Field *field_res1i = m_residuum1[0];

        real *data_mg_tmpi = m_mg_temporal_solution[0]->data;
        real *data_res1i = m_residuum1[0]->data;
        real *data_out = field_out->data;

        size_t size_mg_tmpi = m_mg_temporal_solution[0]->get_size();
        size_t size_res1i = m_residuum1[0]->get_size();
        size_t size_out = field_out->get_size();

#pragma acc data present(data_out[:size_out], data_mg_tmpi[:size_mg_tmpi], data_res1i[:size_res1i])
        {
            Solve(field_out, field_mg_tmpi, field_res1i, max_level, sync);
        }
        return;
    }

//===================== levels going down ====================//
    for (size_t i = 0; i < max_level; ++i) {
        Field *field_res0i = m_residuum0[i];
        Field *field_err1i = m_error1[i];
        Field *field_err1ip = m_error1[i + 1];
        Field *field_mg_tmpi = m_mg_temporal_solution[i];
        Field *field_res1i = m_residuum1[i];
        Field *field_res1ip = m_residuum1[i + 1];

        real *data_res0i = m_residuum0[i]->data;
        real *data_err1i = m_error1[i]->data;
        real *data_err1ip = m_error1[i + 1]->data;
        real *data_mg_tmpi = m_mg_temporal_solution[i]->data;
        real *data_res1i = m_residuum1[i]->data;
        real *data_res1ip = m_residuum1[i + 1]->data;
        real *data_out = field_out->data;

        size_t size_res0i = m_residuum0[i]->get_size();
        size_t size_err1i = m_error1[i]->get_size();
        size_t size_err1ip = m_error1[i + 1]->get_size();
        size_t size_mg_tmpi = m_mg_temporal_solution[i]->get_size();
        size_t size_res1i = m_residuum1[i]->get_size();
        size_t size_res1ip = m_residuum1[i + 1]->get_size();
        size_t size_out = field_out->get_size();

        FieldType type_r0 = field_res0i->get_type();

        FieldType type_r1 = m_residuum1[i + 1]->get_type();

#pragma acc data present(data_res0i[:size_res0i], data_err1i[:size_err1i], data_err1ip[:size_err1ip], \
                         data_mg_tmpi[:size_mg_tmpi], data_res1i[:size_res1i], data_res1ip[:size_res1ip], \
                         data_out[:size_out])
        {
            if (i == 0) {  // use p=field_out on finest grid
                // smooth
                Smooth(field_out, field_mg_tmpi, field_res1i, i, sync);

                // calculate residuum
                Residuum(field_res0i, field_out, field_res1i, i, sync);
                boundary->apply_boundary(data_res0i, i, type_r0, sync); // for m_residuum0 only Dirichlet BC
            } else {
                // smooth
                Smooth(field_err1i, field_mg_tmpi, field_res1i, i, sync);

                // calculate residuum
                Residuum(field_res0i, field_err1i, field_res1i, i, sync);
                boundary->apply_boundary(data_res0i, i, type_r0, sync); // for m_residuum0 only Dirichlet BC
            }

            // restrict
            Restrict(field_res1ip, field_res0i, i, sync);
            boundary->apply_boundary(data_res1ip, i + 1, type_r1, sync); // for res only Dirichlet BC

            // set err to zero at next level
            // TODO(lukas) set everything to 0 ? currently the obstacle cells are excluded
            field_err1ip->set_value(0);

            /*
            // strides (since GPU needs joined list)
            // inner start/ end index of level i + 1
            size_t start_i = boundary->get_inner_list_level_joined_start(i + 1);
            size_t end_i = boundary->get_inner_list_level_joined_end(i + 1) + 1;
            // boundary start/ end index of level i + 1
            size_t start_b = boundary->get_boundary_list_level_joined_start(i + 1);
            size_t end_b = boundary->get_boundary_list_level_joined_end(i + 1) + 1;
            // inner
#pragma acc kernels present(data_err1ip[:size_err1ip], data_inner_list[start_i:(end_i-start_i)]) async
#pragma acc loop independent
            for (size_t j = start_i; j < end_i; ++j) {
                const size_t idx = data_inner_list[j];
                data_err1ip[idx] = 0.0;
            }

            //boundary
#pragma acc kernels present(data_err1ip[:size_err1ip], data_boundary_list[start_b:(end_b-start_b)]) async
#pragma acc loop independent
            for (size_t j = start_b; j < end_b; ++j) {
                const size_t idx = data_boundary_list[j];
                data_err1ip[idx] = 0.0;
            }
             */
        }
    }

//===================== levels going up ====================//
    for (size_t i = max_level; i > 0; --i) {
        Field *field_err0i = m_err0[i];
        Field *field_err1i = m_error1[i];
        Field *field_err1im = m_error1[i - 1];
        Field *field_mg_tmpim = m_mg_temporal_solution[i - 1];
        Field *field_res1im = m_residuum1[i - 1];

        real *data_err0i = m_err0[i]->data;
        real *data_err1i = m_error1[i]->data;
        real *data_err1im = m_error1[i - 1]->data;
        real *data_mg_tmpim = m_mg_temporal_solution[i - 1]->data;
        real *data_res1im = m_residuum1[i - 1]->data;
        real *data_out = field_out->data;

        size_t size_err0i = m_err0[i]->get_size();
        size_t size_err1i = m_error1[i]->get_size();
        size_t size_err1im = m_error1[i - 1]->get_size();
        size_t size_mg_tmpim = m_mg_temporal_solution[i - 1]->get_size();
        size_t size_res1im = m_residuum1[i - 1]->get_size();
        size_t size_out = field_out->get_size();

        FieldType type_e0 = field_err0i->get_type();

        // inner start/end index of level i - 1
        size_t start_i = boundary->get_inner_list_level_joined_start(i - 1);
        size_t end_i = boundary->get_inner_list_level_joined_end(i - 1) + 1;
        // boundary start/end index of level i - 1
        size_t start_b = boundary->get_boundary_list_level_joined_start(i - 1);
        size_t end_b = boundary->get_boundary_list_level_joined_end(i - 1) + 1;

#pragma acc data present(data_err0i[:size_err0i], data_err1i[:size_err1i], data_err1im[:size_err1im], data_mg_tmpim[:size_mg_tmpim], data_res1im[:size_res1im], data_out[:size_out])
        {
            // prolongate
            Prolongate(field_err0i, field_err1i, i, sync);
            boundary->apply_boundary(data_err0i, i - 1, type_e0, sync); // for m_err0 only Dirichlet BC

            // correct
            if (i == 1) {  // use p=field_out on finest grid
                // inner
#pragma acc kernels present(data_err0i[:size_err0i], data_out[:size_out], data_inner_list[start_i:(end_i-start_i)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_i; j < end_i; ++j) {
                        const size_t idx = data_inner_list[j];
                        data_out[idx] += data_err0i[idx];
                    }
                }
                // boundary
#pragma acc kernels present(data_err0i[:size_err0i], data_out[:size_out], data_boundary_list[start_b:(end_b-start_b)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_b; j < end_b; ++j) {
                        const size_t idx = data_boundary_list[j];
                        data_out[idx] += data_err0i[idx];
                    }
                }
                // smooth
                Smooth(field_out, field_mg_tmpim, field_res1im, i - 1, sync);
            } else {
                // correct
                // inner
#pragma acc kernels present(data_err0i[:size_err0i], data_err1im[:size_err1im], data_inner_list[start_i:(end_i-start_i)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_i; j < end_i; ++j) {
                        const size_t idx = data_inner_list[j];
                        data_err1im[idx] += data_err0i[idx];
                    }
                }
                // boundary
#pragma acc kernels present(data_err0i[:size_err0i], data_err1im[:size_err1im], data_boundary_list[start_b:(end_b-start_b)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_b; j < end_b; ++j) {
                        const size_t idx = data_boundary_list[j];
                        data_err1im[idx] += data_err0i[idx];
                    }
                }
                // smooth
                if (i - 1 == m_levels - 1) {
                    Solve(field_err1im, field_mg_tmpim, field_res1im, i - 1, sync);
                } else {
                    Smooth(field_err1im, field_mg_tmpim, field_res1im, i - 1, sync); // for err only Dirichlet BC
                }
            }
        }
    }

    // set boundaries
    real *data_out = field_out->data;
    FieldType type = field_out->get_type();
    boundary->apply_boundary(data_out, type, sync);

    if (sync) {
#pragma acc wait
    }
}

//==================================== Smooth =================================
// *****************************************************************************
/// \brief  Relaxes Ax = b using a diffusion method
/// \param  out         output field (in size of input )
/// \param  tmp         temporary field for JacobiStep
/// \param  b           right hand side
/// \param  level       Multigrid level
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::Smooth(Field *out, Field *tmp, Field *b, size_t level, bool sync) {
    Domain *domain = Domain::getInstance();

    // local variables and parameters for GPU
    const real dx = domain->get_dx(out->get_level());
    const real dy = domain->get_dy(out->get_level());
    const real dz = domain->get_dz(out->get_level());

    tmp->copy_data(*out);
    real *data_out = out->data;
    real *data_tmp = tmp->data;
    real *data_b = b->data;

    Parameters *params = Parameters::getInstance();

    size_t bsize = domain->get_size(out->get_level());
    FieldType type = out->get_type();

    BoundaryController *boundary = BoundaryController::getInstance();

    size_t *data_inner_list = boundary->get_inner_list_level_joined();
    size_t *data_boundary_list = boundary->get_boundary_list_level_joined();

    // start/end index of level
    // inner
    size_t start_i = boundary->get_inner_list_level_joined_start(level);
    size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;
    // boundary
    size_t start_b = boundary->get_boundary_list_level_joined_start(level);
    size_t end_b = boundary->get_boundary_list_level_joined_end(level) + 1;

    // apply boundary: at level 0 apply set BC; else use Dirichlet 0
    boundary->apply_boundary(data_out, level, type, sync);

    // diffuse
    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    // preparation for diffusion step (alpha, beta)
    real alphaX = rdx2;
    real alphaY = rdy2;
    real alphaZ = rdz2;

    const real rbeta = 2. * (alphaX + alphaY + alphaZ);
    const real beta = 1. / rbeta;

    // Diffusion step
    std::string diffusionType = params->get("solver/pressure/diffusion/type");
    if (diffusionType == DiffusionMethods::Jacobi) {
#pragma acc data present(data_out[:bsize], data_tmp[:bsize], data_b[:bsize])
        {
            for (int i = 0; i < m_relaxs; i++) { // fixed iteration number as in xml
                JacobiDiffuse::JacobiStep(level, out, tmp, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
                boundary->apply_boundary(data_out, level, type, sync);

                std::swap(tmp->data, out->data);
                std::swap(data_tmp, data_out);
            }

            if (m_relaxs % 2 != 0) { // swap necessary when odd number of iterations
#pragma acc kernels present(data_out[:bsize], data_tmp[:bsize], data_inner_list[start_i:(end_i-start_i)]) async
#pragma acc loop independent
                // inner
                for (size_t j = start_i; j < end_i; ++j) {
                    const size_t i = data_inner_list[j];
                    data_out[i] = data_tmp[i];
                }
                // boundary
#pragma acc kernels present(data_out[:bsize], data_tmp[:bsize], data_boundary_list[start_b:(end_b-start_b)]) async
#pragma acc loop independent
                for (size_t j = start_b; j < end_b; ++j) {
                    const size_t i = data_boundary_list[j];
                    data_out[i] = data_tmp[i];
                }
            }
        }
    } else if (diffusionType == DiffusionMethods::ColoredGaussSeidel) {

#pragma acc data present(data_out[:bsize], data_tmp[:bsize], data_b[:bsize])
        {
            for (int i=0; i < m_relaxs; i++) {
                ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(out, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
                boundary->apply_boundary(data_out, level, type, sync); // for res/err only Dirichlet BC
            }
        }
    } else {
#ifndef BENCHMARKING
        m_logger->critical("Diffusion method not yet implemented! Simulation stopped!");
#endif
        std::exit(1);
        // TODO Error handling
    }

    if (sync) {
#pragma acc wait
    }
}

//================================== Residuum ===============================
// ************************************************************************
/// \brief  Calculates residuum r = b - Ax
/// \param  out         output field (in size of input field)
/// \param  in          input field
/// \param  b           right hand side field
/// \param  level       Multigrid level
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// ************************************************************************
void VCycleMG::Residuum(Field *out, Field *in, Field *b, size_t level, bool sync) {
    Domain *domain = Domain::getInstance();

    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(in->get_level());
    const size_t Ny = domain->get_Ny(in->get_level());

    const real dx = domain->get_dx(in->get_level());
    const real dy = domain->get_dy(in->get_level());
    const real dz = domain->get_dz(in->get_level());

    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    real *data_out = out->data;
    real *data_in = in->data;
    real *data_b = b->data;

    size_t bsize = domain->get_size(out->get_level());

    auto boundary = BoundaryController::getInstance();
    size_t *data_inner_list = boundary->get_inner_list_level_joined();

    // starts/ends
    // inner
    size_t start_i = boundary->get_inner_list_level_joined_start(level);
    size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;
#pragma acc data present(data_b[:bsize], data_in[:bsize], data_out[:bsize], data_inner_list[start_i:(end_i-start_i)])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = start_i; j < end_i; ++j) {
            const size_t i = data_inner_list[j];
            data_out[i] = data_b[i] - (rdx2 * (data_in[i - 1]       - 2 * data_in[i] + data_in[i + 1])\
                                    +  rdy2 * (data_in[i - Nx]      - 2 * data_in[i] + data_in[i + Nx])\
                                    +  rdz2 * (data_in[i - Nx * Ny] - 2 * data_in[i] + data_in[i + Nx * Ny]));
        }

        if (sync) {
#pragma acc wait
        }
    }//end data region
}

//================================== Restrict ===============================
// *****************************************************************************
/// \brief  Restricts field from fine grid to coarse grid via averaging
/// \param  out         output field (half the size of input vector)
/// \param  in          input field (on fine grid)
/// \param  level       Multigrid level
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::Restrict(Field *out, Field *in, size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    // coarse grid
    const size_t Nx = domain->get_Nx(out->get_level());
    const size_t Ny = domain->get_Ny(out->get_level());

    // fine grid
    const size_t nx = domain->get_Nx(in->get_level());
    const size_t ny = domain->get_Ny(in->get_level());

    auto data_out = out->data;
    auto data_in = in->data;

    size_t bsize_out = domain->get_size(out->get_level());
    size_t bsize_in = domain->get_size(in->get_level());

    auto boundary = BoundaryController::getInstance();
    size_t *data_inner_list = boundary->get_inner_list_level_joined();

    // start/end
    // inner
    size_t start_i = boundary->get_inner_list_level_joined_start(level + 1);
    size_t end_i = boundary->get_inner_list_level_joined_end(level + 1) + 1;

#ifndef BENCHMARKING
    if (end_i == start_i) {
        m_logger->warn("Be cautious: Obstacle might fill up inner cells completely in level {} with nx= {}!",
                       level, domain->get_nx(out->get_level()));
    }
#endif

    // average from eight neighboring cells
    // obstacles not used in fine grid, since coarse grid only obstacle if one of 8 fine grids was an obstacle,
    // thus if coarse cell inner cell, then surrounding fine cells also inner cells!
    size_t i, j, k;

#pragma acc data present(data_in[:bsize_in], data_out[:bsize_out], data_inner_list[start_i:(end_i-start_i)])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = start_i; l < end_i; ++l) {
            const size_t idx = data_inner_list[l];
            k = getCoordinateK(idx, Nx, Ny);
            j = getCoordinateJ(idx, Nx, Ny, k);
            i = getCoordinateI(idx, Nx, Ny, j, k);

            data_out[idx] = 0.125 * (data_in[IX(2 * i - 1, 2 * j - 1, 2 * k - 1, nx, ny)]\
                                   + data_in[IX(2 * i,     2 * j - 1, 2 * k - 1, nx, ny)]\
                                   + data_in[IX(2 * i - 1, 2 * j,     2 * k - 1, nx, ny)]\
                                   + data_in[IX(2 * i,     2 * j,     2 * k - 1, nx, ny)]\
                                   + data_in[IX(2 * i - 1, 2 * j - 1, 2 * k, nx, ny)]\
                                   + data_in[IX(2 * i,     2 * j - 1, 2 * k, nx, ny)]\
                                   + data_in[IX(2 * i - 1, 2 * j,     2 * k, nx, ny)]\
                                   + data_in[IX(2 * i,     2 * j,     2 * k, nx, ny)]);

        }

        if (sync) {
#pragma acc wait
        }
    }
}

//================================== Prolongate ===============================
// *****************************************************************************
/// \brief  Prolongates field from coarse grid to fine grid (trilinear interpolation)
/// \param  out         output field (real the size of input vector)
/// \param  in          input field (on coarse grid)
/// \param  level       Multigrid level
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::Prolongate(Field *out, Field *in, size_t level, bool sync) {
    Domain *domain = Domain::getInstance();

    // local variables and parameters for GPU
    // fine grid
    const size_t nx = domain->get_Nx(out->get_level());
    const size_t ny = domain->get_Ny(out->get_level());

    // coarse grid
    const size_t Nx = domain->get_Nx(in->get_level());
    const size_t Ny = domain->get_Ny(in->get_level());

    real *data_out = out->data;
    real *data_in = in->data;

    size_t bsize_out = domain->get_size(out->get_level());
    size_t bsize_in = domain->get_size(in->get_level());

    BoundaryController *boundary = BoundaryController::getInstance();

    size_t *data_inner_list = boundary->get_inner_list_level_joined();

    // start/end (going backwards)
    // inner
    size_t start_i = boundary->get_inner_list_level_joined_start(level);
    size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;

    // prolongate
    size_t i, j, k;

#pragma acc data present(data_in[:bsize_in], data_out[:bsize_out], data_inner_list[start_i:(end_i-start_i)])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = start_i; l < end_i; ++l) {
            const size_t idx = data_inner_list[l];
            k = idx / (Nx * Ny);
            j = (idx - k * Nx * Ny) / Nx;
            i = idx - k * Nx * Ny - j * Nx;

            data_out[IX(2 * i, 2 * j, 2 * k, nx, ny)] =
                    0.015625 * (27 * data_in[idx]
                               + 9 * data_in[idx + 1]
                               + 9 * data_in[idx + Nx]
                               + 9 * data_in[idx + Nx * Ny]
                               + 3 * data_in[idx + 1 + Nx]
                               + 3 * data_in[idx + 1 + Nx * Ny]
                               + 3 * data_in[idx + Nx + Nx * Ny]
                               +     data_in[idx + 1 + Nx + Nx * Ny]);

            data_out[IX(2 * i, 2 * j, 2 * k - 1, nx, ny)] =
                    0.015625 * (27 * data_in[idx]
                               + 9 * data_in[idx + 1]
                               + 9 * data_in[idx + Nx]
                               + 9 * data_in[idx - Nx * Ny]
                               + 3 * data_in[idx + 1 + Nx]
                               + 3 * data_in[idx + 1 - Nx * Ny]
                               + 3 * data_in[idx + Nx - Nx * Ny]
                               +     data_in[idx + 1 + Nx - Nx * Ny]);

            data_out[IX(2 * i, 2 * j - 1, 2 * k, nx, ny)] =
                    0.015625 * (27 * data_in[idx]
                               + 9 * data_in[idx + 1]
                               + 9 * data_in[idx - Nx]
                               + 9 * data_in[idx + Nx * Ny]
                               + 3 * data_in[idx + 1 - Nx]
                               + 3 * data_in[idx + 1 + Nx * Ny]
                               + 3 * data_in[idx - Nx + Nx * Ny]
                               +     data_in[idx + 1 - Nx + Nx * Ny]);

            data_out[IX(2 * i, 2 * j - 1, 2 * k - 1, nx, ny)] =
                    0.015625 * (27 * data_in[idx]
                               + 9 * data_in[idx + 1]
                               + 9 * data_in[idx - Nx]
                               + 9 * data_in[idx - Nx * Ny]
                               + 3 * data_in[idx + 1 - Nx]
                               + 3 * data_in[idx + 1 - Nx * Ny]
                               + 3 * data_in[idx - Nx - Nx * Ny]
                               +     data_in[idx + 1 - Nx - Nx * Ny]);

            data_out[IX(2 * i - 1, 2 * j, 2 * k, nx, ny)] =
                    0.015625 * (27 * data_in[idx]
                               + 9 * data_in[idx - 1]
                               + 9 * data_in[idx + Nx]
                               + 9 * data_in[idx + Nx * Ny]
                               + 3 * data_in[idx - 1 + Nx]
                               + 3 * data_in[idx - 1 + Nx * Ny]
                               + 3 * data_in[idx + Nx + Nx * Ny]
                               +     data_in[idx - 1 + Nx + Nx * Ny]);

            data_out[IX(2 * i - 1, 2 * j, 2 * k - 1, nx, ny)] =
                    0.015625 * (27 * data_in[idx]
                               + 9 * data_in[idx - 1]
                               + 9 * data_in[idx + Nx]
                               + 9 * data_in[idx - Nx * Ny]
                               + 3 * data_in[idx - 1 + Nx]
                               + 3 * data_in[idx - 1 - Nx * Ny]
                               + 3 * data_in[idx + Nx - Nx * Ny]
                               +     data_in[idx - 1 + Nx - Nx * Ny]);

            data_out[IX(2 * i - 1, 2 * j - 1, 2 * k, nx, ny)] =
                    0.015625 * (27 * data_in[idx]
                               + 9 * data_in[idx - 1]
                               + 9 * data_in[idx - Nx]
                               + 9 * data_in[idx + Nx * Ny]
                               + 3 * data_in[idx - 1 - Nx]
                               + 3 * data_in[idx - 1 + Nx * Ny]
                               + 3 * data_in[idx - Nx + Nx * Ny]
                               +     data_in[idx - 1 - Nx + Nx * Ny]);

            data_out[IX(2 * i - 1, 2 * j - 1, 2 * k - 1, nx, ny)] =
                    0.015625 * (27 * data_in[idx]
                               + 9 * data_in[idx - 1]
                               + 9 * data_in[idx - Nx]
                               + 9 * data_in[idx - Nx * Ny]
                               + 3 * data_in[idx - 1 - Nx]
                               + 3 * data_in[idx - 1 - Nx * Ny]
                               + 3 * data_in[idx - Nx - Nx * Ny]
                               +     data_in[idx - 1 - Nx - Nx * Ny]);
        }

        if (sync) {
#pragma acc wait
        }
    }// end data region
}

//==================================== Smooth =================================
// *****************************************************************************
/// \brief  Solves Ax = b on the lowest level using a diffusion method
/// \param  out         output field (in size of input )
/// \param  tmp         temporary field for JacobiStep
/// \param  b           right hand side
/// \param  level       Multigrid level
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::Solve(Field *out, Field *tmp, Field *b, size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(out->get_level());
    const size_t Ny = domain->get_Ny(out->get_level());

    if (level < m_levels - 1) {
#ifndef BENCHMARKING
        m_logger->warn("Wrong level = {}", level);
#endif
        return;
        // TODO Error handling
    }

    if (Nx <= 4 && Ny <= 4) {
#ifndef BENCHMARKING
        m_logger->warn(" Grid is too coarse with Nx={} and Ny={}. Just smooth here", Nx, Ny);
#endif
        Smooth(out, tmp, b, level, sync);
        return;
        // TODO Error handling
    }

    const real dx = domain->get_dx(out->get_level());
    const real dy = domain->get_dy(out->get_level());
    const real dz = domain->get_dz(out->get_level());

    auto data_out = out->data;
    auto data_tmp = tmp->data;
    auto data_b = b->data;

    Parameters *params = Parameters::getInstance();

    size_t bsize = domain->get_size(out->get_level());
    FieldType type = out->get_type();

    BoundaryController *boundary = BoundaryController::getInstance();

    size_t *data_inner_list = boundary->get_inner_list_level_joined();
    size_t *data_boundary_list = boundary->get_boundary_list_level_joined();
    size_t bsize_b = boundary->get_size_boundary_list_level_joined();

    // start/end
    // inner
    size_t start_i = boundary->get_inner_list_level_joined_start(level);
    size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;

    boundary->apply_boundary(data_out, level, type, sync);

    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    // preparation for diffusion step (alpha, beta)
    real alphaX = rdx2;
    real alphaY = rdy2;
    real alphaZ = rdz2;

    const real rbeta = 2. * (alphaX + alphaY + alphaZ);
    const real beta = 1. / rbeta;

    tmp->copy_data(*out);

// Diffusion step
    std::string diffusionType = params->get("solver/pressure/diffusion/type");
    if (diffusionType == DiffusionMethods::Jacobi) {
#pragma acc data present(data_out[:bsize], data_tmp[:bsize], data_b[:bsize])
        {
            size_t it = 0;
            const size_t max_it = static_cast<size_t> (params->get_int("solver/pressure/diffusion/max_solve"));
            const real tol_res = params->get_real("solver/pressure/diffusion/tol_res");
            real sum;
            real res = 10000.;

            while (res > tol_res && it < max_it) {
                JacobiDiffuse::JacobiStep(level, out, tmp, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
                boundary->apply_boundary(data_out, level, type, sync);

                sum = 0.;

#pragma acc parallel loop independent present(data_out[:bsize], data_tmp[:bsize], data_inner_list[:bsize_i]) async
                for (size_t j = start_i; j < end_i; ++j) {
                    const size_t i = data_inner_list[j];
                    res = data_b[i] - (rdx2 * (data_out[i - 1] - 2 * data_out[i] + data_out[i + 1])\
                                     + rdy2 * (data_out[i - Nx] - 2 * data_out[i] + data_out[i + Nx])\
                                     + rdz2 * (data_out[i - Nx * Ny] - 2 * data_out[i] + data_out[i + Nx * Ny])); //res = rbeta*(data_out[i] - data_tmp[i]);
                    sum += res * res;
                }
                // info: in nvvp profile 8byte size copy from to device to/from pageable due to sum!

#pragma acc wait
                res = sqrt(sum);

                it++;

                std::swap(tmp->data, out->data);
                std::swap(data_tmp, data_out);
            }  // end while

            if (it % 2 != 0) {  // swap necessary when odd number of iterations
                out->copy_data(*tmp);
            }
        } //end data region
    } else if (diffusionType == DiffusionMethods::ColoredGaussSeidel) {
#pragma acc data present(data_out[:bsize], data_tmp[:bsize], data_b[:bsize])
        {
            size_t it = 0;
            const size_t max_it = static_cast<size_t>(params->get_int("solver/diffusion/max_iter"));
            const real tol_res = params->get_real("solver/diffusion/tol_res");
            real sum;
            real res = 10000.;

            while (res > tol_res && it < max_it) {
                ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(out, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
                boundary->apply_boundary(data_out, level, type, sync); // for res/err only Dirichlet BC

                sum = 0.;

#pragma acc parallel loop independent present(data_out[:bsize], data_tmp[:bsize], data_inner_list[:bsize_i]) async
                for (size_t j = start_i; j < end_i; ++j) {
                    const size_t i = data_inner_list[j];
                    res = data_b[i] - (rdx2 * (data_out[i - 1] - 2 * data_out[i] + data_out[i + 1])\
                                     + rdy2 * (data_out[i - Nx] - 2 * data_out[i] + data_out[i + Nx])\
                                     + rdz2 * (data_out[i - Nx * Ny] - 2 * data_out[i] + data_out[i + Nx * Ny]));  //res = rbeta*(data_out[i] - data_tmp[i]);
                    sum += res * res;
                }
                // info: in nvvp profile 8byte size copy from to device to/from pageable due to sum!

#pragma acc wait
                res = sqrt(sum);
                it++;
            }  // end while
        }  // end data region
    } else {
#ifndef BENCHMARKING
        m_logger->critical("Diffusion method not yet implemented! Simulation stopped!");
#endif
        std::exit(1);
        // TODO Error handling
    }

    if (sync) {
#pragma acc wait
    }
}
