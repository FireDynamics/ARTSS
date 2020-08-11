/// \file       VCycleMG.cpp
/// \brief      Defines V-cycle of geometric multigrid method
/// \date       Sep 14, 2016
/// \author     Severt
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#include <cmath>

#include "VCycleMG.h"
#include "../diffusion/JacobiDiffuse.h"
#include "../diffusion/ColoredGaussSeidelDiffuse.h"
#include "../utility/Parameters.h"
#include "../boundary/BoundaryController.h"
#include "../Domain.h"
#include "../solver/SolverSelection.h"
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

    auto params = Parameters::getInstance();
    auto domain = Domain::getInstance();

    levels = domain->get_levels();
    cycles = static_cast<size_t> (params->get_int("solver/pressure/n_cycle"));
    relaxs = static_cast<size_t> (params->get_int("solver/pressure/diffusion/n_relax"));

    m_dsign = -1.;
    m_w = 2. / 3.;
    m_w = params->get_real("solver/pressure/diffusion/w");

    // copies of out and b to prevent aliasing
    auto d_out = out->data;
    auto s_out = domain->get_size(out->GetLevel());
    auto t_out = out->GetType();
    auto d_b = b->data;
    auto s_b = domain->get_size(b->GetLevel());
    auto t_b = b->GetType();

    Field *out_err1 = new Field(t_out, 0.0);
    Field *out_tmp = new Field(t_out, 0.0);
    Field *b_res1 = new Field(t_b, 0.0);

    for (size_t i = 0; i < s_out; ++i) {
        out_err1->data[i] = d_out[i];
        out_tmp->data[i] = d_out[i];
    }
    for (size_t i = 0; i < s_b; ++i) {
        b_res1->data[i] = d_b[i];
    }

// residuum
    residuum1.push_back(b_res1);
    auto data_residuum1 = residuum1.back()->data;
    size_t bsize_residuum1 = domain->get_size(residuum1.back()->GetLevel());

#pragma acc enter data copyin(data_residuum1[:bsize_residuum1])

// error
    error1.push_back(out_err1);
    auto data_err1 = error1.back()->data;
    size_t bsize_err1 = domain->get_size(error1.back()->GetLevel());

#pragma acc enter data copyin(data_err1[:bsize_err1])

// temporal solution
    mg_temporal_solution.push_back(out_tmp);
    auto data_mg_temporal_solution = mg_temporal_solution.back()->data;
    size_t bsize_mg_temporal_solution = domain->get_size(mg_temporal_solution.back()->GetLevel());

#pragma acc enter data copyin(data_mg_temporal_solution[:bsize_mg_temporal_solution])


    //building Fields for level + sending to GPU
    // levels going up
    for (size_t i = 0; i < levels; ++i) {

        // build residuum0
        Field *r0 = new Field(FieldType::P, 0.0, i);
        residuum0.push_back(r0);

        auto data_residuum0 = r0->data;
        size_t bsize_residuum0 = domain->get_size(r0->GetLevel());

#pragma acc enter data copyin(data_residuum0[:bsize_residuum0])

        // build residuum1
        Field *r1 = new Field(FieldType::P, 0.0, i + 1);
        residuum1.push_back(r1);

        auto data_residuum1 = r1->data;
        size_t bsize_residuum1 = domain->get_size(r1->GetLevel());

#pragma acc enter data copyin(data_residuum1[:bsize_residuum1])

        //  build error1
        Field *e1 = new Field(FieldType::P, 0.0, i + 1);
        error1.push_back(e1);

        auto d_err1 = e1->data;
        size_t bsize_err1 = domain->get_size(e1->GetLevel());

#pragma acc enter data copyin(d_err1[:bsize_err1])

        // build mg_temporal_solution
        Field *mg = new Field(FieldType::P, 0.0, i + 1); //new field to prevent aliasing
        mg_temporal_solution.push_back(mg);

        auto data_mg_temporal_solution = mg->data;
        size_t bsize_mg_temporal_solution = domain->get_size(mg->GetLevel());

#pragma acc enter data copyin(data_mg_temporal_solution[:bsize_mg_temporal_solution])
    } // end level loop

//  build err0
    err0.resize(levels + 1);

    Field *e00 = new Field(FieldType::P, 0.0, error1[0]->GetLevel());

    err0[0] = e00;

    auto d_err00 = e00->data;
    size_t bsize_err00 = domain->get_size(e00->GetLevel());

#pragma acc enter data copyin(d_err00[:bsize_err00])

    //building Fields for level
    // levels going down
    for (size_t i = levels; i > 0; --i) {

        // build err0
        Field *e0 = new Field(FieldType::P, 0.0, error1[i - 1]->GetLevel());

        err0[i] = e0;

        auto d_err0 = err0[i]->data;
        size_t bsize_err0 = domain->get_size(err0[i]->GetLevel());

#pragma acc enter data copyin(d_err0[:bsize_err0])
    }
}

VCycleMG::~VCycleMG() {
    auto domain = Domain::getInstance();

    while (residuum0.size() > 0) {
        auto field = residuum0.back();
        auto data = residuum0.back()->data;
        size_t bsize = domain->get_size(residuum0.back()->GetLevel());
#pragma acc exit data delete(data[:bsize])
        delete field;
        field = nullptr;
        residuum0.pop_back();
    }
    while (residuum1.size() > 0) {
        auto field = residuum1.back();
        auto data = residuum1.back()->data;
        size_t bsize = domain->get_size(residuum1.back()->GetLevel());
#pragma acc exit data delete(data[:bsize])
        delete field;
        field = nullptr;
        residuum1.pop_back();
    }
    while (err0.size() > 0) {
        auto field = err0.back();
        auto data = err0.back()->data;
        size_t bsize = domain->get_size(err0.back()->GetLevel());
#pragma acc exit data delete(data[:bsize])
        delete field;
        field = nullptr;
        err0.pop_back();
    }
    while (error1.size() > 0) {
        auto field = error1.back();
        auto data = error1.back()->data;
        size_t bsize = domain->get_size(error1.back()->GetLevel());
#pragma acc exit data delete(data[:bsize])
        delete field;
        field = nullptr;
        error1.pop_back();
    }
    while (mg_temporal_solution.size() > 0) {
        auto field = mg_temporal_solution.back();
        auto data = mg_temporal_solution.back()->data;
        size_t bsize = domain->get_size(mg_temporal_solution.back()->GetLevel());
#pragma acc exit data delete(data[:bsize])
        delete field;
        field = nullptr;
        mg_temporal_solution.pop_back();
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
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    auto d_out = out->data;
    auto s_out = domain->get_size(out->GetLevel());
    auto d_b = b->data;
    auto s_b = domain->get_size(b->GetLevel());

    auto f_err1 = error1[0];
    auto f_mg_tmp = mg_temporal_solution[0];
    auto f_res1 = residuum1[0];

    auto d_err1 = error1[0]->data;
    auto d_mg_tmp = mg_temporal_solution[0]->data;
    auto d_res1 = residuum1[0]->data;

    auto s_err1 = domain->get_size(error1[0]->GetLevel());
    auto s_mg_tmp = domain->get_size(mg_temporal_solution[0]->GetLevel());
    auto s_res1 = domain->get_size(residuum1[0]->GetLevel());

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();
    auto bsize_i = boundary->getSize_innerList();

    // use iList on level 0, since update on level 0
#pragma acc kernels present(d_out[:s_out], d_b[:s_b], d_err1[:s_err1], \
                            d_mg_tmp[:s_mg_tmp], d_res1[:s_res1], d_iList[:bsize_i]) async
    {
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            d_err1[i] = d_out[i];
            d_mg_tmp[i] = d_out[i];
            d_res1[i] = d_b[i];
        }
        if (sync) {
#pragma acc wait
        }
    }//end data region
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
    auto params = Parameters::getInstance();
    auto domain = Domain::getInstance();
    const real dt = params->get_real("physical_parameters/dt");
    const size_t Nt = static_cast<size_t>(std::round(t / dt));

    const size_t set_relaxs = static_cast<size_t>(params->get_int("solver/pressure/diffusion/n_relax"));
    const size_t set_cycles = static_cast<size_t>(params->get_int("solver/pressure/n_cycle"));

    size_t act_cycles = 0;

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

        auto d_out = out->data;
        auto d_b = b->data;

        auto boundary = BoundaryController::getInstance();
        auto bsize_i = boundary->getSize_innerList();
        size_t *d_iList = boundary->get_innerList_level_joined();

        while (r > tol_res &&
                static_cast<int>(act_cycles) < max_cycles &&
                static_cast<int>(relaxs) < max_relaxs) {
            for (size_t i = 0; i < cycles; i++) {
                VCycleMultigrid(out, sync);
                act_cycles++;
            }
            sum = 0.;

            // calculate residuum in inner cells
#pragma acc parallel loop independent present(d_out[:bsize], d_b[:bsize], d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
                r = d_b[i] - (rdx2 * (d_out[i - 1] - 2 * d_out[i] + d_out[i + 1])\
 + rdy2 * (d_out[i - Nx] - 2 * d_out[i] + d_out[i + Nx])\
 + rdz2 * (d_out[i - Nx * Ny] - 2 * d_out[i] + d_out[i + Nx * Ny]));
                sum += r * r;
            }

#pragma acc wait

            r = sqrt(sum);

            relaxs += set_relaxs;
        } //end while
    } // end if

    else { // Nt > 1

        cycles = set_cycles;
        relaxs = set_relaxs;

        for (size_t i = 0; i < cycles; i++)
            VCycleMultigrid(out, sync);

    }

    if (sync) {
#pragma acc wait
    }
}

//==================================== VCycle =================================
// *****************************************************************************
/// \brief  Conducts the V-cycle Multigrid method
/// \param  out         pressure
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::VCycleMultigrid(Field *out, bool sync) {
    size_t max_level = levels;

    auto domain = Domain::getInstance();
    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t *d_bList = boundary->get_boundaryList_level_joined();

//===================== No refinement, when levels=0 =========//
    if (max_level == 0) {
        auto f_mg_tmpi = mg_temporal_solution[0];
        auto f_res1i = residuum1[0];
        auto d_mg_tmpi = mg_temporal_solution[0]->data;
        auto d_res1i = residuum1[0]->data;
        auto d_out = out->data;
        auto s_mg_tmpi = domain->get_size(mg_temporal_solution[0]->GetLevel());
        auto s_res1i = domain->get_size(residuum1[0]->GetLevel());
        auto s_out = domain->get_size();

#pragma acc data present(d_out[:s_out], d_mg_tmpi[:s_mg_tmpi], d_res1i[:s_res1i])
        {
            Solve(out, f_mg_tmpi, f_res1i, max_level, sync);
        }
        return;
    } //end if

//===================== levels going down ====================//
    for (size_t i = 0; i < max_level; ++i) {

        auto f_res0i = residuum0[i];
        auto f_err1i = error1[i];
        auto f_err1ip = error1[i + 1];
        auto f_mg_tmpi = mg_temporal_solution[i];
        auto f_res1i = residuum1[i];
        auto f_res1ip = residuum1[i + 1];

        auto d_res0i = residuum0[i]->data;
        auto d_err1i = error1[i]->data;
        auto d_err1ip = error1[i + 1]->data;
        auto d_mg_tmpi = mg_temporal_solution[i]->data;
        auto d_res1i = residuum1[i]->data;
        auto d_res1ip = residuum1[i + 1]->data;
        auto d_out = out->data;

        auto s_res0i = domain->get_size(residuum0[i]->GetLevel());
        auto s_err1i = domain->get_size(error1[i]->GetLevel());
        auto s_err1ip = domain->get_size(error1[i + 1]->GetLevel());
        auto s_mg_tmpi = domain->get_size(mg_temporal_solution[i]->GetLevel());
        auto s_res1i = domain->get_size(residuum1[i]->GetLevel());
        auto s_res1ip = domain->get_size(residuum1[i + 1]->GetLevel());
        auto s_out = domain->get_size(out->GetLevel());

        FieldType type_r0 = f_res0i->GetType();

        FieldType type_r1 = residuum1[i + 1]->GetType();

        //auto bsize_i = boundary->get_MGiList_Size(i+1);
        //auto bsize_b = boundary->get_MGbList_Size(i+1);

#pragma acc data present(    d_res0i[:s_res0i], d_err1i[:s_err1i], d_err1ip[:s_err1ip], \
                            d_mg_tmpi[:s_mg_tmpi], d_res1i[:s_res1i], d_res1ip[:s_res1ip], \
                            d_out[:s_out])
        {
            if (i == 0) // use p=out on finest grid
            {
                // smooth
                Smooth(out, f_mg_tmpi, f_res1i, i, sync);

                // calculate residuum
                Residuum(f_res0i, out, f_res1i, i, sync);
                boundary->applyBoundary(d_res0i, i, type_r0, sync); // for residuum0 only Dirichlet BC
            } else {
                // smooth
                Smooth(f_err1i, f_mg_tmpi, f_res1i, i, sync);

                // calculate residuum
                Residuum(f_res0i, f_err1i, f_res1i, i, sync);
                boundary->applyBoundary(d_res0i, i, type_r0, sync); // for residuum0 only Dirichlet BC
            }

            // restrict
            Restrict(f_res1ip, f_res0i, i, sync);
            boundary->applyBoundary(d_res1ip, i + 1, type_r1, sync); // for res only Dirichlet BC

            // set err to zero at next level

            // strides (since GPU needs joined list)
            // inner start/ end index of level i + 1
            size_t start_i = boundary->get_innerList_level_joined_start(i + 1);
            size_t end_i = boundary->get_innerList_level_joined_end(i + 1) + 1;
            // boundary start/ end index of level i + 1
            size_t start_b = boundary->get_boundaryList_level_joined_start(i + 1);
            size_t end_b = boundary->get_boundaryList_level_joined_end(i + 1) + 1;
            // inner
#pragma acc kernels present(d_err1ip[:s_err1ip], d_iList[start_i:(end_i-start_i)]) async
#pragma acc loop independent
            for (size_t j = start_i; j < end_i; ++j) {
                const size_t idx = d_iList[j];
                d_err1ip[idx] = 0.0;
            }

            //boundary
#pragma acc kernels present(d_err1ip[:s_err1ip], d_bList[start_b:(end_b-start_b)]) async
#pragma acc loop independent
            for (size_t j = start_b; j < end_b; ++j) {
                const size_t idx = d_bList[j];
                d_err1ip[idx] = 0.0;
            }
        } //end data  region
    } //for (levels going down)

//===================== levels going up ====================//
    for (size_t i = max_level; i > 0; --i) {

        auto f_err0i = err0[i];
        auto f_err1i = error1[i];
        auto f_err1im = error1[i - 1];
        auto f_mg_tmpim = mg_temporal_solution[i - 1];
        auto f_res1im = residuum1[i - 1];

        auto d_err0i = err0[i]->data;
        auto d_err1i = error1[i]->data;
        auto d_err1im = error1[i - 1]->data;
        auto d_mg_tmpim = mg_temporal_solution[i - 1]->data;
        auto d_res1im = residuum1[i - 1]->data;
        auto d_out = out->data;

        auto s_err0i = domain->get_size(err0[i]->GetLevel());
        auto s_err1i = domain->get_size(error1[i]->GetLevel());
        auto s_err1im = domain->get_size(error1[i - 1]->GetLevel());
        auto s_mg_tmpim = domain->get_size(mg_temporal_solution[i - 1]->GetLevel());
        auto s_res1im = domain->get_size(residuum1[i - 1]->GetLevel());
        auto s_out = domain->get_size(out->GetLevel());

        size_t Nx_e0 = domain->get_Nx(f_err0i->GetLevel());
        size_t Ny_e0 = domain->get_Ny(f_err0i->GetLevel());
        size_t Nz_e0 = domain->get_Nz(f_err0i->GetLevel());

        real dx_e0 = domain->get_dx(f_err0i->GetLevel());
        real dy_e0 = domain->get_dy(f_err0i->GetLevel());
        real dz_e0 = domain->get_dz(f_err0i->GetLevel());

        FieldType type_e0 = f_err0i->GetType();

        // inner start/ end index of level i - 1
        size_t start_i = boundary->get_innerList_level_joined_start(i - 1);
        size_t end_i = boundary->get_innerList_level_joined_end(i - 1) + 1;
        // boundary start/ end index of level i - 1
        size_t start_b = boundary->get_boundaryList_level_joined_start(i - 1);
        size_t end_b = boundary->get_boundaryList_level_joined_end(i - 1) + 1;

#pragma acc data present(d_err0i[:s_err0i], d_err1i[:s_err1i], d_err1im[:s_err1im], d_mg_tmpim[:s_mg_tmpim], d_res1im[:s_res1im], d_out[:s_out])
        {
            // prolongate
            Prolongate(f_err0i, f_err1i, i, sync);
            boundary->applyBoundary(d_err0i, i - 1, type_e0, sync); // for err0 only Dirichlet BC

            // correct
            if (i == 1) // use p=out on finest grid
            {
                // inner
#pragma acc kernels present(d_err0i[:s_err0i], d_out[:s_out], d_iList[start_i:(end_i-start_i)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_i; j < end_i; ++j) {
                        const size_t idx = d_iList[j];
                        d_out[idx] += d_err0i[idx];
                    }
                }
                // boundary
#pragma acc kernels present(d_err0i[:s_err0i], d_out[:s_out], d_bList[start_b:(end_b-start_b)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_b; j < end_b; ++j) {
                        const size_t idx = d_bList[j];
                        d_out[idx] += d_err0i[idx];
                    }
                }
                // smooth
                Smooth(out, f_mg_tmpim, f_res1im, i - 1, sync);
            } else {
                // correct
                // inner
#pragma acc kernels present(d_err0i[:s_err0i], d_err1im[:s_err1im], d_iList[start_i:(end_i-start_i)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_i; j < end_i; ++j) {
                        const size_t idx = d_iList[j];
                        d_err1im[idx] += d_err0i[idx];
                    }
                }
                // boundary
#pragma acc kernels present(d_err0i[:s_err0i], d_err1im[:s_err1im], d_bList[start_b:(end_b-start_b)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_b; j < end_b; ++j) {
                        const size_t idx = d_bList[j];
                        d_err1im[idx] += d_err0i[idx];
                    }
                }
                // smooth
                if (i - 1 == levels - 1) Solve(f_err1im, f_mg_tmpim, f_res1im, i - 1, sync);
                else Smooth(f_err1im, f_mg_tmpim, f_res1im, i - 1, sync); // for err only Dirichlet BC
            }
        } //end data region
    }//for (levels going up)

    // set boundaries

    auto d_out = out->data;
    FieldType type = out->GetType();

    boundary->applyBoundary(d_out, type, sync);

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
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    const real dx = domain->get_dx(out->GetLevel());
    const real dy = domain->get_dy(out->GetLevel());
    const real dz = domain->get_dz(out->GetLevel());

    auto d_out = out->data;
    auto d_tmp = tmp->data;
    auto d_b = b->data;

    auto params = Parameters::getInstance();

    size_t bsize = domain->get_size(out->GetLevel());
    FieldType type = out->GetType();

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t *d_bList = boundary->get_boundaryList_level_joined();

    // start/ end index of level
    // inner
    size_t start_i = boundary->get_innerList_level_joined_start(level);
    size_t end_i = boundary->get_innerList_level_joined_end(level) + 1;
    // boundary
    size_t start_b = boundary->get_boundaryList_level_joined_start(level);
    size_t end_b = boundary->get_boundaryList_level_joined_end(level) + 1;

    // apply boundary: at level 0 apply set BC; else use Dirichlet 0
    boundary->applyBoundary(d_out, level, type, sync);

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

#pragma acc data present(d_out[:bsize], d_tmp[:bsize])
    {
        // initialization
        // inner
#pragma acc kernels present(d_iList[start_i:(end_i-start_i)]) async
#pragma acc loop independent
        for (size_t j = start_i; j < end_i; ++j) {
            const size_t i = d_iList[j];
            d_tmp[i] = d_out[i];
        }

        // boundary
#pragma acc kernels present(d_bList[start_b:(end_b-start_b)]) async
#pragma acc loop independent
        for (size_t j = start_b; j < end_b; ++j) {
            const size_t i = d_bList[j];
            d_tmp[i] = d_out[i];
        }
    }//end acc data

    // Diffusion step
    std::string diffusionType = params->get("solver/pressure/diffusion/type");
    if (diffusionType == DiffusionMethods::Jacobi) {

#pragma acc data present(d_out[:bsize], d_tmp[:bsize], d_b[:bsize])
        {
            for (size_t i = 0; i < relaxs; i++) { // fixed iteration number as in xml
                JacobiDiffuse::JacobiStep(level, out, tmp, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
                boundary->applyBoundary(d_out, level, type, sync);

                std::swap(tmp->data, out->data);
                std::swap(d_tmp, d_out);
            }

            if (relaxs % 2 != 0) // swap necessary when odd number of iterations
            {
#pragma acc kernels present(d_out[:bsize], d_tmp[:bsize], d_iList[start_i:(end_i-start_i)]) async
#pragma acc loop independent
                // inner
                for (size_t j = start_i; j < end_i; ++j) {
                    const size_t i = d_iList[j];
                    d_out[i] = d_tmp[i];
                }
                // boundary
#pragma acc kernels present(d_out[:bsize], d_tmp[:bsize], d_bList[start_b:(end_b-start_b)]) async
#pragma acc loop independent
                for (size_t j = start_b; j < end_b; ++j) {
                    const size_t i = d_bList[j];
                    d_out[i] = d_tmp[i];
                }
            }
        }//end data region
    } else if (diffusionType == DiffusionMethods::ColoredGaussSeidel) {

#pragma acc data present(d_out[:bsize], d_tmp[:bsize], d_b[:bsize])
        {
            for (size_t i = 0; i < relaxs; i++) {
                ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(out, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
                boundary->applyBoundary(d_out, level, type, sync); // for res/err only Dirichlet BC
            }
        } //end data region
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
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(in->GetLevel());
    const size_t Ny = domain->get_Ny(in->GetLevel());

    const real dx = domain->get_dx(in->GetLevel());
    const real dy = domain->get_dy(in->GetLevel());
    const real dz = domain->get_dz(in->GetLevel());

    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    auto d_out = out->data;
    auto d_in = in->data;
    auto d_b = b->data;

    size_t bsize = domain->get_size(out->GetLevel());

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();

    // starts/ ends
    // inner
    size_t start_i = boundary->get_innerList_level_joined_start(level);
    size_t end_i = boundary->get_innerList_level_joined_end(level) + 1;
#pragma acc data present(d_b[:bsize], d_in[:bsize], d_out[:bsize], d_iList[start_i:(end_i-start_i)])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = start_i; j < end_i; ++j) {
            const size_t i = d_iList[j];
            d_out[i] = d_b[i] - (rdx2 * (d_in[i - 1] - 2 * d_in[i] + d_in[i + 1])\
 + rdy2 * (d_in[i - Nx] - 2 * d_in[i] + d_in[i + Nx])\
 + rdz2 * (d_in[i - Nx * Ny] - 2 * d_in[i] + d_in[i + Nx * Ny]));
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
    const size_t Nx = domain->get_Nx(out->GetLevel());
    const size_t Ny = domain->get_Ny(out->GetLevel());

    // fine grid
    const size_t nx = domain->get_Nx(in->GetLevel());
    const size_t ny = domain->get_Ny(in->GetLevel());

    auto d_out = out->data;
    auto d_in = in->data;

    size_t bsize_out = domain->get_size(out->GetLevel());
    size_t bsize_in = domain->get_size(in->GetLevel());

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_innerList_level_joined();

    // start/end
    // inner
    size_t start_i = boundary->get_innerList_level_joined_start(level + 1);
    size_t end_i = boundary->get_innerList_level_joined_end(level + 1) + 1;

#ifndef BENCHMARKING
    if (end_i == start_i)
        m_logger->warn("Be cautious: Obstacle might fill up inner cells completely in level {} with nx= {}!",
                level, domain->get_nx(out->GetLevel()));
#endif

    // average from eight neighboring cells
    // obstacles not used in fine grid, since coarse grid only obstacle if one of 8 fine grids was an obstacle,
    // thus if coarse cell inner cell, then surrounding fine cells also inner cells!
    size_t i, j, k;

#pragma acc data present(d_in[:bsize_in], d_out[:bsize_out], d_iList[start_i:(end_i-start_i)])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = start_i; l < end_i; ++l) {
            const size_t idx = d_iList[l];
            k = idx / (Nx * Ny);
            j = (idx - k * Nx * Ny) / Nx;
            i = idx - k * Nx * Ny - j * Nx;

            d_out[idx] = 0.125 * (d_in[IX(2 * i - 1, 2 * j - 1, 2 * k - 1, nx, ny)]\
 + d_in[IX(2 * i, 2 * j - 1, 2 * k - 1, nx, ny)]\
 + d_in[IX(2 * i - 1, 2 * j, 2 * k - 1, nx, ny)]\
 + d_in[IX(2 * i, 2 * j, 2 * k - 1, nx, ny)]\
 + d_in[IX(2 * i - 1, 2 * j - 1, 2 * k, nx, ny)]\
 + d_in[IX(2 * i, 2 * j - 1, 2 * k, nx, ny)]\
 + d_in[IX(2 * i - 1, 2 * j, 2 * k, nx, ny)]\
 + d_in[IX(2 * i, 2 * j, 2 * k, nx, ny)]);

        }

        if (sync) {
#pragma acc wait
        }
    }// end data region
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
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    // fine grid
    const size_t nx = domain->get_Nx(out->GetLevel());
    const size_t ny = domain->get_Ny(out->GetLevel());

    // coarse grid
    const size_t Nx = domain->get_Nx(in->GetLevel());
    const size_t Ny = domain->get_Ny(in->GetLevel());

    auto d_out = out->data;
    auto d_in = in->data;

    size_t bsize_out = domain->get_size(out->GetLevel());
    size_t bsize_in = domain->get_size(in->GetLevel());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();

    // start/end (going backwards)
    // inner
    size_t start_i = boundary->get_innerList_level_joined_start(level);
    size_t end_i = boundary->get_innerList_level_joined_end(level) + 1;

    // prolongate
    size_t i, j, k;

#pragma acc data present(d_in[:bsize_in], d_out[:bsize_out], d_iList[start_i:(end_i-start_i)])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = start_i; l < end_i; ++l) {
            const size_t idx = d_iList[l];
            k = idx / (Nx * Ny);
            j = (idx - k * Nx * Ny) / Nx;
            i = idx - k * Nx * Ny - j * Nx;

            d_out[IX(2 * i, 2 * j, 2 * k, nx, ny)] = 0.015625 * (27 * d_in[idx] + 9 * d_in[idx + 1] + 9 * d_in[idx + Nx] + 9 * d_in[idx + Nx * Ny]\
 + 3 * d_in[idx + 1 + Nx] + 3 * d_in[idx + 1 + Nx * Ny] + 3 * d_in[idx + Nx + Nx * Ny] + d_in[idx + 1 + Nx + Nx * Ny]);
            d_out[IX(2 * i, 2 * j, 2 * k - 1, nx, ny)] = 0.015625 * (27 * d_in[idx] + 9 * d_in[idx + 1] + 9 * d_in[idx + Nx] + 9 * d_in[idx - Nx * Ny]\
 + 3 * d_in[idx + 1 + Nx] + 3 * d_in[idx + 1 - Nx * Ny] + 3 * d_in[idx + Nx - Nx * Ny] + d_in[idx + 1 + Nx - Nx * Ny]);
            d_out[IX(2 * i, 2 * j - 1, 2 * k, nx, ny)] = 0.015625 * (27 * d_in[idx] + 9 * d_in[idx + 1] + 9 * d_in[idx - Nx] + 9 * d_in[idx + Nx * Ny]\
 + 3 * d_in[idx + 1 - Nx] + 3 * d_in[idx + 1 + Nx * Ny] + 3 * d_in[idx - Nx + Nx * Ny] + d_in[idx + 1 - Nx + Nx * Ny]);
            d_out[IX(2 * i, 2 * j - 1, 2 * k - 1, nx, ny)] = 0.015625 * (27 * d_in[idx] + 9 * d_in[idx + 1] + 9 * d_in[idx - Nx] + 9 * d_in[idx - Nx * Ny]\
 + 3 * d_in[idx + 1 - Nx] + 3 * d_in[idx + 1 - Nx * Ny] + 3 * d_in[idx - Nx - Nx * Ny] + d_in[idx + 1 - Nx - Nx * Ny]);
            d_out[IX(2 * i - 1, 2 * j, 2 * k, nx, ny)] = 0.015625 * (27 * d_in[idx] + 9 * d_in[idx - 1] + 9 * d_in[idx + Nx] + 9 * d_in[idx + Nx * Ny]\
 + 3 * d_in[idx - 1 + Nx] + 3 * d_in[idx - 1 + Nx * Ny] + 3 * d_in[idx + Nx + Nx * Ny] + d_in[idx - 1 + Nx + Nx * Ny]);
            d_out[IX(2 * i - 1, 2 * j, 2 * k - 1, nx, ny)] = 0.015625 * (27 * d_in[idx] + 9 * d_in[idx - 1] + 9 * d_in[idx + Nx] + 9 * d_in[idx - Nx * Ny]\
 + 3 * d_in[idx - 1 + Nx] + 3 * d_in[idx - 1 - Nx * Ny] + 3 * d_in[idx + Nx - Nx * Ny] + d_in[idx - 1 + Nx - Nx * Ny]);
            d_out[IX(2 * i - 1, 2 * j - 1, 2 * k, nx, ny)] = 0.015625 * (27 * d_in[idx] + 9 * d_in[idx - 1] + 9 * d_in[idx - Nx] + 9 * d_in[idx + Nx * Ny]\
 + 3 * d_in[idx - 1 - Nx] + 3 * d_in[idx - 1 + Nx * Ny] + 3 * d_in[idx - Nx + Nx * Ny] + d_in[idx - 1 - Nx + Nx * Ny]);
            d_out[IX(2 * i - 1, 2 * j - 1, 2 * k - 1, nx, ny)] = 0.015625 * (27 * d_in[idx] + 9 * d_in[idx - 1] + 9 * d_in[idx - Nx] + 9 * d_in[idx - Nx * Ny]\
 + 3 * d_in[idx - 1 - Nx] + 3 * d_in[idx - 1 - Nx * Ny] + 3 * d_in[idx - Nx - Nx * Ny] + d_in[idx - 1 - Nx - Nx * Ny]);
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
    const size_t Nx = domain->get_Nx(out->GetLevel());
    const size_t Ny = domain->get_Ny(out->GetLevel());
    const size_t Nz = domain->get_Nz(out->GetLevel());

    if (level < levels - 1) {
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

    const real dx = domain->get_dx(out->GetLevel());
    const real dy = domain->get_dy(out->GetLevel());
    const real dz = domain->get_dz(out->GetLevel());

    auto d_out = out->data;
    auto d_tmp = tmp->data;
    auto d_b = b->data;

    auto params = Parameters::getInstance();

    size_t bsize = domain->get_size(out->GetLevel());
    FieldType type = out->GetType();

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_innerList_level_joined();
    size_t *d_bList = boundary->get_boundaryList_level_joined();

    auto bsize_i = boundary->getSize_innerList_level_joined();
    auto bsize_b = boundary->getSize_boundaryList_level_joined();

    // start/ end
    // inner
    size_t start_i = boundary->get_innerList_level_joined_start(level);
    size_t end_i = boundary->get_innerList_level_joined_end(level) + 1;
    size_t start_b = boundary->get_boundaryList_level_joined_start(level);
    size_t end_b = boundary->get_boundaryList_level_joined_end(level) + 1;

    boundary->applyBoundary(d_out, level, type, sync);

    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    // preparation for diffusion step (alpha, beta)
    real alphaX = rdx2;
    real alphaY = rdy2;
    real alphaZ = rdz2;

    const real rbeta = 2. * (alphaX + alphaY + alphaZ);
    const real beta = 1. / rbeta;

#pragma acc data present(d_out[:bsize], d_tmp[:bsize])
    {
        // initialization
        // inner
#pragma acc kernels present(d_iList[start_i:(end_i-start_i)]) async
#pragma acc loop independent
        for (size_t j = start_i; j < end_i; ++j) {
            const size_t i = d_iList[j];
            d_tmp[i] = d_out[i];
        }

        // boundary
#pragma acc kernels present(d_bList[start_b:(end_b-start_b)]) async
#pragma acc loop independent
        for (size_t j = start_b; j < end_b; ++j) {
            const size_t i = d_bList[j];
            d_tmp[i] = d_out[i];
        }
    }  // end data region

// Diffusion step
    std::string diffusionType = params->get("solver/pressure/diffusion/type");
    if (diffusionType == DiffusionMethods::Jacobi) {
#pragma acc data present(d_out[:bsize], d_tmp[:bsize], d_b[:bsize])
        {
            size_t it = 0;
            const size_t max_it = static_cast<size_t> (params->get_int("solver/pressure/diffusion/max_solve"));
            const real tol_res = params->get_real("solver/pressure/diffusion/tol_res");
            real sum;
            real res = 10000.;

            while (res > tol_res && it < max_it) {
                JacobiDiffuse::JacobiStep(level, out, tmp, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
                boundary->applyBoundary(d_out, level, type, sync);

                sum = 0.;

#pragma acc parallel loop independent present(d_out[:bsize], d_tmp[:bsize], d_iList[:bsize_i]) async
                for (size_t j = start_i; j < end_i; ++j) {
                    const size_t i = d_iList[j];
                    res = d_b[i] - (rdx2 * (d_out[i - 1] - 2 * d_out[i] + d_out[i + 1])\
 + rdy2 * (d_out[i - Nx] - 2 * d_out[i] + d_out[i + Nx])\
 + rdz2 * (d_out[i - Nx * Ny] - 2 * d_out[i] + d_out[i + Nx * Ny])); //res = rbeta*(d_out[i] - d_tmp[i]);
                    sum += res * res;
                }
                // info: in nvvp profile 8byte size copy from to device to/from pageable due to sum!

#pragma acc wait
                res = sqrt(sum);

                it++;

                std::swap(tmp->data, out->data);
                std::swap(d_tmp, d_out);
            }  // end while

            if (it % 2 != 0) {  // swap necessary when odd number of iterations
#pragma acc kernels present(d_out[:bsize], d_tmp[:bsize], d_iList[start_i:(end_i-start_i)]) async
#pragma acc loop independent
                // inner
                for (size_t j = start_i; j < end_i; ++j) {
                    const size_t i = d_iList[j];
                    d_out[i] = d_tmp[i];
                }
                // boundary
#pragma acc kernels present(d_out[:bsize], d_tmp[:bsize], d_bList[start_b:(end_b-start_b)]) async
#pragma acc loop independent
                for (size_t j = start_b; j < end_b; ++j) {
                    const size_t i = d_bList[j];
                    d_out[i] = d_tmp[i];
                }
            }
        } //end data region
    } else if (diffusionType == DiffusionMethods::ColoredGaussSeidel) {
#pragma acc data present(d_out[:bsize], d_tmp[:bsize], d_b[:bsize])
        {
            size_t it = 0;
            const size_t max_it = static_cast<size_t>(params->get_int("solver/diffusion/max_iter"));
            const real tol_res = params->get_real("solver/diffusion/tol_res");
            real sum;
            real res = 10000.;

            while (res > tol_res && it < max_it) {
                ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(out, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
                boundary->applyBoundary(d_out, level, type, sync); // for res/err only Dirichlet BC

                sum = 0.;

#pragma acc parallel loop independent present(d_out[:bsize], d_tmp[:bsize], d_iList[:bsize_i]) async
                for (size_t j = start_i; j < end_i; ++j) {
                    const size_t i = d_iList[j];
                    res = d_b[i] - (rdx2 * (d_out[i - 1] - 2 * d_out[i] + d_out[i + 1])\
 + rdy2 * (d_out[i - Nx] - 2 * d_out[i] + d_out[i + Nx])\
 + rdz2 * (d_out[i - Nx * Ny] - 2 * d_out[i] + d_out[i + Nx * Ny])); //res = rbeta*(d_out[i] - d_tmp[i]);
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
