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
VCycleMG::VCycleMG(Field const &out, Field const &b) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif

    auto params = Parameters::getInstance();
    auto domain = Domain::getInstance();

    m_levels = domain->get_levels();
    m_cycles = params->get_int("solver/pressure/n_cycle");
    m_relaxs = params->get_int("solver/pressure/diffusion/n_relax");

    m_dsign = -1.;
    m_w = 2. / 3.;
    m_w = params->get_real("solver/pressure/diffusion/w");

    Field *out_err1 = new Field(out.get_type(), 0.0, 0, domain->get_size());
    Field *out_tmp = new Field(out.get_type(), 0.0, 0, domain->get_size());
    Field *b_res1 = new Field(b.get_type(), 0.0, 0, domain->get_size());

    out_err1->copy_data(out);
    out_tmp->copy_data(out);
    b_res1->copy_data(b);

    // residuum
    b_res1->copyin();
    residuum1.push_back(b_res1);

    // error
    out_err1->copyin();
    error1.push_back(out_err1);

    // temporal solution
    out_tmp->copyin();
    mg_temporal_solution.push_back(out_tmp);

    // building Fields for level + sending to GPU
    // levels going up
    for (size_t i = 0; i < m_levels; ++i) {

        // build residuum0
        Field *r0 = new Field(FieldType::P, 0.0, i, domain->get_size(i));
        r0->copyin();
        residuum0.push_back(r0);

        // build residuum1
        Field *r1 = new Field(FieldType::P, 0.0, i + 1, domain->get_size(i + 1));
        r1->copyin();
        residuum1.push_back(r1);

        //  build error1
        Field *e1 = new Field(FieldType::P, 0.0, i + 1, domain->get_size(i + 1));
        e1->copyin();
        error1.push_back(e1);

        // build mg_temporal_solution
        Field *mg = new Field(FieldType::P, 0.0, i + 1, domain->get_size(i + 1));  // new field to prevent aliasing
        mg->copyin();
        mg_temporal_solution.push_back(mg);
    }

    //  build err0
    err0.resize(m_levels + 1);

    auto level = error1[0]->get_level();
    Field *e00 = new Field(FieldType::P, 0.0, level, domain->get_size(level));
    e00->copyin();

    err0[0] = e00;

    // building Fields for level
    // levels going down
    for (int i = m_levels; i > 0; --i) {
        // build err0
        level = error1[i - 1]->get_level();
        Field *e0 = new Field(FieldType::P, 0.0, level, domain->get_size(level));
        e0->copyin();

        err0[i] = e0;
    }
}

VCycleMG::~VCycleMG() {
    for (size_t i = 0; i < residuum0.size(); i++) {
        delete residuum0[i];
    }
    for (size_t i = 0; i < residuum1.size(); i++) {
        delete residuum1[i];
    }

    for (size_t i = 0; i < err0.size(); i++) {
        delete err0[i];
    }
    for (size_t i = 0; i < error1.size(); i++) {
        delete error1[i];
    }

    for (size_t i = 0; i < mg_temporal_solution.size(); i++) {
        delete mg_temporal_solution[i];
    }
}

// =============================== Update ===============================
// ************************************************************************
/// \brief  Update input
/// \param  out     pressure
/// \param  b       rhs
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ************************************************************************
void VCycleMG::UpdateInput(Field &out, Field const &b, bool sync) {
    // local variables and parameters for GPU
    auto err1 = error1[0];
    auto mg_tmp = mg_temporal_solution[0];
    auto res1 = residuum1[0];

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_inner_list_level_joined();
    auto bsize_i = boundary->get_size_inner_list();

    // use iList on level 0, since update on level 0
#pragma acc kernels present(out, b, err1, mg_tmp, res1, d_iList[:bsize_i]) async
    {
#pragma acc loop independent
        for (size_t j = 0; j < bsize_i; ++j) {
            const size_t i = d_iList[j];
            (*err1)[i] = out[i];
            (*mg_tmp)[i] = out[i];
            (*res1)[i] = b[i];
        }
        if (sync) {
#pragma acc wait
        }
    }
}

//==================================== Pressure =================================
// *****************************************************************************
/// \brief  solves Poisson equation \f$ \nabla^2 p = rhs\f$ via geometric multigrid (VCycle)
/// \param  out         pressure
/// \param  b           rhs
/// \param  t           current time
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::pressure(Field &out, Field const &b, real t, bool sync) {
    // Update first
    UpdateInput(out, b);

    // solve more accurately, in first time step
    auto params = Parameters::getInstance();
    auto domain = Domain::getInstance();
    const real dt = params->get_real("physical_parameters/dt");
    const size_t Nt = static_cast<size_t>(std::round(t / dt));

    const int set_relaxs = params->get_int("solver/pressure/diffusion/n_relax");
    const int set_cycles = params->get_int("solver/pressure/n_cycle");

    size_t act_cycles = 0;

    if (Nt == 1) {
        const size_t max_cycles = params->get_int("solver/pressure/max_cycle");
        const size_t max_relaxs = params->get_int("solver/pressure/diffusion/max_solve");

        real r = 10000.;
        real sum = 0;
        const real tol_res = params->get_real("solver/pressure/tol_res");

        const size_t Nx = domain->get_Nx();
        const size_t Ny = domain->get_Ny();

        const real dx = domain->get_dx();
        const real dy = domain->get_dy();
        const real dz = domain->get_dz();

        const real rdx2 = 1. / (dx * dx);
        const real rdy2 = 1. / (dy * dy);
        const real rdz2 = 1. / (dz * dz);

        auto boundary = BoundaryController::getInstance();
        auto bsize_i = boundary->get_size_inner_list();
        size_t *d_iList = boundary->get_inner_list_level_joined();

        while (r > tol_res &&
                act_cycles < max_cycles &&
                m_relaxs < max_relaxs) {
            for (size_t i = 0; i < m_cycles; i++) {
                VCycleMultigrid(out, sync);
                act_cycles++;
            }
            sum = 0.;

            // calculate residuum in inner cells
#pragma acc parallel loop independent present(out, b, d_iList[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = d_iList[j];
                r = b[i] - (rdx2 * (out[i - 1] - 2 * out[i] + out[i + 1])
                     + rdy2 * (out[i - Nx] - 2 * out[i] + out[i + Nx])
                     + rdz2 * (out[i - Nx * Ny] - 2 * out[i] + out[i + Nx * Ny]));
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
/// \param  out         pressure
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::VCycleMultigrid(Field &out, bool sync) {
    size_t max_level = m_levels;
    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_inner_list_level_joined();
    size_t *d_bList = boundary->get_boundary_list_level_joined();

//===================== No refinement, when levels=0 =========//
    if (max_level == 0) {
        auto f_mg_tmpi = mg_temporal_solution[0];
        auto f_res1i = residuum1[0];

#pragma acc data present(out, mg_tmpi, res1i)
        {
            Solve(out, *f_mg_tmpi, *f_res1i, max_level, sync);
        }
        return;
    }

    //===================== levels going down ====================//
    for (size_t i = 0; i < max_level; ++i) {
        auto f_res0i = residuum0[i];
        auto f_err1i = error1[i];
        auto f_mg_tmpi = mg_temporal_solution[i];
        auto f_res1i = residuum1[i];
        auto f_res1ip = residuum1[i + 1];

        auto d_res0i = residuum0[i]->data;
        auto d_err1ip = error1[i + 1]->data;
        auto d_res1ip = residuum1[i + 1]->data;

        FieldType type_r0 = f_res0i->get_type();
        FieldType type_r1 = residuum1[i + 1]->get_type();

#pragma acc data present(res0i, err1i, err1ip, mg_tmpi, res1i, res1ip, out)
        {
            if (i == 0) {  // use p=out on finest grid
                // smooth
                Smooth(out, *f_mg_tmpi, *f_res1i, i, sync);

                // calculate residuum
                Residuum(*f_res0i, out, *f_res1i, i, sync);
                boundary->apply_boundary(d_res0i, i, type_r0, sync);  // for residuum0 only Dirichlet BC
            } else {
                // smooth
                Smooth(*f_err1i, *f_mg_tmpi, *f_res1i, i, sync);

                // calculate residuum
                Residuum(*f_res0i, *f_err1i, *f_res1i, i, sync);
                boundary->apply_boundary(d_res0i, i, type_r0, sync);  // for residuum0 only Dirichlet BC
            }

            // restrict
            Restrict(*f_res1ip, *f_res0i, i, sync);
            boundary->apply_boundary(d_res1ip, i + 1, type_r1, sync);  // for res only Dirichlet BC

            // set err to zero at next level

            // strides (since GPU needs joined list)
            // inner start/ end index of level i + 1
            size_t start_i = boundary->get_inner_list_level_joined_start(i + 1);
            size_t end_i = boundary->get_inner_list_level_joined_end(i + 1) + 1;
            // boundary start/ end index of level i + 1
            size_t start_b = boundary->get_boundary_list_level_joined_start(i + 1);
            size_t end_b = boundary->get_boundary_list_level_joined_end(i + 1) + 1;
            // inner
#pragma acc kernels present(err1ip, d_iList[start_i:(end_i-start_i)]) async
#pragma acc loop independent
            for (size_t j = start_i; j < end_i; ++j) {
                const size_t idx = d_iList[j];
                d_err1ip[idx] = 0.0;
            }

            // boundary
#pragma acc kernels present(err1ip, d_bList[start_b:(end_b-start_b)]) async
#pragma acc loop independent
            for (size_t j = start_b; j < end_b; ++j) {
                const size_t idx = d_bList[j];
                d_err1ip[idx] = 0.0;
            }
        }
    }

//===================== levels going up ====================//
    for (size_t i = max_level; i > 0; --i) {
        Field &f_err0i = *err0[i];
        Field &f_err1i = *error1[i];
        Field &f_err1im = *error1[i - 1];
        Field &f_mg_tmpim = *mg_temporal_solution[i - 1];
        Field &f_res1im = *residuum1[i - 1];

        FieldType type_e0 = f_err0i.get_type();

        // inner start/ end index of level i - 1
        size_t start_i = boundary->get_inner_list_level_joined_start(i - 1);
        size_t end_i = boundary->get_inner_list_level_joined_end(i - 1) + 1;
        // boundary start/ end index of level i - 1
        size_t start_b = boundary->get_boundary_list_level_joined_start(i - 1);
        size_t end_b = boundary->get_boundary_list_level_joined_end(i - 1) + 1;

#pragma acc data present(f_err0i, f_err1i, f_err1im, f_mg_tmpim, f_res1im, out)
        {
            // prolongate
            Prolongate(f_err0i, f_err1i, i, sync);
            boundary->apply_boundary(f_err0i.data, i - 1, type_e0, sync);  // for err0 only Dirichlet BC

            // correct
            if (i == 1) {  // use p=out on finest grid
                // inner
#pragma acc kernels present(f_err0i, out, d_iList[start_i:(end_i-start_i)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_i; j < end_i; ++j) {
                        const size_t idx = d_iList[j];
                        out[idx] += f_err0i[idx];
                    }
                }
                // boundary
#pragma acc kernels present(f_err0i, out, d_bList[start_b:(end_b-start_b)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_b; j < end_b; ++j) {
                        const size_t idx = d_bList[j];
                        out[idx] += f_err0i[idx];
                    }
                }
                // smooth
                Smooth(out, f_mg_tmpim, f_res1im, i - 1, sync);
            } else {
                // correct
                // inner
#pragma acc kernels present(f_err0i, f_err1im, d_iList[start_i:(end_i-start_i)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_i; j < end_i; ++j) {
                        const size_t idx = d_iList[j];
                        f_err1im[idx] += f_err0i[idx];
                    }
                }
                // boundary
#pragma acc kernels present(err0i, err1im, d_bList[start_b:(end_b-start_b)]) async
                {
#pragma acc loop independent
                    for (size_t j = start_b; j < end_b; ++j) {
                        const size_t idx = d_bList[j];
                        f_err1im[idx] += f_err0i[idx];
                    }
                }
                // smooth
                if (i == m_levels) {
                    Solve(f_err1im, f_mg_tmpim, f_res1im, i - 1, sync);
                } else {
                    Smooth(f_err1im, f_mg_tmpim, f_res1im, i - 1, sync);  // for err only Dirichlet BC
                }
            }
        }
    }

    // set boundaries

    auto d_out = out.data;
    FieldType type = out.get_type();

    boundary->apply_boundary(d_out, type, sync);

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
void VCycleMG::Smooth(Field &out, Field &tmp, Field const &b, size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    const real dx = domain->get_dx(out.get_level());
    const real dy = domain->get_dy(out.get_level());
    const real dz = domain->get_dz(out.get_level());

    auto params = Parameters::getInstance();

    FieldType type = out.get_type();

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_inner_list_level_joined();
    size_t *d_bList = boundary->get_boundary_list_level_joined();

    // start/ end index of level
    // inner
    size_t start_i = boundary->get_inner_list_level_joined_start(level);
    size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;
    // boundary
    size_t start_b = boundary->get_boundary_list_level_joined_start(level);
    size_t end_b = boundary->get_boundary_list_level_joined_end(level) + 1;

    // apply boundary: at level 0 apply set BC; else use Dirichlet 0
    boundary->apply_boundary(out.data, level, type, sync);

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

#pragma acc data present(out, tmp)
    {
        // initialization
        // inner
#pragma acc kernels present(d_iList[start_i:(end_i-start_i)]) async
#pragma acc loop independent
        for (size_t j = start_i; j < end_i; ++j) {
            const size_t i = d_iList[j];
            tmp[i] = out[i];
        }

        // boundary
#pragma acc kernels present(d_bList[start_b:(end_b-start_b)]) async
#pragma acc loop independent
        for (size_t j = start_b; j < end_b; ++j) {
            const size_t i = d_bList[j];
            tmp[i] = out[i];
        }
    }

    // Diffusion step
    std::string diffusionType = params->get("solver/pressure/diffusion/type");
    if (diffusionType == DiffusionMethods::Jacobi) {
#pragma acc data present(out, tmp, b)
        {
            for (size_t i = 0; i < m_relaxs; i++) {  // fixed iteration number as in xml
                JacobiDiffuse::JacobiStep(level,
                        out, tmp, b,
                        alphaX, alphaY, alphaZ,
                        beta, m_dsign, m_w, sync);
                boundary->apply_boundary(out.data, level, type, sync);

                Field::swap(tmp, out);
            }

            if (m_relaxs % 2 != 0) {  // swap necessary when odd number of iterations
#pragma acc kernels present(out, tmp, d_iList[start_i:(end_i-start_i)]) async
#pragma acc loop independent
                // inner
                for (size_t j = start_i; j < end_i; ++j) {
                    const size_t i = d_iList[j];
                    out[i] = tmp[i];
                }
                // boundary
#pragma acc kernels present(out, tmp, d_bList[start_b:(end_b-start_b)]) async
#pragma acc loop independent
                for (size_t j = start_b; j < end_b; ++j) {
                    const size_t i = d_bList[j];
                    out[i] = tmp[i];
                }
            }
        }
    } else if (diffusionType == DiffusionMethods::ColoredGaussSeidel) {
#pragma acc data present(out, tmp, b)
        {
            for (size_t i = 0; i < m_relaxs; i++) {
                ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(out, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
                boundary->apply_boundary(out.data, level, type, sync);  // for res/err only Dirichlet BC
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
void VCycleMG::Residuum(Field &out, Field const &in, Field const &b, size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(in.get_level());
    const size_t Ny = domain->get_Ny(in.get_level());

    const real dx = domain->get_dx(in.get_level());
    const real dy = domain->get_dy(in.get_level());
    const real dz = domain->get_dz(in.get_level());

    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_inner_list_level_joined();

    // starts/ ends
    // inner
    size_t start_i = boundary->get_inner_list_level_joined_start(level);
    size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;
#pragma acc data present(d_b[:bsize], d_in[:bsize], d_out[:bsize], d_iList[start_i:(end_i-start_i)])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = start_i; j < end_i; ++j) {
            const size_t i = d_iList[j];
            out[i] = b[i] - (rdx2 * (in[i - 1] - 2 * in[i] + in[i + 1])
                         + rdy2 * (in[i - Nx] - 2 * in[i] + in[i + Nx])
                         + rdz2 * (in[i - Nx * Ny] - 2 * in[i] + in[i + Nx * Ny]));
        }

        if (sync) {
#pragma acc wait
        }
    }
}

//================================== Restrict ===============================
// *****************************************************************************
/// \brief  Restricts field from fine grid to coarse grid via averaging
/// \param  out         output field (half the size of input vector)
/// \param  in          input field (on fine grid)
/// \param  level       Multigrid level
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::Restrict(Field &out, Field const &in, size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    // coarse grid
    const size_t Nx = domain->get_Nx(out.get_level());
    const size_t Ny = domain->get_Ny(out.get_level());

    // fine grid
    const size_t nx = domain->get_Nx(in.get_level());
    const size_t ny = domain->get_Ny(in.get_level());

    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_inner_list_level_joined();

    // start/end
    // inner
    size_t start_i = boundary->get_inner_list_level_joined_start(level + 1);
    size_t end_i = boundary->get_inner_list_level_joined_end(level + 1) + 1;

#ifndef BENCHMARKING
    if (end_i == start_i) {
        m_logger->warn("Be cautious: Obstacle might fill up inner cells completely in level {} with nx= {}!",
                       level, domain->get_nx(out.get_level()));
    }
#endif

    // average from eight neighboring cells
    // obstacles not used in fine grid, since coarse grid only obstacle if one of 8 fine grids was an obstacle,
    // thus if coarse cell inner cell, then surrounding fine cells also inner cells!
    size_t i, j, k;

#pragma acc data present(in, out, d_iList[start_i:(end_i-start_i)])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = start_i; l < end_i; ++l) {
            const size_t idx = d_iList[l];
            k = idx / (Nx * Ny);
            j = (idx - k * Nx * Ny) / Nx;
            i = idx - k * Nx * Ny - j * Nx;

            out[idx] = 0.125 * (in[IX(2 * i - 1, 2 * j - 1, 2 * k - 1, nx, ny)]
                         + in[IX(2 * i, 2 * j - 1, 2 * k - 1, nx, ny)]
                         + in[IX(2 * i - 1, 2 * j, 2 * k - 1, nx, ny)]
                         + in[IX(2 * i, 2 * j, 2 * k - 1, nx, ny)]
                         + in[IX(2 * i - 1, 2 * j - 1, 2 * k, nx, ny)]
                         + in[IX(2 * i, 2 * j - 1, 2 * k, nx, ny)]
                         + in[IX(2 * i - 1, 2 * j, 2 * k, nx, ny)]
                         + in[IX(2 * i, 2 * j, 2 * k, nx, ny)]);
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
void VCycleMG::Prolongate(Field &out, Field const &in, size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    // fine grid
    const size_t nx = domain->get_Nx(out.get_level());
    const size_t ny = domain->get_Ny(out.get_level());

    // coarse grid
    const size_t Nx = domain->get_Nx(in.get_level());
    const size_t Ny = domain->get_Ny(in.get_level());

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_inner_list_level_joined();

    // start/end (going backwards)
    // inner
    size_t start_i = boundary->get_inner_list_level_joined_start(level);
    size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;

    // prolongate
    size_t i, j, k;

#pragma acc data present(in, out, d_iList[start_i:(end_i-start_i)])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = start_i; l < end_i; ++l) {
            const size_t idx = d_iList[l];
            k = idx / (Nx * Ny);
            j = (idx - k * Nx * Ny) / Nx;
            i = idx - k * Nx * Ny - j * Nx;

            out[IX(2 * i, 2 * j, 2 * k, nx, ny)] = 0.015625 * (27 * in[idx] + 9 * in[idx + 1] + 9 * in[idx + Nx] + 9 * in[idx + Nx * Ny]\
 + 3 * in[idx + 1 + Nx] + 3 * in[idx + 1 + Nx * Ny] + 3 * in[idx + Nx + Nx * Ny] + in[idx + 1 + Nx + Nx * Ny]);
            out[IX(2 * i, 2 * j, 2 * k - 1, nx, ny)] = 0.015625 * (27 * in[idx] + 9 * in[idx + 1] + 9 * in[idx + Nx] + 9 * in[idx - Nx * Ny]\
 + 3 * in[idx + 1 + Nx] + 3 * in[idx + 1 - Nx * Ny] + 3 * in[idx + Nx - Nx * Ny] + in[idx + 1 + Nx - Nx * Ny]);
            out[IX(2 * i, 2 * j - 1, 2 * k, nx, ny)] = 0.015625 * (27 * in[idx] + 9 * in[idx + 1] + 9 * in[idx - Nx] + 9 * in[idx + Nx * Ny]\
 + 3 * in[idx + 1 - Nx] + 3 * in[idx + 1 + Nx * Ny] + 3 * in[idx - Nx + Nx * Ny] + in[idx + 1 - Nx + Nx * Ny]);
            out[IX(2 * i, 2 * j - 1, 2 * k - 1, nx, ny)] = 0.015625 * (27 * in[idx] + 9 * in[idx + 1] + 9 * in[idx - Nx] + 9 * in[idx - Nx * Ny]\
 + 3 * in[idx + 1 - Nx] + 3 * in[idx + 1 - Nx * Ny] + 3 * in[idx - Nx - Nx * Ny] + in[idx + 1 - Nx - Nx * Ny]);
            out[IX(2 * i - 1, 2 * j, 2 * k, nx, ny)] = 0.015625 * (27 * in[idx] + 9 * in[idx - 1] + 9 * in[idx + Nx] + 9 * in[idx + Nx * Ny]\
 + 3 * in[idx - 1 + Nx] + 3 * in[idx - 1 + Nx * Ny] + 3 * in[idx + Nx + Nx * Ny] + in[idx - 1 + Nx + Nx * Ny]);
            out[IX(2 * i - 1, 2 * j, 2 * k - 1, nx, ny)] = 0.015625 * (27 * in[idx] + 9 * in[idx - 1] + 9 * in[idx + Nx] + 9 * in[idx - Nx * Ny]\
 + 3 * in[idx - 1 + Nx] + 3 * in[idx - 1 - Nx * Ny] + 3 * in[idx + Nx - Nx * Ny] + in[idx - 1 + Nx - Nx * Ny]);
            out[IX(2 * i - 1, 2 * j - 1, 2 * k, nx, ny)] = 0.015625 * (27 * in[idx] + 9 * in[idx - 1] + 9 * in[idx - Nx] + 9 * in[idx + Nx * Ny]\
 + 3 * in[idx - 1 - Nx] + 3 * in[idx - 1 + Nx * Ny] + 3 * in[idx - Nx + Nx * Ny] + in[idx - 1 - Nx + Nx * Ny]);
            out[IX(2 * i - 1, 2 * j - 1, 2 * k - 1, nx, ny)] = 0.015625 * (27 * in[idx] + 9 * in[idx - 1] + 9 * in[idx - Nx] + 9 * in[idx - Nx * Ny]\
 + 3 * in[idx - 1 - Nx] + 3 * in[idx - 1 - Nx * Ny] + 3 * in[idx - Nx - Nx * Ny] + in[idx - 1 - Nx - Nx * Ny]);
        }

        if (sync) {
#pragma acc wait
        }
    }
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
void VCycleMG::Solve(Field &out, Field &tmp, Field const &b, size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(out.get_level());
    const size_t Ny = domain->get_Ny(out.get_level());

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

    const real dx = domain->get_dx(out.get_level());
    const real dy = domain->get_dy(out.get_level());
    const real dz = domain->get_dz(out.get_level());

    auto params = Parameters::getInstance();

    FieldType type = out.get_type();

    auto boundary = BoundaryController::getInstance();

    size_t *d_iList = boundary->get_inner_list_level_joined();
    size_t *d_bList = boundary->get_boundary_list_level_joined();

    auto bsize_i = boundary->get_size_inner_list_level_joined();
    auto bsize_b = boundary->get_size_boundary_list_level_joined();

    // start/ end
    // inner
    size_t start_i = boundary->get_inner_list_level_joined_start(level);
    size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;
    size_t start_b = boundary->get_boundary_list_level_joined_start(level);
    size_t end_b = boundary->get_boundary_list_level_joined_end(level) + 1;

    boundary->apply_boundary(out.data, level, type, sync);

    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    // preparation for diffusion step (alpha, beta)
    real alphaX = rdx2;
    real alphaY = rdy2;
    real alphaZ = rdz2;

    const real rbeta = 2. * (alphaX + alphaY + alphaZ);
    const real beta = 1. / rbeta;

#pragma acc data present(out, tmp)
    {
        // initialization
        // inner
#pragma acc kernels present(d_iList[start_i:(end_i-start_i)]) async
#pragma acc loop independent
        for (size_t j = start_i; j < end_i; ++j) {
            const size_t i = d_iList[j];
            tmp[i] = out[i];
        }

        // boundary
#pragma acc kernels present(d_bList[start_b:(end_b-start_b)]) async
#pragma acc loop independent
        for (size_t j = start_b; j < end_b; ++j) {
            const size_t i = d_bList[j];
            tmp[i] = out[i];
        }
    }  // end data region

// Diffusion step
    std::string diffusionType = params->get("solver/pressure/diffusion/type");
    if (diffusionType == DiffusionMethods::Jacobi) {
#pragma acc data present(out, tmp, b)
        {
            size_t it = 0;
            const size_t max_it = static_cast<size_t> (params->get_int("solver/pressure/diffusion/max_solve"));
            const real tol_res = params->get_real("solver/pressure/diffusion/tol_res");
            real sum;
            real res = 10000.;

            while (res > tol_res && it < max_it) {
                JacobiDiffuse::JacobiStep(level, out, tmp, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
                boundary->apply_boundary(out.data, level, type, sync);

                sum = 0.;

#pragma acc parallel loop independent present(out, tmp, d_iList[:bsize_i]) async
                for (size_t j = start_i; j < end_i; ++j) {
                    const size_t i = d_iList[j];
                    res = b[i] - (rdx2 * (out[i - 1] - 2 * out[i] + out[i + 1])\
 + rdy2 * (out[i - Nx] - 2 * out[i] + out[i + Nx])\
 + rdz2 * (out[i - Nx * Ny] - 2 * out[i] + out[i + Nx * Ny])); //res = rbeta*(out[i] - tmp[i]);
                    sum += res * res;
                }
                // info: in nvvp profile 8byte size copy from to device to/from pageable due to sum!

#pragma acc wait
                res = sqrt(sum);

                it++;

                Field::swap(tmp, out);
            }

            if (it % 2 != 0) {  // swap necessary when odd number of iterations
#pragma acc kernels present(out, tmp, d_iList[start_i:(end_i-start_i)]) async
#pragma acc loop independent
                // inner
                for (size_t j = start_i; j < end_i; ++j) {
                    const size_t i = d_iList[j];
                    out[i] = tmp[i];
                }
                // boundary
#pragma acc kernels present(out, tmp, d_bList[start_b:(end_b-start_b)]) async
#pragma acc loop independent
                for (size_t j = start_b; j < end_b; ++j) {
                    const size_t i = d_bList[j];
                    out[i] = tmp[i];
                }
            }
        } //end data region
    } else if (diffusionType == DiffusionMethods::ColoredGaussSeidel) {
#pragma acc data present(out, tmp, b)
        {
            size_t it = 0;
            const size_t max_it = static_cast<size_t>(params->get_int("solver/diffusion/max_iter"));
            const real tol_res = params->get_real("solver/diffusion/tol_res");
            real sum;
            real res = 10000.;

            while (res > tol_res && it < max_it) {
                ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(out, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
                boundary->apply_boundary(out.data, level, type, sync); // for res/err only Dirichlet BC

                sum = 0.;

#pragma acc parallel loop independent present(out, tmp, d_iList[:bsize_i]) async
                for (size_t j = start_i; j < end_i; ++j) {
                    const size_t i = d_iList[j];
                    res = b[i] - (rdx2 * (out[i - 1] - 2 * out[i] + out[i + 1])\
 + rdy2 * (out[i - Nx] - 2 * out[i] + out[i + Nx])\
 + rdz2 * (out[i - Nx * Ny] - 2 * out[i] + out[i + Nx * Ny])); //res = rbeta*(out[i] - tmp[i]);
                    sum += res * res;
                }
                // info: in nvvp profile 8byte size copy from to device to/from pageable due to sum!

#pragma acc wait
                res = sqrt(sum);
                it++;
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
