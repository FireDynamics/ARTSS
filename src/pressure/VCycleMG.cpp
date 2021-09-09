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


// =============================== Constructor ===============================
// *****************************************************************************
/// \brief  Constructor
/// \param  out     pressure
/// \param  b       rhs
// *****************************************************************************
VCycleMG::VCycleMG(Field const &out, Field const &b) :
        m_levels(Domain::getInstance()->get_levels()),
        m_n_cycle(Parameters::getInstance()->get_int("solver/pressure/n_cycle")),
        m_n_relax(Parameters::getInstance()->get_int("solver/pressure/diffusion/n_relax")),
        m_dt(Parameters::getInstance()->get_real("physical_parameters/dt")),
        m_w(Parameters::getInstance()->get_real("solver/pressure/diffusion/w"))
{
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    Parameters *params = Parameters::getInstance();

    std::string diffusion_type = params->get("solver/pressure/diffusion/type");
    if (diffusion_type == DiffusionMethods::Jacobi) {
        m_diffusion_max_iter = static_cast<size_t>(params->get_int("solver/pressure/diffusion/max_solve"));
        m_smooth_diffusion_function = &VCycleMG::call_smooth_jacobi;
        m_solve_diffusion_function = &VCycleMG::call_solve_jacobi;
    } else if (diffusion_type == DiffusionMethods::ColoredGaussSeidel) {
        m_diffusion_max_iter = static_cast<size_t>(params->get_int("solver/pressure/diffusion/max_iter"));
        m_smooth_diffusion_function = &VCycleMG::call_smooth_colored_gauss_seidel;
        m_solve_diffusion_function = &VCycleMG::call_solve_colored_gauss_seidel;
    } else {
#ifndef BENCHMARKING
        m_logger->critical("Diffusion method not yet implemented! Simulation stopped!");
#endif
        // TODO Error handling
    }
    m_diffusion_tol_res = params->get_real("solver/pressure/diffusion/tol_res");

    m_residuum0 = new Field*[m_levels];
    m_residuum1 = new Field*[m_levels + 1];
    m_error0 = new Field*[m_levels];
    m_error1 = new Field*[m_levels + 1];
    m_mg_temporal_solution = new Field*[m_levels + 1];

    // copies of out and b to prevent aliasing
    // residuum
    auto *b_res1 = new Field(b);
    b_res1->copyin();
    *(m_residuum1) = b_res1;

    // error
    auto *out_err1 = new Field(out);
    out_err1->copyin();
    *(m_error1) = out_err1;

    // temporal solution
    auto *out_tmp = new Field(out);
    out_tmp->copyin();
    *(m_mg_temporal_solution) = out_tmp;

    // building fields for level + sending to GPU
    // level going up
    for (size_t level = 0; level < m_levels; ++level) {
        // build m_residuum0
        auto *r0 = new Field(FieldType::P, 0.0, level);
        r0->copyin();
        *(m_residuum0 + level) = r0;

        // build m_residuum1
        auto *r1 = new Field(FieldType::P, 0.0, level + 1);
        r1->copyin();
        *(m_residuum1 + level + 1) = r1;

        //  build m_error0
        auto *e0 = new Field(FieldType::P, 0.0, level);
        e0->copyin();
        *(m_error0 + level) = e0;

        //  build m_error1
        auto *e1 = new Field(FieldType::P, 0.0, level + 1);
        e1->copyin();
        *(m_error1 + level + 1) = e1;

        // build m_mg_temporal_solution
        auto *mg = new Field(FieldType::P, 0.0, level + 1);
        mg->copyin();
        *(m_mg_temporal_solution + level + 1) = mg;
    }
}

VCycleMG::~VCycleMG() {
    for (size_t i = 0; i < m_levels; i++) {
        delete m_residuum0[i];
        delete m_error0[i];
    }
    for (size_t i = 0; i < m_levels + 1; i++) {
        delete m_residuum1[i];
        delete m_error1[i];
        delete m_mg_temporal_solution[i];
    }
}

// =============================== Update ===============================
// ************************************************************************
/// \brief  Update input
/// \param  out     pressure
/// \param  b       rhs
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ************************************************************************
void VCycleMG::UpdateInput(Field &out, const Field &b, bool sync) {
    m_error1[0]->copy_data(out);
    m_mg_temporal_solution[0]->copy_data(out);
    m_residuum1[0]->copy_data(b);

    if (sync) {
#pragma acc wait
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
    UpdateInput(out, b, sync);

    // solve more accurately, in first time step
    Parameters *params = Parameters::getInstance();
    Domain *domain = Domain::getInstance();
    const auto Nt = static_cast<size_t>(std::round(t / m_dt));
    size_t act_cycles = 0;

    if (Nt == 1) {
        const size_t max_cycles = params->get_int("solver/pressure/max_cycle");
        const size_t max_relaxs = params->get_int("solver/pressure/diffusion/max_solve");
        size_t set_relax = m_n_relax;
        real r = 10000.;
        real sum;
        const real tol_res = params->get_real("solver/pressure/tol_res");

        const size_t Nx = domain->get_Nx();
        const size_t Ny = domain->get_Ny();

        const real dx = domain->get_dx();
        const real dy = domain->get_dy();
        const real dz = domain->get_dz();

        const real rdx2 = 1. / (dx * dx);
        const real rdy2 = 1. / (dy * dy);
        const real rdz2 = 1. / (dz * dz);

        BoundaryController *boundary = BoundaryController::getInstance();
        const size_t bsize_i = boundary->get_size_inner_list();
        size_t *data_inner_list = boundary->get_inner_list_level_joined();

        const size_t neighbour_i = 1;
        const size_t neighbour_j = Nx;
        const size_t neighbour_k = Nx * Ny;
        while (r > tol_res &&
                act_cycles < max_cycles &&
                set_relax < max_relaxs) {
            for (int i = 0; i < m_n_cycle; i++) {
                VCycleMultigrid(out, sync);
                act_cycles++;
            }
            sum = 0.;

            // calculate residuum in inner cells
#pragma acc parallel loop independent present(out, b, data_inner_list[:bsize_i]) async
            for (size_t j = 0; j < bsize_i; ++j) {
                const size_t i = data_inner_list[j];
                r = b[i] - (rdx2 * (out[i - neighbour_i] - 2 * out[i] + out[i + neighbour_i])
                         +  rdy2 * (out[i - neighbour_j] - 2 * out[i] + out[i + neighbour_j])
                         +  rdz2 * (out[i - neighbour_k] - 2 * out[i] + out[i + neighbour_k]));
                sum += r * r;
            }
#pragma acc wait
            r = sqrt(sum);
            set_relax += m_n_relax;
        }
    } else {  // Nt > 1
        for (int i = 0; i < m_n_cycle; i++) {
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
    auto boundary = BoundaryController::getInstance();
//===================== No refinement, when max_level = 0 =========//
    if (m_levels == 0) {
        Field *field_mg_temporal_level = m_mg_temporal_solution[0];
        Field *field_residuum1_level = m_residuum1[0];
#pragma acc data present(out, field_mg_temporal_level, field_residuum1_level)
        {
            Solve(out, *field_mg_temporal_level, *field_residuum1_level, m_levels, sync);
        }
    } else {
        m_error1[0]->copy_data(out);
        //===================== levels going down ====================//
        for (size_t level = 0; level < m_levels; ++level) {
            Field *field_residuum0_level = *(m_residuum0 + level);
            Field *field_error1_level = *(m_error1 + level);
            Field *field_error1_level_plus_1 = *(m_error1 + level + 1);
            Field *field_mg_temporal_level = *(m_mg_temporal_solution + level);
            Field *field_residuum1_level = *(m_residuum1 + level);
            Field *field_residuum1_level_plus_1 = *(m_residuum1 + level + 1);

            real *data_residuum0_level = field_residuum0_level->data;
            real *data_residuum1_level_plus_1 = field_residuum1_level_plus_1->data;

            FieldType type_r0 = field_residuum0_level->get_type();
            FieldType type_r1 = field_residuum1_level_plus_1->get_type();

#pragma acc data present(residuum0_level, error1_level, error1_level_plus_1, \
                         mg_temporal_level, residuum1_level, residuum1_level_plus_1, \
                         out)
            {
                // smooth
                Smooth(*field_error1_level, *field_mg_temporal_level, *field_residuum1_level, level, sync);

                // calculate residuum
                Residuum(*field_residuum0_level, *field_error1_level, *field_residuum1_level, level, sync);
                // for m_residuum0 only Dirichlet BC
                boundary->apply_boundary(data_residuum0_level, level, type_r0, sync);

                // restrict
                Restrict(*field_residuum1_level_plus_1, *field_residuum0_level, level, sync);
                // for res only Dirichlet BC
                boundary->apply_boundary(data_residuum1_level_plus_1, level + 1, type_r1, sync);

                // set err to zero at next level
                field_error1_level_plus_1->set_value(0);
            }
        }

        //===================== levels going up ====================//
        for (size_t level = m_levels; level > 0; --level) {
            Field *field_error0_level = *(m_error0 + level - 1);
            Field *field_error1_level = *(m_error1 + level);
            Field *field_error1_level_minus_1 = *(m_error1 + level - 1);
            Field *field_mg_temporal_solution_level_minus_1 = *(m_mg_temporal_solution + level - 1);
            Field *field_residuum1_level_minus_1 = *(m_residuum1 + level - 1);

            real *data_error0_level = field_error0_level->data;

            FieldType type_e0 = field_error0_level->get_type();

#pragma acc data present(error0_level, error1_level, error1_level_minus_1, mg_temporal_solution_level_minus_1, residuum1_level_minus_1, out)
            {
                // prolongate
                Prolongate(*field_error0_level, *field_error1_level, level, sync);
                // for m_error0 only Dirichlet BC
                boundary->apply_boundary(data_error0_level, level - 1, type_e0, sync);

                // correct
                *field_error1_level_minus_1 += *field_error0_level;
                // smooth
                if (level == m_levels) {
                    Solve(*field_error1_level_minus_1, *field_mg_temporal_solution_level_minus_1, *field_residuum1_level_minus_1, level - 1, sync);
                } else {
                    // for err only Dirichlet BC
                    Smooth(*field_error1_level_minus_1, *field_mg_temporal_solution_level_minus_1, *field_residuum1_level_minus_1, level - 1, sync);
                }
            }
        }
        // set boundaries
        real *data_out = out.data;
        FieldType type = out.get_type();
        boundary->apply_boundary(data_out, type, sync);

        if (sync) {
#pragma acc wait
        }
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
void VCycleMG::Smooth(Field &out, Field &tmp, Field const &b, const size_t level, bool sync) {
    real *data_out = out.data;
    FieldType type = out.get_type();

    BoundaryController *boundary = BoundaryController::getInstance();
    // apply boundary: at level 0 apply set BC; else use Dirichlet 0
    boundary->apply_boundary(data_out, level, type, sync);

    // Diffusion step
    (this->*m_smooth_diffusion_function)(out, tmp, b, level, sync);

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
void VCycleMG::Residuum(Field &out, Field const &in, Field const &b, const size_t level, bool sync) {
    Domain *domain = Domain::getInstance();

    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(level);
    const size_t Ny = domain->get_Ny(level);

    const real dx = domain->get_dx(level);
    const real dy = domain->get_dy(level);
    const real dz = domain->get_dz(level);

    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    BoundaryController *boundary = BoundaryController::getInstance();
    size_t *data_inner_list = boundary->get_inner_list_level_joined();

    // starts/ends
    // inner
    const size_t start_i = boundary->get_inner_list_level_joined_start(level);
    const size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;

    // neighbour cells, i/j/k represent the directions
    const size_t neighbour_cell_i = 1;
    const size_t neighbour_cell_j = Nx;
    const size_t neighbour_cell_k = Nx * Ny;
#pragma acc data present(b, in, out, inner_list)
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t j = start_i; j < end_i; ++j) {
            const size_t i = data_inner_list[j];
            real weighted = (rdx2 * (in[i - neighbour_cell_i] - 2 * in[i] + in[i + neighbour_cell_i])
                           + rdy2 * (in[i - neighbour_cell_j] - 2 * in[i] + in[i + neighbour_cell_j])
                           + rdz2 * (in[i - neighbour_cell_k] - 2 * in[i] + in[i + neighbour_cell_k]));
            out[i] = b[i] - weighted;
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
void VCycleMG::Restrict(Field &out, Field const &in, const size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    // coarse grid
    const size_t Nx_coarse = domain->get_Nx(out.get_level());
    const size_t Ny_coarse = domain->get_Ny(out.get_level());

    // fine grid
    const size_t Nx_fine = domain->get_Nx(in.get_level());
    const size_t Ny_fine = domain->get_Ny(in.get_level());

    auto boundary = BoundaryController::getInstance();
    size_t *data_inner_list = boundary->get_inner_list_level_joined();

    // start/end
    // inner
    const size_t start_i = boundary->get_inner_list_level_joined_start(level + 1);
    const size_t end_i = boundary->get_inner_list_level_joined_end(level + 1) + 1;

#ifndef BENCHMARKING
    if (end_i == start_i) {
        m_logger->warn("Be cautious: Obstacle might fill up inner cells completely in level {} with Nx_fine= {}!",
                       level + 1, domain->get_nx(out.get_level()));
    }
#endif

    // neighbour cells, i/j/k represent the directions
    const size_t neighbour_cell_i = 1;
    const size_t neighbour_cell_j = 1;
    const size_t neighbour_cell_k = 1;
    // average from eight neighboring cells
    // obstacles not used in fine grid, since coarse grid only obstacle if one of 8 fine grids was an obstacle,
    // thus if coarse cell inner cell, then surrounding fine cells also inner cells!
#pragma acc data present(in, out, data_inner_list[start_i:(end_i-start_i)])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = start_i; l < end_i; ++l) {
            const size_t idx = data_inner_list[l];
            const size_t k = getCoordinateK(idx, Nx_coarse, Ny_coarse);
            const size_t j = getCoordinateJ(idx, Nx_coarse, Ny_coarse, k);
            const size_t i = getCoordinateI(idx, Nx_coarse, Ny_coarse, j, k);

            out[idx] = 0.125
                   * (in[IX(2 * i - neighbour_cell_i, 2 * j - neighbour_cell_j, 2 * k - neighbour_cell_k, Nx_fine, Ny_fine)]
                   +  in[IX(2 * i,                    2 * j - neighbour_cell_j, 2 * k - neighbour_cell_k, Nx_fine, Ny_fine)]
                   +  in[IX(2 * i - neighbour_cell_i, 2 * j,                    2 * k - neighbour_cell_k, Nx_fine, Ny_fine)]
                   +  in[IX(2 * i,                    2 * j,                    2 * k - neighbour_cell_k, Nx_fine, Ny_fine)]
                   +  in[IX(2 * i - neighbour_cell_i, 2 * j - neighbour_cell_j, 2 * k,                    Nx_fine, Ny_fine)]
                   +  in[IX(2 * i,                    2 * j - neighbour_cell_j, 2 * k,                    Nx_fine, Ny_fine)]
                   +  in[IX(2 * i - neighbour_cell_i, 2 * j,                    2 * k,                    Nx_fine, Ny_fine)]
                   +  in[IX(2 * i,                    2 * j,                    2 * k,                    Nx_fine, Ny_fine)]);

        }

        if (sync) {
#pragma acc wait
        }
    }
}

//================================== Prolongate ================================
// *****************************************************************************
/// \brief  Prolongates field from coarse grid to fine grid (trilinear interpolation)
/// \param  out         output field (real the size of input vector)
/// \param  in          input field (on coarse grid)
/// \param  level       Multigrid level
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::Prolongate(Field &out, Field const &in, const size_t level, bool sync) {
    Domain *domain = Domain::getInstance();

    // local variables and parameters for GPU
    // fine grid
    const size_t Nx_fine = domain->get_Nx(out.get_level());
    const size_t Ny_fine = domain->get_Ny(out.get_level());

    // coarse grid
    const size_t Nx_coarse = domain->get_Nx(in.get_level());
    const size_t Ny_coarse = domain->get_Ny(in.get_level());

    BoundaryController *boundary = BoundaryController::getInstance();
    size_t *data_inner_list = boundary->get_inner_list_level_joined();

    // start/end (going backwards)
    // inner
    const size_t start_i = boundary->get_inner_list_level_joined_start(level);
    const size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;

    // neighbour cells, i/j/k represent the directions
    const size_t neighbour_cell_i = 1;
    const size_t neighbour_cell_j = Nx_coarse;
    const size_t neighbour_cell_k = Nx_coarse * Ny_coarse;

    // prolongate
#pragma acc data present(in, out, data_inner_list[start_i:(end_i-start_i)])
    {
#pragma acc kernels async
#pragma acc loop independent
        for (size_t l = start_i; l < end_i; ++l) {
            const size_t idx = data_inner_list[l];
            const size_t k = getCoordinateK(idx, Nx_coarse, Ny_coarse);
            const size_t j = getCoordinateJ(idx, Nx_coarse, Ny_coarse, k);
            const size_t i = getCoordinateI(idx, Nx_coarse, Ny_coarse, j, k);

            // x + 1 and y + 1 and z + 1
            out[IX(2 * i, 2 * j, 2 * k, Nx_fine, Ny_fine)] =
                    0.015625 * (27 * in[idx]
                               + 9 * in[idx + neighbour_cell_i]
                               + 9 * in[idx + neighbour_cell_j]
                               + 9 * in[idx + neighbour_cell_k]
                               + 3 * in[idx + neighbour_cell_i + neighbour_cell_j]
                               + 3 * in[idx + neighbour_cell_i + neighbour_cell_k]
                               + 3 * in[idx + neighbour_cell_j + neighbour_cell_k]
                               +     in[idx + neighbour_cell_i + neighbour_cell_j + neighbour_cell_k]);
            // x + 1 and y + 1 and z - 1
            out[IX(2 * i, 2 * j, 2 * k - 1, Nx_fine, Ny_fine)] =
                    0.015625 * (27 * in[idx]
                               + 9 * in[idx + neighbour_cell_i]
                               + 9 * in[idx + neighbour_cell_j]
                               + 9 * in[idx - neighbour_cell_k]
                               + 3 * in[idx + neighbour_cell_i + neighbour_cell_j]
                               + 3 * in[idx + neighbour_cell_i - neighbour_cell_k]
                               + 3 * in[idx + neighbour_cell_j - neighbour_cell_k]
                               +     in[idx + neighbour_cell_i + neighbour_cell_j - neighbour_cell_k]);
            // x + 1 and y - 1 and k + 1
            out[IX(2 * i, 2 * j - 1, 2 * k, Nx_fine, Ny_fine)] =
                    0.015625 * (27 * in[idx]
                               + 9 * in[idx + neighbour_cell_i]
                               + 9 * in[idx - neighbour_cell_j]
                               + 9 * in[idx + neighbour_cell_k]
                               + 3 * in[idx + neighbour_cell_i - neighbour_cell_j]
                               + 3 * in[idx + neighbour_cell_i + neighbour_cell_k]
                               + 3 * in[idx - neighbour_cell_j + neighbour_cell_k]
                               +     in[idx + neighbour_cell_i - neighbour_cell_j + neighbour_cell_k]);
            // x + 1 and y - 1 and z - 1
            out[IX(2 * i, 2 * j - 1, 2 * k - 1, Nx_fine, Ny_fine)] =
                    0.015625 * (27 * in[idx]
                               + 9 * in[idx + neighbour_cell_i]
                               + 9 * in[idx - neighbour_cell_j]
                               + 9 * in[idx - neighbour_cell_k]
                               + 3 * in[idx + neighbour_cell_i - neighbour_cell_j]
                               + 3 * in[idx + neighbour_cell_i - neighbour_cell_k]
                               + 3 * in[idx - neighbour_cell_j - neighbour_cell_k]
                               +     in[idx + neighbour_cell_i - neighbour_cell_j - neighbour_cell_k]);
            // x - 1 and y + 1 and z + 1
            out[IX(2 * i - 1, 2 * j, 2 * k, Nx_fine, Ny_fine)] =
                    0.015625 * (27 * in[idx]
                               + 9 * in[idx - neighbour_cell_i]
                               + 9 * in[idx + neighbour_cell_j]
                               + 9 * in[idx + neighbour_cell_k]
                               + 3 * in[idx - neighbour_cell_i + neighbour_cell_j]
                               + 3 * in[idx - neighbour_cell_i + neighbour_cell_k]
                               + 3 * in[idx + neighbour_cell_j + neighbour_cell_k]
                               +     in[idx - neighbour_cell_i + neighbour_cell_j + neighbour_cell_k]);
            // x - 1 and j + 1 and z - 1
            out[IX(2 * i - 1, 2 * j, 2 * k - 1, Nx_fine, Ny_fine)] =
                    0.015625 * (27 * in[idx]
                               + 9 * in[idx - neighbour_cell_i]
                               + 9 * in[idx + neighbour_cell_j]
                               + 9 * in[idx - neighbour_cell_k]
                               + 3 * in[idx - neighbour_cell_i + neighbour_cell_j]
                               + 3 * in[idx - neighbour_cell_i - neighbour_cell_k]
                               + 3 * in[idx + neighbour_cell_j - neighbour_cell_k]
                               +     in[idx - neighbour_cell_i + neighbour_cell_j - neighbour_cell_k]);
            // x - 1 and j - 1 and z + 1
            out[IX(2 * i - 1, 2 * j - 1, 2 * k, Nx_fine, Ny_fine)] =
                    0.015625 * (27 * in[idx]
                               + 9 * in[idx - neighbour_cell_i]
                               + 9 * in[idx - neighbour_cell_j]
                               + 9 * in[idx + neighbour_cell_k]
                               + 3 * in[idx - neighbour_cell_i - neighbour_cell_j]
                               + 3 * in[idx - neighbour_cell_i + neighbour_cell_k]
                               + 3 * in[idx - neighbour_cell_j + neighbour_cell_k]
                               +     in[idx - neighbour_cell_i - neighbour_cell_j + neighbour_cell_k]);
            // x - 1 and y - 1 and z - 1
            out[IX(2 * i - 1, 2 * j - 1, 2 * k - 1, Nx_fine, Ny_fine)] =
                    0.015625 * (27 * in[idx]
                               + 9 * in[idx - neighbour_cell_i]
                               + 9 * in[idx - neighbour_cell_j]
                               + 9 * in[idx - neighbour_cell_k]
                               + 3 * in[idx - neighbour_cell_i - neighbour_cell_j]
                               + 3 * in[idx - neighbour_cell_i - neighbour_cell_k]
                               + 3 * in[idx - neighbour_cell_j - neighbour_cell_k]
                               +     in[idx - neighbour_cell_i - neighbour_cell_j - neighbour_cell_k]);
        }
        if (sync) {
#pragma acc wait
        }
    }
}

//==================================== Smooth ==================================
// *****************************************************************************
/// \brief  Solves Ax = b on the lowest level (m_levels - 1) using a diffusion method
/// \param  out         output field (in size of input )
/// \param  tmp         temporary field for JacobiStep
/// \param  b           right hand side
/// \param  level       Multigrid level
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// *****************************************************************************
void VCycleMG::Solve(Field &out, Field &tmp, Field const &b, const size_t level, bool sync) {
    auto domain = Domain::getInstance();

#ifndef BENCHMARKING
    if (level < m_levels - 1) {
        m_logger->warn("Trying to solve on level {}, but should be {}", level, m_levels - 1);
        return;
        // TODO Error handling
    }
#endif

    const size_t Nx = domain->get_Nx(out.get_level());
    const size_t Ny = domain->get_Ny(out.get_level());

    if (Nx <= 4 && Ny <= 4) {
#ifndef BENCHMARKING
        m_logger->warn(" Grid is too coarse with Nx={} and Ny={}. Just smooth here", Nx, Ny);
#endif
        Smooth(out, tmp, b, level, sync);
        return;
        // TODO Error handling
    }

// Diffusion step
    (this->*m_solve_diffusion_function)(out, tmp, b, level, sync);
    if (sync) {
#pragma acc wait
    }
}

void VCycleMG::call_smooth_colored_gauss_seidel(Field &out, Field &tmp, Field const &b, const size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    const real dx = domain->get_dx(level);
    const real dy = domain->get_dy(level);
    const real dz = domain->get_dz(level);

    FieldType type = out.get_type();
    real *data_out = out.data;

    BoundaryController *boundary = BoundaryController::getInstance();
    boundary->apply_boundary(data_out, level, type, sync);

    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    // preparation for diffusion step (alpha, beta)
    const real alphaX = rdx2;
    const real alphaY = rdy2;
    const real alphaZ = rdz2;

    const real rbeta = 2. * (alphaX + alphaY + alphaZ);
    const real beta = 1. / rbeta;

    tmp.copy_data(out);
#pragma acc data present(out, tmp, b)
    {
        for (int i = 0; i < m_n_relax; i++) {
            ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(out, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
            // for res/err only Dirichlet BC
            boundary->apply_boundary(data_out, level, type, sync);
        }
    }

}
void VCycleMG::call_solve_colored_gauss_seidel(Field &out, Field &tmp, Field const &b, const size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(level);
    const size_t Ny = domain->get_Ny(level);

    const real dx = domain->get_dx(level);
    const real dy = domain->get_dy(level);
    const real dz = domain->get_dz(level);

    real *data_out = out.data;
    FieldType type = out.get_type();

    BoundaryController *boundary = BoundaryController::getInstance();

    size_t *data_inner_list = boundary->get_inner_list_level_joined();
    const size_t bsize_i = boundary->get_size_inner_list_level_joined();

    // start/end
    // inner
    const size_t start_i = boundary->get_inner_list_level_joined_start(level);
    const size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;

    boundary->apply_boundary(data_out, level, type, sync);

    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    // preparation for diffusion step (alpha, beta)
    const real alphaX = rdx2;
    const real alphaY = rdy2;
    const real alphaZ = rdz2;

    const real rbeta = 2. * (alphaX + alphaY + alphaZ);
    const real beta = 1. / rbeta;

    const size_t max_it = m_diffusion_max_iter;
    const real tol_res = m_diffusion_tol_res;

    tmp.copy_data(out);
#pragma acc data present(out, tmp, b)
    {
        size_t it = 0;
        real sum;
        real res = 10000.;

        const size_t neighbour_i = 1;
        const size_t neighbour_j = Nx;
        const size_t neighbour_k = Nx * Ny;
        while (res > tol_res && it < max_it) {
            ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(out, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
            // for res/err only Dirichlet BC
            boundary->apply_boundary(data_out, level, type, sync);

            sum = 0.;
#pragma acc parallel loop independent present(out, tmp, data_inner_list[:bsize_i]) async
            for (size_t j = start_i; j < end_i; ++j) {
                const size_t i = data_inner_list[j];
                res = b[i] - (rdx2 * (data_out[i - neighbour_i] - 2 * data_out[i] + data_out[i + neighbour_i])
                           +  rdy2 * (data_out[i - neighbour_j] - 2 * data_out[i] + data_out[i + neighbour_j])
                           +  rdz2 * (data_out[i - neighbour_k] - 2 * data_out[i] + data_out[i + neighbour_k]));  //res = rbeta*(data_out[i] - data_tmp[i]);
                sum += res * res;
            }
            // info: in nvvp profile 8byte size copy from to device to/from pageable due to sum!

#pragma acc wait
            res = sqrt(sum);
            it++;
        }  // end while
    }  // end data region

}
void VCycleMG::call_smooth_jacobi(Field &out, Field &tmp, Field const &b, const size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    const real dx = domain->get_dx(level);
    const real dy = domain->get_dy(level);
    const real dz = domain->get_dz(level);

    real *data_out = out.data;
    real *data_tmp = tmp.data;

    FieldType type = out.get_type();

    BoundaryController *boundary = BoundaryController::getInstance();
    boundary->apply_boundary(data_out, level, type, sync);

    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    // preparation for diffusion step (alpha, beta)
    const real alphaX = rdx2;
    const real alphaY = rdy2;
    const real alphaZ = rdz2;

    const real rbeta = 2. * (alphaX + alphaY + alphaZ);
    const real beta = 1. / rbeta;

    tmp.copy_data(out);
#pragma acc data present(out, tmp, b)
    {
        size_t it = 0;
        for (int i = 0; i < m_n_relax; i++){
            JacobiDiffuse::JacobiStep(level, out, tmp, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
            boundary->apply_boundary(data_out, level, type, sync);
            Field::swap(tmp, out);
            std::swap(data_tmp, data_out);
            it++;
        }

        if (it % 2 != 0) {  // get data from tmp field when number of iterations is odd
            out.copy_data(tmp);
        }
    }
}
void VCycleMG::call_solve_jacobi(Field &out, Field &tmp, Field const &b, const size_t level, bool sync) {
    auto domain = Domain::getInstance();

    // local variables and parameters for GPU
    const size_t Nx = domain->get_Nx(level);
    const size_t Ny = domain->get_Ny(level);

    const real dx = domain->get_dx(level);
    const real dy = domain->get_dy(level);
    const real dz = domain->get_dz(level);

    real *data_out = out.data;
    real *data_tmp = tmp.data;

    FieldType type = out.get_type();

    BoundaryController *boundary = BoundaryController::getInstance();

    size_t *data_inner_list = boundary->get_inner_list_level_joined();
    const size_t bsize_i = boundary->get_size_inner_list_level_joined();

    // start/end
    // inner
    const size_t start_i = boundary->get_inner_list_level_joined_start(level);
    const size_t end_i = boundary->get_inner_list_level_joined_end(level) + 1;

    boundary->apply_boundary(data_out, level, type, sync);

    const real rdx2 = 1. / (dx * dx);
    const real rdy2 = 1. / (dy * dy);
    const real rdz2 = 1. / (dz * dz);

    // preparation for diffusion step (alpha, beta)
    const real alphaX = rdx2;
    const real alphaY = rdy2;
    const real alphaZ = rdz2;

    const real rbeta = 2. * (alphaX + alphaY + alphaZ);
    const real beta = 1. / rbeta;

    tmp.copy_data(out);
#pragma acc data present(out, tmp, b)
    {
        size_t it = 0;
        real sum;
        real res = 10000.;

        const size_t neighbour_i = 1;
        const size_t neighbour_j = Nx;
        const size_t neighbour_k = Nx * Ny;
        while (res > m_diffusion_tol_res && it < m_diffusion_max_iter) {
            JacobiDiffuse::JacobiStep(level, out, tmp, b, alphaX, alphaY, alphaZ, beta, m_dsign, m_w, sync);
            boundary->apply_boundary(data_out, level, type, sync);

            sum = 0.;

#pragma acc parallel loop independent present(out, tmp, data_inner_list[:bsize_i]) async
            for (size_t j = start_i; j < end_i; ++j) {
                const size_t i = data_inner_list[j];
                res = b[i] - (rdx2 * (out[i - neighbour_i] - 2 * out[i] + data_out[i + neighbour_i])
                           +  rdy2 * (out[i - neighbour_j] - 2 * out[i] + data_out[i + neighbour_j])
                           +  rdz2 * (out[i - neighbour_k] - 2 * out[i] + data_out[i + neighbour_k])); //res = rbeta*(data_out[i] - data_tmp[i]);
                sum += res * res;
            }
            // info: in nvvp profile 8byte size copy from to device to/from pageable due to sum!

#pragma acc wait
            res = sqrt(sum);

            it++;

            Field::swap(tmp, out);
            std::swap(data_tmp, data_out);
        }

        if (it % 2 != 0) {  // get data from tmp field when number of iterations is odd
            out.copy_data(tmp);
        }
    }
}
