/// \file         ColoredGaussSeidelDiffuse.cpp
/// \brief        solves Diffusion equation \f$ \partial_t \phi_2 = \nu \nabla^2 \phi_2 \f$
///               via calculated iterations of CGS step (dependent on residual/ maximal 1000)
/// \date         June 21, 2016
/// \author       My Linh Wuerzburger
/// \copyright    <2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.


#include "ColoredGaussSeidelDiffuse.h"

#include <cmath>

#include "../domain/DomainController.h"
#include "../domain/DomainData.h"


// ========================== Constructor =================================
ColoredGaussSeidelDiffuse::ColoredGaussSeidelDiffuse(Settings::Settings const &settings) :
        m_settings(settings),
        m_dsign(1),
        m_w(m_settings.get_real("solver/diffusion/w")),
        m_max_iter(m_settings.get_size_t("solver/diffusion/max_iter")),
        m_tol_res(m_settings.get_real("solver/diffusion/tol_res")) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
#endif
    create_red_black_lists(0, odd_indices, even_indices);
}


// ========================== Diffuse =================================
// ************************************************************************
/// \brief  solves Diffusion equation \f$ \partial_t \phi_2 = \nu \ nabla^2 \phi_2 \f$
///     via calculated iterations of CGS step (dependent on residual/ maximal iterations)
/// \param  out     output pointer
/// \param  in      input pointer
/// \param  b       source pointer
/// \param  D     diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  sync    synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::diffuse(
        Field &out, const Field &in, Field const &b,
        const real D, bool sync) {
    out.copy_data(in);  // cgs only calculates on out
    auto domain_data = DomainData::getInstance();
    auto domain_controller = DomainController::getInstance();

    auto size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();

#pragma acc data present(out, b, in)
    {
        const real dx = domain_data->get_spacing(CoordinateAxis::X);  // due to unnecessary parameter passing of *this
        const real dy = domain_data->get_spacing(CoordinateAxis::Y);
        const real dz = domain_data->get_spacing(CoordinateAxis::Z);

        const real reciprocal_dx = 1. / dx;
        const real reciprocal_dy = 1. / dy;
        const real reciprocal_dz = 1. / dz;

        const real dt = m_settings.get_real("physical_parameters/dt");
        const real alpha_x = D * dt * reciprocal_dx * reciprocal_dx;  // due to better pgi handling of scalars (instead of arrays)
        const real alpha_y = D * dt * reciprocal_dy * reciprocal_dy;
        const real alpha_z = D * dt * reciprocal_dz * reciprocal_dz;

        const real reciprocal_beta = (1. + 2. * (alpha_x + alpha_y + alpha_z));
        const real beta = 1. / reciprocal_beta;

        size_t it = 0;

        real sum;
        real res = 1.;
        while (res > m_tol_res && it < m_max_iter) {
            in.copy_data(out);  // necessary for calculation of residuum
            colored_gauss_seidel_step(out, b, alpha_x, alpha_y, alpha_z, beta, m_dsign, m_w, odd_indices, even_indices, sync);
            domain_controller->apply_boundary(out, sync);

            sum = 0;

#pragma acc parallel loop independent present(out, in, domain_inner_list[:size_domain_inner_list]) async
            for (size_t j = 0; j < size_domain_inner_list; ++j) {
                const size_t index = domain_inner_list[j];
                res = reciprocal_beta * (out[index] - in[index]);
                sum += res * res;
            }

#pragma acc wait
            res = sqrt(sum);
            it++;
        } //end while

        if (sync) {
#pragma acc wait
        }

#ifndef BENCHMARKING
        m_logger->info("Number of iterations: {}", it);
        m_logger->info("Colored Gauss-Seidel ||res|| = {:.5e}", res);
#endif
    } //end data region
}

// ===================== Turbulent version ============================
// ************************************************************************
/// \brief  solves Diffusion equation \f$ \partial_t \phi_2 = \nu_{eff} \ nabla^2 \phi_2 \f$
///         via calculated iterations of CGS step (dependent on residual/ maximal iterations)
/// \param  out         output pointer
/// \param  in          input pointer
/// \param  b           source pointer
/// \param  D           diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  EV          eddy viscosity (nu_turb)
/// \param  sync        synchronization boolean (true=sync (default), false=async)
// ************************************************************************
void ColoredGaussSeidelDiffuse::diffuse(
        Field &out, const Field &in, const Field &b,
        const real D, const Field &EV, bool sync) {
    out.copy_data(in);  // cgs only calculates on out
    auto domain_data = DomainData::getInstance();
    auto domain_controller = DomainController::getInstance();

    auto size_domain_inner_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();

#pragma acc data present(out, b, EV)
    {
        const size_t Nx = domain_data->get_Nx();
        const size_t Ny = domain_data->get_Ny();

        const real dx = domain_data->get_dx();
        const real dy = domain_data->get_dy();
        const real dz = domain_data->get_dz();

        const real reciprocal_dx = 1. / dx;
        const real reciprocal_dy = 1. / dy;
        const real reciprocal_dz = 1. / dz;

        real dt = m_settings.get_real("physical_parameters/dt");

        real alpha_x, alpha_y, alpha_z, reciprocal_beta;  // calculated in colored_gauss_seidel_step!

        size_t it = 0;
        real sum;
        real res = 1.;

        const size_t neighbour_i = 1;
        const size_t neighbour_j = Nx;
        const size_t neighbour_k = Nx * Ny;
        while (res > m_tol_res && it < m_max_iter) {
            in.copy_data(out);  // necessary for calculation of residuum
            colored_gauss_seidel_step(out, b, m_dsign, m_w, D, EV, dt, odd_indices, even_indices, sync);
            domain_controller->apply_boundary(out, sync);

            sum = 0;

#pragma acc parallel loop independent present(out, b, EV, domain_inner_list[:size_domain_inner_list]) async
            for (size_t j = 0; j < size_domain_inner_list; ++j) {
                const size_t i = domain_inner_list[j];
                alpha_x = (D + EV[i]) * dt * reciprocal_dx * reciprocal_dx;
                alpha_y = (D + EV[i]) * dt * reciprocal_dy * reciprocal_dy;
                alpha_z = (D + EV[i]) * dt * reciprocal_dz * reciprocal_dz;
                reciprocal_beta = (1. + 2. * (alpha_x + alpha_y + alpha_z));

                res = reciprocal_beta * (out[i] - in[i]);
                sum += res * res;
            }

#pragma acc wait
            res = sqrt(sum);
            it++;

        } //end while

        if (sync) {
#pragma acc wait
        }

#ifndef BENCHMARKING
        m_logger->info("Number of iterations: {}", it);
        m_logger->info("Colored Gauss-Seidel ||res|| = {.5e}", res);
#endif
    }
}

//========================== Iteration step =================================
// *****************************************************************************
/// \brief  applies single CGS step on red-black grid
/// \param  out      output pointer
/// \param  b        source pointer
/// \param  alpha    2-dimensional array;
///                      \f$ (reciprocal_dx^2, reciprocal_dy^2)\f$ for pressure,
///                      \f$ (\nu\cdot dt\cdot reciprocal_dx^2, \nu\cdot dt\cdot reciprocal_dy^2)\f$  for velocity,
///                      \f$ (\kappa\cdot dt\cdot reciprocal_dx^2, \kappa\cdot dt\cdot reciprocal_dy^2)\f$ for temperature
/// \param  beta     \f$ 1./(2\cdot(\alpha_0 + \alpha_1) + 0.5)\f$ for pressure
///                  \f$ 1/(2\cdot(\alpha_0 + \alpha_1) + 1)\f$ for velocity
/// \param  dsign    sign (\a -1. for pressure, \a 1. else)
/// \param  w        weight (1. - diffusion, 2./3. - multigrid)
/// \param  sync     synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(
        Field &out, const Field &b,
        const real alpha_x, const real alpha_y, const real alpha_z,
        const real beta, const real dsign, const real w,
        const std::vector<size_t> &odd,
        const std::vector<size_t> &even,
        bool) {

    auto domain_data = DomainData::getInstance();
    // local parameters for GPU
    const size_t Nx = domain_data->get_number_of_cells(CoordinateAxis::X);
    const size_t Ny = domain_data->get_number_of_cells(CoordinateAxis::Y);

    auto d_out = out.data;
    auto d_b = b.data;

    const size_t *data_odd = odd.data();
    size_t size_odd = odd.size();
    const size_t *data_even = even.data();
    size_t size_even = even.size();
    // red
#pragma acc parallel loop independent present(out, b, data_even[:size_even])
    for (size_t i = 0; i < size_even; i++) {
        colored_gauss_seidel_stencil(data_even[i], d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, Nx, Ny);
    }

#pragma acc wait
    // black
#pragma acc parallel loop independent present(out, b, data_odd[:size_odd])
    for (size_t i = 0; i < size_odd; i++) {
        colored_gauss_seidel_stencil(data_odd[i], d_out, d_b, alpha_x, alpha_y, alpha_z, dsign, beta, w, Nx, Ny);
    }
#pragma acc wait
}

// ============== Turbulent version of iteration step =================
// *****************************************************************************
/// \brief  applies single CGS step on red-black grid
/// \param  out      output pointer
/// \param  b        source pointer
/// \param  dsign    sign (\a -1. for pressure, \a 1. else)
/// \param  w        weight (1. - diffusion, 2./3. - multigrid)
/// \param  D        diffusion coefficient (nu - velocity, kappa - temperature)
/// \param  EV       eddy viscosity (nu_turb)
/// \param  dt       time step
/// \param  sync     synchronization boolean (true=sync (default), false=async)
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::colored_gauss_seidel_step(
        Field &out, const Field &b,
        const real dsign, const real w, const real D,
        const Field &EV,
        const real dt,
        const std::vector<size_t> &odd,
        const std::vector<size_t> &even,
        bool) {

    auto domain_data = DomainData::getInstance();
    // local parameters for GPU
    const size_t Nx = domain_data->get_Nx(out.get_level());
    const size_t Ny = domain_data->get_Ny(out.get_level());

    const real dx = domain_data->get_dx(out.get_level());  // due to unnecessary parameter passing of *this
    const real dy = domain_data->get_dy(out.get_level());
    const real dz = domain_data->get_dz(out.get_level());

    const real reciprocal_dx = 1. / dx;  // due to unnecessary parameter passing of *this
    const real reciprocal_dy = 1. / dy;
    const real reciprocal_dz = 1. / dz;

    real aX, aY, aZ, bb, rb;  // multipliers calculated

    auto d_out = out.data;
    auto d_b = b.data;

    const size_t *data_odd = odd.data();
    size_t size_odd = odd.size();
    const size_t *data_even = even.data();
    size_t size_even = even.size();
    // red
#pragma acc parallel loop independent present(out, b, data_even[:size_even], EV)
    for (size_t i = 0; i < size_even; i++) {
        size_t index = data_even[i];
        aX = (D + EV[index]) * dt * reciprocal_dx * reciprocal_dx;
        aY = (D + EV[index]) * dt * reciprocal_dy * reciprocal_dy;
        aZ = (D + EV[index]) * dt * reciprocal_dz * reciprocal_dz;

        rb = (1. + 2. * (aX + aY + aZ));
        bb = 1. / rb;

        colored_gauss_seidel_stencil(index, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
    }
    // black
#pragma acc parallel loop independent present(out, b, data_odd[:size_odd], EV)
    for (size_t i = 0; i < size_odd; i++) {
        size_t index = data_odd[i];
        aX = (D + EV[index]) * dt * reciprocal_dx * reciprocal_dx;
        aY = (D + EV[index]) * dt * reciprocal_dy * reciprocal_dy;
        aZ = (D + EV[index]) * dt * reciprocal_dz * reciprocal_dz;

        rb = (1. + 2. * (aX + aY + aZ));
        bb = 1. / rb;

        colored_gauss_seidel_stencil(index, d_out, d_b, aX, aY, aZ, dsign, bb, w, Nx, Ny);
    }
}

// ========================= CGS stencil ==============================
// ************************************************************************
/// \brief  single CGS step
/// \param  i        index in x-direction
/// \param  j        index in y-direction
/// \param  k        index in z-direction
/// \param  out      output pointer
/// \param  b        source pointer
/// \param  alpha    2-dimensional array;
///                      \f$ (reciprocal_dx^2, reciprocal_dy^2)\f$ for pressure,
///                      \f$ (\nu\cdot dt\cdot reciprocal_dx^2, \nu\cdot dt\cdot reciprocal_dy^2)\f$  for velocity,
///                      \f$ (\kappa\cdot dt\cdot reciprocal_dx^2, \kappa\cdot dt\cdot reciprocal_dy^2)\f$ for temperature
/// \param  beta     \f$ 1./(2\cdot(\alpha_0 + \alpha_1) + 0.5)\f$ for pressure
///                  \f$ 1/(2\cdot(\alpha_0 + \alpha_1) + 1)\f$ for velocity
/// \param  dsign    sign (\a -1. for pressure, \a 1. else)
/// \param  w        weight (1. - diffusion, 2./3. - multigrid)
/// \param  Nx       number of cells in x-direction of computational domain
/// \param  Ny       number of cells in y-direction
// ***************************************************************************************
void ColoredGaussSeidelDiffuse::colored_gauss_seidel_stencil(
        const size_t index,
        real *out, const real *b,
        const real alpha_x, const real alpha_y, const real alpha_z,
        const real dsign, const real beta, const real w,
        const size_t Nx, const size_t Ny) {
    size_t neighbour_i = 1;
    size_t neighbour_j = Nx;
    size_t neighbour_k = Nx * Ny;

    real d_out_x = *(out + index + neighbour_i); // per value (not access) necessary due to performance issues
    real d_out_x2 = *(out + index - neighbour_i);
    real d_out_y = *(out + index + neighbour_j);
    real d_out_y2 = *(out + index - neighbour_j);
    real d_out_z = *(out + index + neighbour_k);
    real d_out_z2 = *(out + index - neighbour_k);
    real d_b = *(b + index);
    real r_out = *(out + index);

    real out_h = beta * (dsign * d_b
                         + alpha_x * (d_out_x + d_out_x2)
                         + alpha_y * (d_out_y + d_out_y2)
                         + alpha_z * (d_out_z + d_out_z2));

    *(out + index) = (1 - w) * r_out + w * out_h;
}

void ColoredGaussSeidelDiffuse::create_red_black_lists(
        size_t level,
        std::vector<size_t> &odd_indices,
        std::vector<size_t> &even_indices) {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_number_of_cells(CoordinateAxis::X, level);
    size_t Ny = domain_data->get_number_of_cells(CoordinateAxis::Y, level);
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    size_t start = domain_controller->get_domain_inner_list_level_joined_start(level);
    size_t end = domain_controller->get_domain_inner_list_level_joined_end(level);
    for (size_t c = start; c <= end; c++) {
        size_t index = domain_inner_list[c];
        size_t k = getCoordinateK(index, Nx, Ny);
        size_t j = getCoordinateJ(index, Nx, Ny, k);
        size_t i = getCoordinateI(index, Nx, Ny, j, k);

        if ((i + j + k) % 2 == 0) {
            even_indices.emplace_back(index);
        } else {
            odd_indices.emplace_back(index);
        }
    }
    size_t *odd __attribute__((unused)) = odd_indices.data();
    size_t size_odd __attribute__((unused)) = odd_indices.size();
    size_t *even __attribute__((unused)) = even_indices.data();
    size_t size_even __attribute__((unused)) = even_indices.size();
#pragma acc enter data create(odd[:size_odd])
#pragma acc enter data create(even[:size_even])
}
