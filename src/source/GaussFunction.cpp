/// \file       GaussFunction.cpp
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "GaussFunction.h"
#include "../utility/Parameters.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"

GaussFunction::GaussFunction(
        real HRR, real cp,
        real x0, real y0, real z0,
        real sigma_x, real sigma_y, real sigma_z,
        real tau) : m_tau(tau) {
    m_field_spatial_values = new Field(FieldType::RHO, 0.);
    create_spatial_values(HRR, cp, x0, y0, z0, sigma_x, sigma_y, sigma_z);
}

GaussFunction::~GaussFunction() {
    auto data_spatial = m_field_spatial_values->data;
    size_t size = Domain::getInstance()->get_size();
#pragma acc exit data delete(data_spatial[:size])
    delete m_field_spatial_values;
}

void GaussFunction::update_source(Field *out, real t_cur) {
    size_t size = Domain::getInstance()->get_size();
    auto data_out = out->data;
    auto data_spatial = m_field_spatial_values->data;
    auto time_val = get_time_value(t_cur);

#pragma acc data present(data_out[:size], data_spatial[:size])
    {
#pragma acc parallel loop independent present(data_out[:size], data_spatial[:size]) async
        for (size_t i = 0; i < size; i++) {
            data_out[i] = data_spatial[i] * time_val;  // TODO * random value (5%)
        }
#pragma acc wait
    }
}

// ***************************************************************************************
/// \brief  Volumetric Gaussian temperature source in energy equation
/// \param  out   energy source
/// \param  HRR   total heat release rate
/// \param  cp    heat capacity
/// \param  x0    center of Gaussian (x-direction)
/// \param  y0    center of Gaussian (y-direction)
/// \param  z0    center of Gaussian (z-direction)
/// \param  sigma Radius of Gaussian
// ***************************************************************************************
void GaussFunction::create_spatial_values(real HRR, real cp,
        real x0, real y0, real z0,
        real sigma_x, real sigma_y, real sigma_z) {
    auto domain = Domain::getInstance();
    auto bsize = domain->get_size();
    // local variables and parameters for GPU
    auto d_out = m_field_spatial_values->data;
    auto level = m_field_spatial_values->get_level();

    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);

    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    real dx = domain->get_dx(level);
    real dy = domain->get_dy(level);
    real dz = domain->get_dz(level);

    // get parameters for Gaussian
    real sigma_x_2 = 2 * sigma_x * sigma_x;
    real r_sigma_x_2 = 1. / sigma_x_2;
    real sigma_y_2 = 2 * sigma_y * sigma_y;
    real r_sigma_y_2 = 1. / sigma_y_2;
    real sigma_z_2 = 2 * sigma_z * sigma_z;
    real r_sigma_z_2 = 1. / sigma_z_2;

    // set Gaussian to cells
    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_inner_list_level_joined();

    const auto multigrid = boundary->get_multigrid();
    const auto obst_size = multigrid.get_size_obstacle_list();
    auto bsize_i = boundary->get_size_inner_list();

    auto i0 = (x0 - X1) / dx;
    auto j0 = (y0 - Y1) / dy;
    auto k0 = (z0 - Z1) / dz;

    real V = 0.;
    for (size_t l = 0; l < bsize_i; ++l) {
        const size_t idx = d_iList[l];
        int k = getCoordinateK(idx, Nx, Ny);
        int j = getCoordinateJ(idx, Nx, Ny, k);
        int i = getCoordinateI(idx, Nx, Ny, j, k);

        real di = (i0 - i);
        real dj = (j0 - j);
        real dk = (k0 - k);

        for (auto obst_id=0; obst_id < obst_size; ++obst_id) {
            auto obst = multigrid.get_obstacle(level, obst_id);
            bool blocked = obst->line_crosses(i0, j0, k0, i, j, k);

            if (blocked) {
                d_out[idx] = 0.0;
                continue;
            }
        }

        real x_i = dx * di;
        real y_j = dy * dj;
        real z_k = dz * dk;

        real expr = std::exp(-(r_sigma_x_2 * (x_i * x_i)
                                + r_sigma_y_2 * (y_j * y_j)
                                + r_sigma_z_2 * (z_k * z_k)));
        d_out[idx] = expr;

        V += expr * dx * dy * dz;
    }

    real HRRrV = HRR / V;  // in case of concentration Ys*HRR
    real rcp = 1. / cp;    // to get [K/s] for energy equation (d_t T), rho:=1, otherwise *1/rho; in case of concentration 1/Hc to get kg/m^3s

    for (size_t l = 0; l < bsize_i; ++l) {
        const size_t idx = d_iList[l];
        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);

        auto x_i = (xi(i, X1, dx) - x0);
        auto y_j = (yj(j, Y1, dy) - y0);
        auto z_k = (zk(k, Z1, dz) - z0);
        real expr = std::exp(-(r_sigma_x_2 * x_i * x_i + r_sigma_y_2 * y_j * y_j + r_sigma_z_2 * z_k * z_k));
        real tmp = HRRrV * rcp * expr;
        if (tmp < 0.01) {
            tmp = 0;
        }
        d_out[idx] = tmp;
    }

#pragma acc enter data copyin(d_out[:bsize])
}

// ============================= Ramp up function for HRR source =========================
// ***************************************************************************************
/// \brief  Ramp up function (in time) for Gaussian source in energy equation
/// \param  t time
// ***************************************************************************************
real GaussFunction::get_time_value(real t_cur) {
    return tanh(t_cur / m_tau);
}
