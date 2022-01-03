/// \file       GaussFunction.cpp
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "GaussFunction.h"

#ifdef _OPENACC
#include "accel.h"
#endif

#include "../domain/DomainData.h"
#include "../domain/DomainController.h"

GaussFunction::GaussFunction(const Settings::solver::sources::gauss &settings) :
        m_field_spatial_values(FieldType::RHO),
        m_settings(settings) {
    create_spatial_values();
}

void GaussFunction::update_source(Field &out, real t_cur) {
    auto time_val = get_time_value(t_cur);
    out.copy_data(m_field_spatial_values);
    out *= time_val;

    if (m_has_noise) {
        out *= m_noise_maker->random_field(out.get_size());
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
void GaussFunction::create_spatial_values() {
    auto domain_data = DomainData::getInstance();
    // local variables and parameters for GPU
    size_t level = m_field_spatial_values.get_level();

    size_t Nx = domain_data->get_Nx(level);
    size_t Ny = domain_data->get_Ny(level);

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    real dx = domain_data->get_dx(level);
    real dy = domain_data->get_dy(level);
    real dz = domain_data->get_dz(level);

    real sigma_x_2 = 2 * m_settings.dimension[CoordinateAxis::X] * m_settings.dimension[CoordinateAxis::X];
    real r_sigma_x_2 = 1. / sigma_x_2;
    real sigma_y_2 = 2 * m_settings.dimension[CoordinateAxis::Y] * m_settings.dimension[CoordinateAxis::Y];
    real r_sigma_y_2 = 1. / sigma_y_2;
    real sigma_z_2 = 2 * m_settings.dimension[CoordinateAxis::Z] * m_settings.dimension[CoordinateAxis::Z];
    real r_sigma_z_2 = 1. / sigma_z_2;

    // set Gaussian to cells
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();

    auto size_domain_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    real HRRrV;

    real V = 0.;
    for (size_t l = 0; l < size_domain_list; ++l) {
        const size_t idx = domain_inner_list[l];
        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);

        auto x_i = xi(i, X1, dx) - m_settings.position[CoordinateAxis::X];
        auto y_j = yj(j, Y1, dy) - m_settings.position[CoordinateAxis::Y];
        auto z_k = zk(k, Z1, dz) - m_settings.position[CoordinateAxis::Z];
        real expr = std::exp(-(r_sigma_x_2 * (x_i * x_i) + r_sigma_y_2 * (y_j * y_j) + r_sigma_z_2 * (z_k * z_k)));
        V += expr * dx * dy * dz;
    }

    HRRrV = m_settings.heat_release_rate / V;       // in case of concentration Ys*HRR
    real rcp = 1. / m_settings.heat_capacity;    // to get [K/s] for energy equation (d_t T), rho:=1, otherwise *1/rho; in case of concentration 1/Hc to get kg/m^3s

    for (size_t l = 0; l < size_domain_list; ++l) {
        const size_t idx = domain_inner_list[l];
        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);

        auto x_i = (xi(i, X1, dx) - m_settings.position[CoordinateAxis::X]);
        auto y_j = (yj(j, Y1, dy) - m_settings.position[CoordinateAxis::Y]);
        auto z_k = (zk(k, Z1, dz) - m_settings.position[CoordinateAxis::Z]);
        real expr = std::exp(-(r_sigma_x_2 * x_i * x_i + r_sigma_y_2 * y_j * y_j + r_sigma_z_2 * z_k * z_k));
        m_field_spatial_values[idx] = HRRrV * rcp * expr;
    }
    m_field_spatial_values.update_dev();
}

// ============================= Ramp up function for HRR source =========================
// ***************************************************************************************
/// \brief  Ramp up function (in time) for Gaussian source in energy equation
/// \param  t time
// ***************************************************************************************
real GaussFunction::get_time_value(real t_cur) {
    return tanh(t_cur / m_settings.tau);
}
