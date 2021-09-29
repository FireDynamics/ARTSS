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
        real sigma_x, real sigma_y, real sigma_z, real tau) :
    m_field_spatial_values(FieldType::RHO),
    m_tau(tau),
    m_HRR(HRR),
    m_cp(cp),
    m_x0(x0),
    m_y0(y0),
    m_z0(z0),
    m_sigma_x(sigma_x),
    m_sigma_y(sigma_y),
    m_sigma_z(sigma_z) {
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
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    size_t level = m_field_spatial_values.get_level();

    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);

    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    real dx = domain->get_dx(level);
    real dy = domain->get_dy(level);
    real dz = domain->get_dz(level);

    real sigma_x_2 = 2 * m_sigma_x * m_sigma_x;
    real r_sigma_x_2 = 1. / sigma_x_2;
    real sigma_y_2 = 2 * m_sigma_y * m_sigma_y;
    real r_sigma_y_2 = 1. / sigma_y_2;
    real sigma_z_2 = 2 * m_sigma_z * m_sigma_z;
    real r_sigma_z_2 = 1. / sigma_z_2;

    // set Gaussian to cells
    auto boundary = BoundaryController::getInstance();
    size_t *d_iList = boundary->get_inner_list_level_joined();

    auto bsize_i = boundary->get_size_inner_list();

    real HRRrV;

    real V = 0.;
    for (size_t l = 0; l < bsize_i; ++l) {
        const size_t idx = d_iList[l];
        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);

        auto x_i = xi(i, X1, dx) - m_x0;
        auto y_j = yj(j, Y1, dy) - m_y0;
        auto z_k = zk(k, Z1, dz) - m_z0;
        real expr = std::exp(-(r_sigma_x_2 * (x_i * x_i) + r_sigma_y_2 * (y_j * y_j) + r_sigma_z_2 * (z_k * z_k)));
        V += expr * dx * dy * dz;
    }

    HRRrV = m_HRR / V;       // in case of concentration Ys*HRR
    real rcp = 1. / m_cp;    // to get [K/s] for energy equation (d_t T), rho:=1, otherwise *1/rho; in case of concentration 1/Hc to get kg/m^3s

    for (size_t l = 0; l < bsize_i; ++l) {
        const size_t idx = d_iList[l];
        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);

        auto x_i = (xi(i, X1, dx) - m_x0);
        auto y_j = (yj(j, Y1, dy) - m_y0);
        auto z_k = (zk(k, Z1, dz) - m_z0);
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
    return tanh(t_cur / m_tau);
}

std::string GaussFunction::write_header_part() {
#ifdef ASSIMILATION
    std::string header_part = fmt::format("###coordinates;{};{};{}", m_x0, m_y0, m_z0);
    header_part.append(fmt::format("###volume;{};{};{}", m_sigma_x, m_sigma_y, m_sigma_z));
    header_part.append(fmt::format("###HRR;{};cp;{}", m_HRR, m_cp));
    return header_part;
#else
    return {};
#endif
}

void GaussFunction::read_header_part(std::string &header) {
#ifdef ASSIMILATION
    std::vector<std::string> lines = Utility::split(header, '\n');

    std::vector<std::string> coordinates = Utility::split(lines[0], ';');
    m_x0 = std::stod(coordinates[1]);
    m_y0 = std::stod(coordinates[2]);
    m_z0 = std::stod(coordinates[3]);

    std::vector<std::string> volume = Utility::split(lines[1], ';');
    m_sigma_x = std::stod(volume[1]);
    m_sigma_y = std::stod(volume[2]);
    m_sigma_z = std::stod(volume[3]);

    std::vector<std::string> heat_source = Utility::split(lines[2], ';');
    m_HRR = std::stod(heat_source[1]);
    m_cp = std::stod(heat_source[3]);

    create_spatial_values();
#endif
}