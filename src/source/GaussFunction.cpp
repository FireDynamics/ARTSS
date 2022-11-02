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
        m_field_spatial_values(FieldType::RHO, 0),
        m_settings(settings) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
    m_logger->debug("create gauss function with parameters dimension: {}, position: {}, hc: {}, HRR: {}, tau: {}",
                    settings.dimension, settings.position, settings.heat_capacity, settings.heat_release_rate, settings.tau);
#endif
    create_spatial_values();
}

void GaussFunction::update_source(Field &out, real t_cur) {
    auto time_val = get_time_value(t_cur);
    out.copy_data(m_field_spatial_values);
    out *= time_val;
    if (m_has_noise) {
        if (m_absolute) {
            out += m_noise_maker->random_field(out.get_size());
        } else {
            auto noise = m_noise_maker->random_field(out.get_size());
            noise += 1;
            out *= noise;
        }
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
    const Coordinate<size_t> number_of_cells = domain_data->get_number_of_cells();
    const Coordinate<real> start_coord_PD = domain_data->get_start_coord_PD();
    const Coordinate<real> spacing = domain_data->get_spacing();

    Coordinate<real> r_sigma = m_settings.dimension;
    r_sigma *= m_settings.dimension;
    r_sigma *= 2;
    for (size_t axis = 0; axis < number_of_axes; axis++) {
        r_sigma[axis] = 1. / r_sigma[axis];
    }

    // set Gaussian to cells
    auto domain_controller = DomainController::getInstance();
    size_t *domain_inner_list = domain_controller->get_domain_inner_list_level_joined();
    auto size_domain_list = domain_controller->get_size_domain_inner_list_level_joined(0);

    real HRRrV;

    Coordinate<size_t> start_coordinates = Utility::get_index(m_settings.position);

    real V = 0.;
    for (size_t l = 0; l < size_domain_list; ++l) {
        const size_t idx = domain_inner_list[l];
        auto index_components = Utility::get_coordinates(idx, number_of_cells);

        if (DomainController::getInstance()->is_blocked_by_obstacle(start_coordinates, index_components)){
            continue;
        }

        auto midpoints = Utility::get_physical_coords_midpoint(start_coord_PD, index_components, spacing);
        Coordinate<real> delta(midpoints[X] - m_settings.position[X], midpoints[Y] - m_settings.position[Y], midpoints[Z] - m_settings.position[Z]);
        delta *= delta;
        delta *= r_sigma;
        real expr = std::exp( -(delta[X] + delta[Y] + delta[Z]));
        V += expr * spacing[CoordinateAxis::X] * spacing[CoordinateAxis::Y] * spacing[CoordinateAxis::Z];
    }

    HRRrV = m_settings.heat_release_rate / V;       // in case of concentration Ys*HRR
    real rcp = 1. / m_settings.heat_capacity;    // to get [K/s] for energy equation (d_t T), rho:=1, otherwise *1/rho; in case of concentration 1/Hc to get kg/m^3s

    for (size_t l = 0; l < size_domain_list; ++l) {
        const size_t idx = domain_inner_list[l];
        Coordinate<size_t> index_components = Utility::get_coordinates(idx, number_of_cells);

        if (DomainController::getInstance()->is_blocked_by_obstacle(start_coordinates, index_components)){
            continue;
        }
        auto midpoints = Utility::get_physical_coords_midpoint(start_coord_PD, index_components, spacing);
        Coordinate<real> delta(midpoints[X] - m_settings.position[X], midpoints[Y] - m_settings.position[Y], midpoints[Z] - m_settings.position[Z]);
        delta *= delta;
        delta *= r_sigma;
        real expr = std::exp( -(delta[X] + delta[Y] + delta[Z]));
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
