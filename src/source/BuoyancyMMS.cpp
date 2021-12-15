/// \file       BuoyancyMMS.cpp
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "BuoyancyMMS.h"
#include "../boundary/DomainController.h"
#include "../boundary/DomainData.h"

BuoyancyMMS::BuoyancyMMS(Settings::Settings const &settings) :
        m_settings(settings),
        m_source_field(FieldType::RHO) {
    set_up();
}

//===================================== Energy Source ====================================
// ***************************************************************************************
/// \brief  Manufactured (MMS) energy source in energy equation
/// \param  out   energy source
/// \param  t   time
// ***************************************************************************************
void BuoyancyMMS::set_up() {
    auto domain_data = DomainData::getInstance();
    // local variables and parameters for GPU
    size_t Nx = domain_data->get_Nx();
    size_t Ny = domain_data->get_Ny();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();

    real nu = m_settings.get_real("physical_parameters/nu");
    real beta = m_settings.get_real("physical_parameters/beta");
    real kappa = m_settings.get_real("physical_parameters/kappa");
    real g = m_settings.get_real("physical_parameters/g");
    real rhoa = m_settings.get_real("initial_conditions/rhoa");
    real rbeta = 1. / beta;
    real rg = 1. / g;
    real c_nu = 2 * nu * M_PI * M_PI - 1;
    real c_kappa = 2 * kappa * M_PI * M_PI - 1;

    auto domain_controller = DomainController::getInstance();

    size_t *domain_list = domain_controller->get_domain_list_level_joined();
    size_t size_domain_list = domain_controller->get_slice_size_domain_list_level_joined(0);

    // inner cells
#pragma acc parallel loop independent present(m_source_field, domain_list[:size_domain_list]) async
    for (size_t l = 0; l < size_domain_list; ++l) {
        const size_t idx = domain_list[l];
        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);
        m_source_field[idx] = rhoa * rbeta * rg * 2 * c_nu * c_kappa * std::sin(M_PI * (xi(i, X1, dx) + yj(j, Y1, dy)));
    }
}

void BuoyancyMMS::update_source(Field &out, real t_cur) {
    out.copy_data(m_source_field);
    out *= exp(-t_cur);
}
