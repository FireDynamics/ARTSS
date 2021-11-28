/// \file       BuoyancyMMS.cpp
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "BuoyancyMMS.h"
#include "../boundary/BoundaryController.h"
#include "../Domain.h"

BuoyancyMMS::BuoyancyMMS(Settings const &settings) :
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
    auto domain = Domain::getInstance();
    // local variables and parameters for GPU
    size_t Nx = domain->get_Nx();
    size_t Ny = domain->get_Ny();

    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();

    real dx = domain->get_dx();
    real dy = domain->get_dy();

    real nu = m_settings.get_real("physical_parameters/nu");
    real beta = m_settings.get_real("physical_parameters/beta");
    real kappa = m_settings.get_real("physical_parameters/kappa");
    real g = m_settings.get_real("physical_parameters/g");
    real rhoa = m_settings.get_real("initial_conditions/rhoa");
    real rbeta = 1. / beta;
    real rg = 1. / g;
    real c_nu = 2 * nu * M_PI * M_PI - 1;
    real c_kappa = 2 * kappa * M_PI * M_PI - 1;

    auto boundary = BoundaryController::getInstance();

    size_t *d_inner_list = boundary->get_inner_list_level_joined();
    size_t *d_bList = boundary->get_boundary_list_level_joined();

    auto bsize_i = boundary->get_size_inner_list();
    auto bsize_b = boundary->get_size_boundary_list();

    // inner cells
    for (size_t l = 0; l < bsize_i; ++l) {
        const size_t idx = d_inner_list[l];
        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);
        m_source_field[idx] = rhoa * rbeta * rg * 2 * c_nu * c_kappa * std::sin(M_PI * (xi(i, X1, dx) + yj(j, Y1, dy)));
    }

    // boundary cells
    for (size_t l = 0; l < bsize_b; ++l) {
        const size_t idx = d_bList[l];
        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);
        m_source_field[idx] = rhoa * rbeta * rg * 2 * c_nu * c_kappa * std::sin(M_PI * (xi(i, X1, dx) + yj(j, Y1, dy)));
    }

    m_source_field.update_dev();
}

void BuoyancyMMS::update_source(Field &out, real t_cur) {
    out.copy_data(m_source_field);
    out *= exp(-t_cur);
}
