/// \file       BuoyancyMMS.cpp
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "BuoyancyMMS.h"
#include "../boundary/BoundaryController.h"
#include "../utility/Parameters.h"
#include "../DomainData.h"

BuoyancyMMS::BuoyancyMMS() : m_source_field(FieldType::RHO) {
    set_up();
}

//===================================== Energy Source ====================================
// ***************************************************************************************
/// \brief  Manufactured (MMS) energy source in energy equation
/// \param  out   energy source
/// \param  t   time
// ***************************************************************************************
void BuoyancyMMS::set_up() {
    auto domain = DomainData::getInstance();
    // local variables and parameters for GPU
    size_t Nx = domain->get_Nx();
    size_t Ny = domain->get_Ny();

    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();

    real dx = domain->get_dx();
    real dy = domain->get_dy();

    auto params = Parameters::getInstance();

    real nu = params->get_real("physical_parameters/nu");
    real beta = params->get_real("physical_parameters/beta");
    real kappa = params->get_real("physical_parameters/kappa");
    real g = params->get_real("physical_parameters/g");
    real rhoa = params->get_real("initial_conditions/rhoa");
    real rbeta = 1. / beta;
    real rg = 1. / g;
    real c_nu = 2 * nu * M_PI * M_PI - 1;
    real c_kappa = 2 * kappa * M_PI * M_PI - 1;

    auto boundary = BoundaryController::getInstance();

    size_t *domain_list = boundary->get_domain_list_level_joined();
    auto size_domain_list = boundary->get_slice_size_domain_list_level_joined(0);

    // inner cells
    for (size_t l = 0; l < size_domain_list; ++l) {
        const size_t idx = domain_list[l];
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
