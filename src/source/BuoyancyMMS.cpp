/// \file       BuoyancyMMS.cpp
/// \brief      
/// \date       Aug 11, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "BuoyancyMMS.h"
#include "../boundary/BoundaryController.h"
#include "../utility/Parameters.h"
#include "../Domain.h"

BuoyancyMMS::BuoyancyMMS() {
    m_source_field = new Field(FieldType::RHO, 0);
    set_up();
}

BuoyancyMMS::~BuoyancyMMS() {
    auto data_source = m_source_field->data;
    size_t size = Domain::getInstance()->get_size();
#pragma acc exit data delete(data_source[:size])
    delete m_source_field;
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
    auto bsize = domain->get_size();
    auto d_out = m_source_field->data;
    auto level = m_source_field->get_level();

    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);

    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();

    real dx = domain->get_dx(level);
    real dy = domain->get_dy(level);

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

    size_t *d_iList = boundary->get_inner_list_level_joined();
    size_t *d_bList = boundary->get_boundary_list_level_joined();

    auto bsize_i = boundary->get_size_inner_list();
    auto bsize_b = boundary->get_size_boundary_list();

    // inner cells
    for (size_t l = 0; l < bsize_i; ++l) {
        const size_t idx = d_iList[l];
        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);
        d_out[idx] = rhoa * rbeta * rg * 2 * c_nu * c_kappa * std::sin(M_PI * (xi(i, X1, dx) + yj(j, Y1, dy)));
    }

    // boundary cells
    for (size_t l = 0; l < bsize_b; ++l) {
        const size_t idx = d_bList[l];
        size_t k = getCoordinateK(idx, Nx, Ny);
        size_t j = getCoordinateJ(idx, Nx, Ny, k);
        size_t i = getCoordinateI(idx, Nx, Ny, j, k);
        d_out[idx] = rhoa * rbeta * rg * 2 * c_nu * c_kappa * std::sin(M_PI * (xi(i, X1, dx) + yj(j, Y1, dy)));
    }

#pragma acc enter data copyin(d_out[:bsize])
}

void BuoyancyMMS::update_source(Field *out, real t_cur) {
    auto boundary = BoundaryController::getInstance();
    size_t size = Domain::getInstance()->get_size();

    auto d_out = out->data;
    auto d_source = m_source_field->data;

#pragma acc data present(d_out[:size], d_source[:size])
    {
        size_t *d_iList = boundary->get_inner_list_level_joined();
        size_t *d_bList = boundary->get_boundary_list_level_joined();

        auto bsize_i = boundary->get_size_inner_list();
        auto bsize_b = boundary->get_size_boundary_list();

#pragma acc parallel loop independent present(d_out[:size], d_source[:size]) async
        // inner cells
        for (size_t l = 0; l < bsize_i; ++l) {
            const size_t idx = d_iList[l];
            d_out[idx] = d_source[idx] * exp(-t_cur);
        }

        // boundary cells
        for (size_t l = 0; l < bsize_b; ++l) {
            const size_t idx = d_bList[l];
            d_out[idx] = d_source[idx] * exp(-t_cur);
        }
    }
}

void BuoyancyMMS::read_header_part(std::string &header) {

}

std::string BuoyancyMMS::write_header_part() {
    return std::string();
}