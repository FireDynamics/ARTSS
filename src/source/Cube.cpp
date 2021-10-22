/// \file       Cube.cpp
/// \brief      
/// \date       Sep 29, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "Cube.h"
#include "../DomainData.h"
#include "../boundary/BoundaryController.h"

void Cube::update_source(Field &out, real) {
    auto boundary = BoundaryController::getInstance();
    out.copy_data(m_source_field);
    if (m_has_noise) {
        out *= m_noise_maker->random_field(out.get_size());
    }
}

Cube::Cube(
        real value,
        real x_start, real y_start, real z_start,
        real x_end, real y_end, real z_end) :
    m_source_field(FieldType::T) {
    set_up(value, x_start, y_start, z_start, x_end, y_end, z_end);
}

void Cube::set_up(
        real value,
        real x_start, real y_start, real z_start,
        real x_end, real y_end, real z_end) {
    DomainData *domain = DomainData::getInstance();
    size_t Nx = domain->get_Nx();
    size_t Ny = domain->get_Ny();

    real X1 = domain->get_X1();
    real Y1 = domain->get_Y1();
    real Z1 = domain->get_Z1();

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    size_t index_start_x = Utility::get_index(x_start, dx, X1);
    size_t index_end_x = Utility::get_index(x_end, dx, X1);
    size_t index_start_y = Utility::get_index(y_start, dy, Y1);
    size_t index_end_y = Utility::get_index(y_end, dy, Y1);
    size_t index_start_z = Utility::get_index(z_start, dz, Z1);
    size_t index_end_z = Utility::get_index(z_end, dz, Z1);

    for (size_t i = index_start_x; i <= index_end_x; i++) {
        for (size_t j = index_start_y; j <= index_end_y; j++) {
            for (size_t k = index_start_z; k <= index_end_z; k++) {
                size_t index = IX(i, j, k, Nx, Ny);
                m_source_field[index] = value;
            }
        }
    }
    m_source_field.update_dev();
}
