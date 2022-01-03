/// \file       Cube.cpp
/// \brief      
/// \date       Sep 29, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "Cube.h"

#include "../domain/DomainData.h"

void Cube::update_source(Field &out, real) {
    out.copy_data(m_source_field);
    if (m_has_noise) {
        out *= m_noise_maker->random_field(out.get_size());
    }
}

Cube::Cube(const Settings::solver::sources::cube &settings) :
        m_settings(settings),
        m_source_field(FieldType::T) {
    set_up();
}

void Cube::set_up() {
    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_number_of_cells(CoordinateAxis::X);
    size_t Ny = domain_data->get_number_of_cells(CoordinateAxis::Y);

    CoordinateAxis axis = CoordinateAxis::X;
    size_t index_start_x = Utility::get_index(m_settings.coords_start[axis],
                                              domain_data->get_spacing(axis),
                                              domain_data->get_start_coord_PD(axis));
    size_t index_end_x = Utility::get_index(m_settings.coords_end[axis],
                                            domain_data->get_spacing(axis),
                                            domain_data->get_end_coord_PD(axis));

    axis = CoordinateAxis::Y;
    size_t index_start_y = Utility::get_index(m_settings.coords_start[axis],
                                              domain_data->get_spacing(axis),
                                              domain_data->get_start_coord_PD(axis));
    size_t index_end_y = Utility::get_index(m_settings.coords_end[axis],
                                            domain_data->get_spacing(axis),
                                            domain_data->get_end_coord_PD(axis));

    axis = CoordinateAxis::Z;
    size_t index_start_z = Utility::get_index(m_settings.coords_start[axis],
                                              domain_data->get_spacing(axis),
                                              domain_data->get_start_coord_PD(axis));
    size_t index_end_z = Utility::get_index(m_settings.coords_end[axis],
                                              domain_data->get_spacing(axis),
                                              domain_data->get_end_coord_PD(axis));

    for (size_t i = index_start_x; i <= index_end_x; i++) {
        for (size_t j = index_start_y; j <= index_end_y; j++) {
            for (size_t k = index_start_z; k <= index_end_z; k++) {
                size_t index = IX(i, j, k, Nx, Ny);
                m_source_field[index] = m_settings.value;
            }
        }
    }
    m_source_field.update_dev();
}

void Cube::read_header_part(std::string &header) {

}

std::string Cube::write_header_part() {
    return std::string();
}