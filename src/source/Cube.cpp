/// \file       Cube.cpp
/// \brief      
/// \date       Sep 29, 2020
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich All rights reserved.
//
#include "Cube.h"
#include "../Domain.h"
#include "../boundary/BoundaryController.h"

void Cube::update_source(Field *out, real t_cur) {
    size_t size = Domain::getInstance()->get_size();
    auto data_out = out->data;
    auto data_spatial = m_source_field->data;
    std::uniform_int_distribution<int> dist(-m_steps, m_steps);

#pragma acc data present(data_out[:size], data_spatial[:size])
    {
#pragma acc parallel loop independent present(data_out[:size], data_spatial[:size]) async
        for (size_t i = 0; i < size; i++) {
            double no = dist(m_mt) * m_step_size;
            // ((1 + no) - no * (1 - m_has_noise)) = 1 + no - no + no * m_has_noise
            double noise = 1 + no * m_has_noise;  // * 1 if no noise, * (1-no) else
            data_out[i] = data_spatial[i] * noise;
        }
#pragma acc wait
    }
}

Cube::Cube(real value, real x_start, real y_start, real z_start, real x_end, real y_end, real z_end) {
    m_source_field = new Field(FieldType::T);
    set_up(value, x_start, y_start, z_start, x_end, y_end, z_end);
}

Cube::~Cube() {
    real *data = m_source_field->data;
    size_t size = Domain::getInstance()->get_size();
#pragma acc exit data delete(data[:size])
    delete m_source_field;
}

void Cube::set_up(real value, real x_start, real y_start, real z_start, real x_end, real y_end, real z_end) {
    Domain *domain = Domain::getInstance();
    real *data = m_source_field->data;
    size_t size = domain->get_size();

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
                data[index] = value;
            }
        }
    }
#pragma acc enter data copyin(data[:size])
}

void Cube::set_noise(real range, int seed, real step_size) {
    m_has_noise = true;
    m_steps = range / step_size;
    m_step_size = step_size;

    if (seed > 0) {
        m_mt = std::mt19937(seed);
    } else {
        std::random_device rd;
        m_mt = std::mt19937(rd());
    }
}
