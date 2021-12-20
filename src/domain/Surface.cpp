/// \file       Surface.cpp
/// \brief      Data class of surface object
/// \date       Oct 01, 2019
/// \author     My Linh Wuerzburger
/// \copyright  <2015-2021> Forschungszentrum Juelich GmbH. All rights reserved.

#include "Surface.h"


Surface::Surface(real x1, real x2,
                 real y1, real y2,
                 real z1, real z2,
                 const std::string &name,
                 Patch patch) :
        m_patch(patch),
        m_name(name),
        m_level(0) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
    m_logger->debug("################ SURFACE ################");
#endif
    auto domain_data = DomainData::getInstance();
    Coordinate<real> tmp_start(x1, y1, z1);
    Coordinate<real> tmp_end(x2, y2, z2);
    for (size_t axis = 0; axis < number_of_axes; axis++) {
        auto coord_axis = CoordinateAxis(axis);
        m_start[axis] = get_matching_index(tmp_start[axis],
                                           domain_data->get_spacing(coord_axis, m_level),
                                           domain_data->get_start_coord_PD(coord_axis)) + 1;
        m_end[axis] = get_matching_index(tmp_end[axis],
                                         domain_data->get_spacing(coord_axis, m_level),
                                         domain_data->get_start_coord_PD(coord_axis));
    }

#ifndef BENCHMARKING
    m_logger->debug("read surface '{}' for patch {}", m_name, Mapping::get_patch_name(m_patch));
#endif
    init();
#ifndef BENCHMARKING
    print();
    m_logger->info("----------------END SURFACE ----------------");
#endif
}

Surface::Surface(Coordinate<size_t> &coords_start, Coordinate<size_t> &coords_end,
                 size_t level,
                 const std::string &name,
                 Patch patch) :
        m_patch(patch),
        m_name(name),
        m_start(coords_start), m_end(coords_end),
        m_level(level) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
    m_logger->info("################ SURFACE ################");
#endif

    init();
#ifndef BENCHMARKING
    print();
    m_logger->info("----------------END SURFACE ----------------");
#endif
}

void Surface::print() {
#ifndef BENCHMARKING
    m_logger->info("Surface name {} on Patch {}", m_name, Mapping::get_patch_name(m_patch));
    m_logger->info("strides: X: {}, Y: {}, Z:{}",
                   get_stride(X),
                   get_stride(Y),
                   get_stride(Z));
    m_logger->info("size of Surface: {}", m_size_surface_list);
    m_logger->info("coords: ({}|{}) ({}|{}) ({}|{})",
                   m_start[X], m_end[X],
                   m_start[Y], m_end[Y],
                   m_start[Z], m_end[Z]);
#endif
}

void Surface::init() {
    for (size_t axis = 0; axis < number_of_axes; axis++) {
        if (m_end[axis] != m_start[axis]) {
            m_strides[axis] = m_end[axis] - m_start[axis] + 1;
        }
    }
    m_size_surface_list = get_stride(X) * get_stride(Y) * get_stride(Z);
    m_surface_list.resize(m_size_surface_list);

    auto domain_data = DomainData::getInstance();
    size_t Nx = domain_data->get_number_of_inner_cells(CoordinateAxis::X);
    size_t Ny = domain_data->get_number_of_inner_cells(CoordinateAxis::Y);

    size_t counter = 0;
#ifndef BENCHMARKING
    m_logger->info("list size of sList: {}", m_size_surface_list);
#endif
    // fill sList with corresponding indices
    for (size_t k = m_start[Z]; k < m_end[Z]; ++k) {
        for (size_t j = m_start[Y]; j < m_end[Y]; ++j) {
            for (size_t i = m_start[X]; i < m_end[X]; ++i) {
                size_t idx = IX(i, j, k, Nx, Ny);
                m_surface_list[counter++] = idx;
            }
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("control create Surface ({}|{})", counter,
                    m_size_surface_list);
    m_logger->debug("end of creating sList");
#endif
}

size_t Surface::get_matching_index(real surface_coordinate, real spacing, real start_coordinate) {
    return static_cast<int>(round((-start_coordinate + surface_coordinate) / spacing));
}
