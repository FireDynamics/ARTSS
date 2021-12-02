//
// Created by linh on 01.10.19.
//

#include "Surface.h"


// TODO(issue 15): surface implementing
//  - create file description
//  - remove surface cells from domain domain boundary cells
//  - develop a (better) concept for boundary conditions
//  - make sure that surfaces does not extend to corner or edge cells
//  - consider moving parsing to BoundaryController.cpp same as for obstacles,
//    may be not possible/practical depending on the concept for BC
Surface::Surface(Settings::Settings const &settings,
                 real x1, real x2,
                 real y1, real y2,
                 real z1, real z2,
                 const std::string &name) :
        m_settings(settings),
        m_name(name),
        m_level(0) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(m_settings, typeid(this).name());
    m_logger->debug("################ SURFACE ################");
#endif
    DomainData *domain_data = DomainData::getInstance();

    real dx = domain_data->get_dx();
    real dy = domain_data->get_dy();
    real dz = domain_data->get_dz();

    real X1 = domain_data->get_X1();
    real Y1 = domain_data->get_Y1();
    real Z1 = domain_data->get_Z1();

    size_t i1 = get_matching_index(x1, dx, X1) + 1;  // plus 1 for ghost cell
    size_t j1 = get_matching_index(y1, dy, Y1) + 1;
    size_t k1 = get_matching_index(z1, dz, Z1) + 1;
    m_start.set_coordinate(i1, j1, k1);

    size_t i2 = get_matching_index(x2, dx, X1);
    size_t j2 = get_matching_index(y2, dy, Y1);
    size_t k2 = get_matching_index(z2, dz, Z1);
    m_end.set_coordinate(i2, j2, k2);

    for (size_t axis = 0; axis < number_of_axis; axis++) {
        if (m_start[axis] == domain_data->get_start_index_CD(CoordinateAxis(axis)) +-1) {
            m_patch = Patch(axis * 2);
        }
        if (m_end[axis] == domain_data->get_end_index_CD(CoordinateAxis(axis)) + 1) {
            m_patch = Patch(axis * 2 + 1);
        }
    }
#ifndef BENCHMARKING
    m_logger->debug("read surface '{}' for patch {}", m_name, PatchObject::get_patch_name(m_patch));
#endif
    init();
#ifndef BENCHMARKING
    print();
    m_logger->info("----------------END SURFACE ----------------");
#endif
}

Surface::Surface(Settings::Settings const &settings,
                 Coordinate<size_t> &coords_start, Coordinate<size_t> &coords_end,
                 size_t level,
                 const std::string &name,
                 Patch patch) :
        m_settings(settings),
        m_patch(patch),
        m_name(name),
        m_start(coords_start), m_end(coords_end),
        m_level(level) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(settings, typeid(this).name());
    m_logger->info("################ SURFACE ################");
#endif

    init();
#ifndef BENCHMARKING
    print();
    m_logger->info("----------------END SURFACE ----------------");
#endif
}

Surface::~Surface() {
    delete (m_surface_list);
}

void Surface::print() {
#ifndef BENCHMARKING
    m_logger->info("Surface name {} on Patch {}", m_name, PatchObject::get_patch_name(m_patch));
    m_logger->info("strides: X: {}, Y: {}, Z:{}",
                   get_stride(Z),
                   get_stride(Y),
                   get_stride(X));
    m_logger->info("size of Surface: {}", m_size_surface_list);
    m_logger->info("coords: ({}|{}) ({}|{}) ({}|{})",
                   m_start[X], m_end[X],
                   m_start[Y], m_end[Y],
                   m_start[Z], m_end[Z]);
#endif
}

void Surface::init() {
    m_size_surface_list = get_stride(X) * get_stride(Y) * get_stride(Z);
    m_surface_list = new size_t[m_size_surface_list];

    DomainData *domain_data = DomainData::getInstance();
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
