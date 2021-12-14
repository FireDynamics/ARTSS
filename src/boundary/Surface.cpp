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
        m_name(name),
        m_patch(patch),
        m_level(0) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
    m_logger->debug("################ SURFACE ################");
#endif
#ifndef BENCHMARKING
    m_logger->debug("read surface '{}' for patch {}", m_name, PatchObject::get_patch_name(m_patch));
#endif


    //TODO
    auto domain = DomainData::getInstance();
    size_t Nx = domain->get_Nx(m_level);
    size_t Ny = domain->get_Ny(m_level);
    createSurface(Nx, Ny);
    print();
#ifndef BENCHMARKING
    m_logger->info("----------------END SURFACE ----------------");
#endif
}

Surface::Surface(const std::string &name, Patch patch, Coordinate<size_t> &start, Coordinate<size_t> &end, size_t level) :
        m_patch(patch), m_name(name), m_start(start), m_end(end) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
    m_logger->info("################ SURFACE ################");
#endif

    DomainData *domain = DomainData::getInstance();
    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);

    init(Nx, Ny);
#ifndef BENCHMARKING
    m_logger->info("----------------END SURFACE ----------------");
#endif
}

Surface::~Surface() {
    for (BoundaryData *bd: dataList) {
        delete (bd);
    }
    delete (m_surface_list);
}

void Surface::print() {
#ifndef BENCHMARKING
    m_logger->info("Surface name {}", m_name);
    m_logger->info("strides: X: {}, Y: {}, Z:{}",
                   get_stride_z(),
                   get_stride_y(),
                   get_stride_x());
    m_logger->info("size of Surface: {}", m_size_surfaceList);
    m_logger->info("coords: ({}|{}) ({}|{}) ({}|{})",
                   m_start[X], m_end[X],
                   m_start[Y], m_end[Y],
                   m_start[Z], m_end[Z]);
#endif
}

void Surface::init(size_t Nx, size_t Ny) {
    m_size_surfaceList = get_stride_x() * get_stride_y() * get_stride_z();
    m_surface_list = new size_t[m_size_surfaceList];

    createSurface(Nx, Ny);
    print();
}

void Surface::createSurface(size_t Nx, size_t Ny) {
    size_t counter = 0;
#ifndef BENCHMARKING
    m_logger->info("list size of sList: {}", m_size_surfaceList);
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
                   m_size_surfaceList);
    m_logger->debug("end of creating sList");
#endif
}

void Surface::apply_boundary_conditions(Field &field, bool sync) {
    // m_bdc_boundary->apply_boundary_condition(dataField, indexFields, patch_starts, patch_ends, fieldType, level, sync);
    //TODO(issue 5)
}

size_t Surface::get_matching_index(real surface_coordinate, real spacing, real start_coordinate) {
    return static_cast<int>(round((-start_coordinate + surface_coordinate) / spacing));
}
