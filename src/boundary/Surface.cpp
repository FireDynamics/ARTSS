//
// Created by linh on 01.10.19.
//
// TODO(linh): fix file header

#include "Surface.h"


// TODO(linh): duplicates ?
Surface::Surface(tinyxml2::XMLElement* element) {
    m_boundaryDataController = new BoundaryDataController();
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
    m_logger->info("################ SURFACE ################");
#endif
    auto domain = Domain::getInstance();

    real dx = domain->get_dx();
    real dy = domain->get_dy();
    real dz = domain->get_dz();

    real rdx = 1./dx;
    real rdy = 1./dy;
    real rdz = 1./dz;

    real X1 = domain->get_X1();
    real Y1 = domain->get_X1();
    real Z1 = domain->get_X1();

    size_t Nx = domain->get_Nx();
    size_t Ny = domain->get_Ny();

    m_surfaceID = element->IntAttribute("ID");

#ifndef BENCHMARKING
    m_logger->info("surface ID {}", m_surfaceID);
#endif

    real sx1 = element->DoubleAttribute("sx1");
    real sx2 = element->DoubleAttribute("sx2");
    real sy1 = element->DoubleAttribute("sy1");
    real sy2 = element->DoubleAttribute("sy2");
    real sz1 = element->DoubleAttribute("sz1");
    real sz2 = element->DoubleAttribute("sz2");
    real lsx = fabs(sx2 - sx1);
    real lsy = fabs(sy2 - sy1);
    real lsz = fabs(sz2 - sz1);
    m_strideX = static_cast<size_t> (floor(lsx * rdx + 1 + 0.5));
    m_strideY = static_cast<size_t> (floor(lsy * rdy + 1 + 0.5));
    m_strideZ = static_cast<size_t> (floor(lsz * rdz + 1 + 0.5));
    // round due to cells at boundary (sx as midpoint of cells in xml)
    m_i1 = static_cast<size_t> (round(fabs(sx1 - (X1-0.5*dx)) * rdx));
    m_j1 = static_cast<size_t> (round(fabs(sy1 - (Y1-0.5*dy)) * rdy));
    m_k1 = static_cast<size_t> (round(fabs(sz1 - (Z1-0.5*dz)) * rdz));

    createSurface(Nx, Ny);
    print();
#ifndef BENCHMARKING
    m_logger->info("----------------END SURFACE ----------------");
#endif
}

Surface::Surface(size_t surfaceID, size_t startIndex, size_t strideX, size_t strideY, size_t strideZ, size_t level) {
#ifndef BENCHMARKING
    m_logger = Utility::create_logger(typeid(this).name());
    m_logger->info("################ SURFACE ################");
#endif
    m_surfaceID = surfaceID;

    m_strideX = strideX;
    m_strideY = strideY;
    m_strideZ = strideZ;

    Domain* domain = Domain::getInstance();
    size_t Nx = domain->get_Nx(level);
    size_t Ny = domain->get_Ny(level);
    std::vector<size_t> coords = Utility::coordinateFromLinearIndex(startIndex, Nx, Ny);
    m_i1 = coords[0];
    m_j1 = coords[1];
    m_k1 = coords[2];

    init(Nx, Ny);
#ifndef BENCHMARKING
    m_logger->info("----------------END SURFACE ----------------");
#endif
}

Surface::~Surface() {
    for (BoundaryData* bd : dataList) {
        delete(bd);
    }
    delete(m_surfaceList);
}

void Surface::print() {
#ifdef BENCHMARKING
    return;
#else
    size_t i2 = get_i2();
    size_t j2 = get_j2();
    size_t k2 = get_k2();
    m_logger->info("Surface ID {}", m_surfaceID);
    m_logger->info("strides: X: {}, Y: {}, Z:{}", m_strideZ,
                                                  m_strideY,
                                                  m_strideX);
    m_logger->info("size of Surface: {}", m_size_surfaceList);
    m_logger->info("coords: ({}|{}) ({}|{}) ({}|{})", m_i1, i2,
                                                      m_j1, j2,
                                                      m_k1, k2);
#endif
}

void Surface::init(size_t Nx, size_t Ny) {
    m_size_surfaceList = m_strideX * m_strideY * m_strideZ;
    m_surfaceList = new size_t[m_size_surfaceList];

    createSurface(Nx, Ny);
    print();
}

void Surface::createSurface(size_t Nx, size_t Ny) {
    size_t counter = 0;
#ifndef BENCHMARKING
    m_logger->info("list size of sList: {}", m_size_surfaceList);
#endif

    // fill sList with corresponding indices
    for (size_t k = m_k1; k < m_k1 + m_strideZ; ++k) {
        for (size_t j = m_j1; j < m_j1 + m_strideY; ++j) {
            for (size_t i = m_i1; i < m_i1 + m_strideX; ++i) {
                size_t idx = (size_t)(IX(i,j,k,Nx,Ny));
                *(m_surfaceList + counter) = idx ;
                counter++;
            }
        }
    }
#ifndef BENCHMARKING
    m_logger->info("control create Surface ({}|{})", counter,
                                                     m_size_surfaceList);
    m_logger->info("end of creating sList");
#endif
}

void Surface::setBoundaryConditions(tinyxml2::XMLElement *xmlElement) {
    m_boundaryDataController->addBoundaryData(xmlElement);
}

size_t Surface::get_i2() {
    return m_i1 + m_strideX - 1;
}
size_t Surface::get_j2() {
    return m_j1 + m_strideY - 1;
}
size_t Surface::get_k2() {
    return m_k1 + m_strideZ - 1;
}

void Surface::applyBoundaryConditions(real *dataField, FieldType fieldType, size_t level, bool sync) {
    // TODO(linh)
    // m_bdc_boundary->applyBoundaryCondition(dataField, indexFields, patch_starts, patch_ends, fieldType, level, sync);
}
