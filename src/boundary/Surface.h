//
// Created by linh on 01.10.19.
//

#ifndef ARTSS_BOUNDARY_SURFACE_H_
#define ARTSS_BOUNDARY_SURFACE_H_

#include <vector>
#include "BoundaryData.h"

#include "../Domain.h"
#include "../utility/Utility.h"
#include "../utility/tinyxml2.h"

#include "BoundaryDataController.h"

class Surface {
 public:
    Surface(Settings const &settings, SurfaceSetting const &surfsettings);
    Surface(Settings const &settings, size_t surfaceID, size_t startIndex, size_t strideX, size_t strideY, size_t strideZ, size_t level);
    ~Surface();
    size_t* getSurfaceList() {return m_surfaceList;}
    size_t getSize_surfaceList() {return m_size_surfaceList;}

    size_t getStrideX() { return m_strideX;}
    size_t getStrideY() { return m_strideY;}
    size_t getStrideZ() { return m_strideZ;}

    size_t getSurfaceID() { return m_surfaceID;}

    void setBoundaryConditions(BoundarySetting const &boundary);

    void applyBoundaryConditions(real *dataField, FieldType FieldType, size_t level, bool sync);

    void print();

 private:
    BoundaryDataController* m_boundaryDataController;

    size_t m_surfaceID;

    size_t m_i1;
    size_t m_j1;
    size_t m_k1;

    size_t *m_surfaceList;  // indices of surface
    size_t m_size_surfaceList;

    size_t m_strideX;
    size_t m_strideY;
    size_t m_strideZ;

    std::vector<BoundaryData*> dataList;

    void init(size_t Nx, size_t Ny);
    void createSurface(size_t Nx, size_t Ny);

    size_t get_i2();
    size_t get_j2();
    size_t get_k2();

#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
};


#endif /* ARTSS_BOUNDARY_SURFACE_H_ */
