/// \file 		BoundaryController.h
/// \brief 		Controll class for boundary
/// \date 		Oct 01, 2020
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_BOUNDARYCONTROLLER_H_
#define ARTSS_BOUNDARY_BOUNDARYCONTROLLER_H_


#include "../Field.h"
#include "Multigrid.h"
#include "Surface.h"
#include "Obstacle.h"
#include "Boundary.h"
#include "BoundaryData.h"
#include "../utility/tinyxml2.h"
#include "BoundaryDataController.h"

class BoundaryController {
public:
    static BoundaryController* getInstance();
    ~BoundaryController();

    void applyBoundary(real *d, FieldType f, bool sync = true);
    void applyBoundary(real *d, size_t level, FieldType f, bool sync = true);
    //void applyBoundary(real *d, size_t level, FieldType f, real* val, bool sync = true); // for non-const BC

    void printBoundaries();
    void updateLists();

    size_t getSize_innerList();
    size_t getSize_boundaryList();
    size_t* get_obstacleList();
    size_t getSize_obstacleList();

    size_t* get_innerList_level_joined();
    size_t getSize_innerList_level_joined(); //TODO necessary?
    size_t get_innerList_level_joined_start(size_t level);
    size_t get_innerList_level_joined_end(size_t level);

    size_t* get_boundaryList_level_joined();
    size_t getSize_boundaryList_level_joined(); //TODO necessary?
    size_t get_boundaryList_level_joined_start(size_t level);
    size_t get_boundaryList_level_joined_end(size_t level);

    size_t getSize_surfaceList() {return m_size_sList;};

    size_t getObstacleStrideX(size_t id, size_t level);
    size_t getObstacleStrideY(size_t id, size_t level);
    size_t getObstacleStrideZ(size_t id, size_t level);

    std::vector<FieldType> get_used_fields();

private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    static BoundaryController* singleton;

    BoundaryDataController *m_bdc_boundary;
    BoundaryDataController **m_bdc_obstacles;
    Multigrid* m_multigrid;

    Surface** m_surfaceList;
    size_t m_numberOfSurfaces = 0;
    Obstacle** m_obstacleList;
    size_t m_numberOfObstacles = 0;

    size_t m_size_sList = 0;

    bool m_hasObstacles;
    bool m_hasSurfaces;

    BoundaryController();
    void readXML();
    void parseBoundaryParameter(tinyxml2::XMLElement *xmlParameter);
    void parseObstacleParameter(tinyxml2::XMLElement *xmlParameter);
    void parseSurfaceParameter(tinyxml2::XMLElement *xmlParameter);
};


#endif /* ARTSS_BOUNDARY_BOUNDARYCONTROLLER_H_ */
