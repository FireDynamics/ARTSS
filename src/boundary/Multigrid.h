/// \file 		Multigrid.h
/// \brief 		Creats all lists needed for multigrid
/// \date 		Oct 01, 2019
/// \author 	My Linh Würzburger
/// \copyright 	<2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_MULTIGRID_H_
#define ARTSS_BOUNDARY_MULTIGRID_H_

#include "Boundary.h"
#include "Surface.h"
#include "Obstacle.h"
#include "../Field.h"
#include "../Domain.h"
#include "../utility/Utility.h"
#include "BoundaryDataController.h"
#include <vector>


class Multigrid {
public:
    Multigrid(size_t numberOfSurfaces, Surface** surfaceList, size_t numberOfObstacles, Obstacle** obstacleList, BoundaryDataController* bdc_boundary, BoundaryDataController **bdc_obstacles);
    explicit Multigrid(BoundaryDataController *bdc_boundary);
    ~Multigrid();

    size_t getSize_innerList(size_t level = 0);
    size_t getSize_boundaryList(size_t level = 0);
    size_t getSize_obstacleList();
    size_t *get_obstacleList();

    size_t* getInnerList_level_joined() { return m_data_MG_iList_level_joined; };
    size_t getSize_innerList_level_joined() { return *(m_size_MG_iList_level + m_levels + 1); };
    size_t getInnerList_level_joined_start(size_t level);
    size_t getInnerList_level_joined_end(size_t level);

    size_t* getBoundaryList_level_joined() { return m_data_MG_bList_level_joined; };
    size_t getSize_boundaryList_level_joined() { return *(m_size_MG_bList_level + m_levels + 1); };
    size_t getBoundaryList_level_joined_start(size_t level);
    size_t getBoundaryList_level_joined_end(size_t level);

    void updateLists();

    void applyBoundaryCondition(real* d, size_t level, FieldType f, bool sync = false);

    size_t getObstacleStrideX(size_t id, size_t level);
    size_t getObstacleStrideY(size_t id, size_t level);
    size_t getObstacleStrideZ(size_t id, size_t level);

private:
#ifndef BENCHMARKING
    std::shared_ptr<spdlog::logger> m_logger;
#endif
    std::vector<BoundaryData*> m_boundaryData;
    size_t m_levels;
    // all surfaces divided by level
    Surface*** m_MG_surfaceList; //m_MG_surfaceList[level][surfaceID]
    size_t m_numberOfSurfaces = 0;
    // all obstacles divided by level
    Obstacle*** m_MG_obstacleList; //m_MG_obstacleList[level][obstacleID]
    size_t m_numberOfObstacles = 0;
    // boundary for each level
    Boundary** m_MG_boundaryList; //m_MG_boundaryList[level]

    //obstacle indices divided by level
    size_t** m_MG_oList;
    //surface indices divided by level
    size_t** m_MG_sList;

    //----- patches divided by level (obstacle indices)
    size_t** m_MG_oFront;
    size_t** m_MG_oBack;
    size_t** m_MG_oBottom;
    size_t** m_MG_oTop;
    size_t** m_MG_oLeft;
    size_t** m_MG_oRight;

    // start index of each level in level joined list
    size_t* m_size_MG_iList_level;
    size_t* m_size_MG_bList_level;

    // start index of each boundary object in level joined list (SliceZ = Front/Back, SliceY = Bottom/Top, SliceX = Left/Right)
    size_t* m_size_MG_bSliceZ_level;
    size_t* m_size_MG_bSliceY_level;
    size_t* m_size_MG_bSliceX_level;

    // start index of each obstacle object in level joined list
    size_t* m_size_MG_oFront_level;
    size_t* m_size_MG_oBack_level;
    size_t* m_size_MG_oTop_level;
    size_t* m_size_MG_oBottom_level;
    size_t* m_size_MG_oLeft_level;
    size_t* m_size_MG_oRight_level;


    size_t* m_size_MG_sList_level;
    size_t* m_size_MG_oList_level;

    //---- all level joined / arrays for GPU -----
    size_t* m_data_MG_iList_level_joined;
    size_t* m_data_MG_bList_level_joined;
    size_t* m_data_MG_sList_level_joined;

    size_t* m_data_MG_bFront_level_joined;
    size_t* m_data_MG_bBack_level_joined;
    size_t* m_data_MG_bTop_level_joined;
    size_t* m_data_MG_bBottom_level_joined;
    size_t* m_data_MG_bLeft_level_joined;
    size_t* m_data_MG_bRight_level_joined;

    size_t* m_data_MG_oFront_level_joined;
    size_t* m_data_MG_oBack_level_joined;
    size_t* m_data_MG_oTop_level_joined;
    size_t* m_data_MG_oBottom_level_joined;
    size_t* m_data_MG_oLeft_level_joined;
    size_t* m_data_MG_oRight_level_joined;
    size_t* m_data_MG_oList_zero_joined;

    size_t getSize_oList(size_t level);
    size_t getLastIndex_oFront( size_t level, size_t id);
    size_t getLastIndex_oBack( size_t level, size_t id);
    size_t getLastIndex_oBottom( size_t level, size_t id);
    size_t getLastIndex_oTop( size_t level, size_t id);
    size_t getLastIndex_oLeft( size_t level, size_t id);
    size_t getLastIndex_oRight( size_t level, size_t id);
    size_t getFirstIndex_oFront(size_t level, size_t id);
    size_t getFirstIndex_oBack( size_t level, size_t id);
    size_t getFirstIndex_oBottom( size_t level, size_t id);
    size_t getFirstIndex_oTop(size_t level, size_t id);
    size_t getFirstIndex_oLeft( size_t level, size_t id);
    size_t getFirstIndex_oRight( size_t level, size_t id);
    size_t getLen_oFront(size_t level);
    size_t getLen_oBack(size_t level);
    size_t getLen_oBottom(size_t level);
    size_t getLen_oTop(size_t level);
    size_t getLen_oLeft(size_t level);
    size_t getLen_oRight(size_t level);

    // get length of from/for joined array
    size_t getLen_oFront_joined();
    size_t getLen_oBack_joined();
    size_t getLen_oBottom_joined();
    size_t getLen_oTop_joined();
    size_t getLen_oLeft_joined();
    size_t getLen_oRight_joined();

    // get length of bSlice from/for joined array
    size_t getLen_bSliceZ_joined();
    size_t getLen_bSliceY_joined();
    size_t getLen_bSliceX_joined();

    size_t getFirstIndex_bSliceZ( size_t level);
    size_t getFirstIndex_bSliceX( size_t level);
    size_t getFirstIndex_bSliceY( size_t level);
    size_t getLastIndex_bSliceZ( size_t level);
    size_t getLastIndex_bSliceX( size_t level);
    size_t getLastIndex_bSliceY( size_t level);

    size_t getLen_iList_joined();
    size_t getLen_bList_joined();
    size_t getLen_sList_joined();
    size_t getFirstIndex_iList(size_t level);
    size_t getFirstIndex_bList(size_t level);
    size_t getFirstIndex_sList(size_t level);
    size_t getLastIndex_iList(size_t level);
    size_t getLastIndex_bList(size_t level);
    size_t getLastIndex_sList(size_t level);

    void init();
    void addMGLists();
    void calcObstacles(Obstacle** obstacleList);
    void calcSurfaces(Surface** surfaceList);
    void sendListsToGPU();
    void sendBoundaryListsToGPU();
    void sendSurfaceListsToGPU();
    void sendObstacleListsToGPU();

    void surfaceDominantRestriction(size_t level);
    Obstacle** obstacleDominantRestriction(size_t level);

    void control();
    void print();

    size_t** m_data_boundary_patches_joined;
    //size_t** m_data_surfaces_patches_joined;
    size_t** m_data_obstacles_patches_joined;
    BoundaryDataController* m_bdc_boundary;
    BoundaryDataController** m_bdc_obstacle;

    void removeBoundaryListsFromGPU();
};


#endif /* ARTSS_BOUNDARY_MULTIGRID_H_*/
