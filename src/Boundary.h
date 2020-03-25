/*
 * Boundary.h
 *
 *  Created on: May 25, 2016
 *      Author: Severt
 */

/// \file 		Boundary.h
/// \brief 		Parses boundary conditions from xml, builds index lists, provides functions to apply BC
/// \date 		May 25, 2016
/// \author 	Severt
/// \copyright 	<2015-2018> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef SRC_BOUNDARY_H_
#define SRC_BOUNDARY_H_

#include <vector>

#include "GlobalMacrosTypes.h"
#include "tinyxml2.h"

const size_t nTypes = 6;
enum class Types : size_t { RHO = 0, U = 1, V = 2, W = 3, P = 4, T = 5 };

const size_t nPatches = 6;
enum class Patches : size_t {FRONT = 0, BACK = 1, TOP = 2, BOTTOM = 3, LEFT = 4, RIGHT = 5};

const size_t nBCTypes = 3;
enum class BCTypes : size_t {NEUMANN = 0, DIRICHLET = 1, PERIODIC = 2};

struct BoundaryInfo {
    Types field;
    Patches patchID;
    BCTypes type;
    real value;
};

struct SurfaceInfo {
    size_t SurfaceID;
    Types field;
    Patches patchID;
    BCTypes type;
    real value;
};

struct ObstacleInfo {
    size_t ObstacleID;
    Types field;
    Patches patchID;
    BCTypes type;
    real value;
};

class Boundary {

public:
    Boundary();
    Boundary(tinyxml2::XMLElement xmlParameter);
    ~Boundary();

    static Boundary* getInstance();
    static void resetInstance();

    void parseParameter(tinyxml2::XMLElement *xmlParameter);
    void printBoundaries();

    std::string getTypeName(Types f);

    void updateListsDynamical(bool sync = true);		// for dynamical domain extension
	void addCGSMGLists();		// for colored gauss seidel (multigrid)
	void createCGSList(size_t Nx, size_t Ny); // for colored gauss seidel (divide in red and black list)

    std::vector<size_t> coordinateFromLinearIndex(size_t idx, size_t Nx, size_t Ny);

    void applyBoundary(real *d, size_t Nx, size_t Ny, size_t Nz, Types f, real dx, real dy, real dz, bool sync = true);
    void applyBoundary(real *d, size_t level, size_t Nx, size_t Ny, size_t Nz, Types f, real dx, real dy, real dz, bool sync = true); // for MG
    void applyBoundary(real *d, size_t Nx, size_t Ny, size_t Nz, Types f, real dx, real dy, real dz, real* val, bool sync = true); // for non-const BC

// GETTER
    // get lists/ vectors
    // boundary
    std::vector<std::vector<size_t>> get_bList() const {return this->m_bList;};
    std::vector<size_t> get_bList_joined() const {return this->m_bList_joined;};
    std::vector<size_t> get_bFront() const {return this->m_bList[0];};
    std::vector<size_t> get_bBack() const {return this->m_bList[1];};
    std::vector<size_t> get_bTop() const {return this->m_bList[2];};
    std::vector<size_t> get_bBottom() const {return this->m_bList[3];};
    std::vector<size_t> get_bLeft() const {return this->m_bList[4];};
    std::vector<size_t> get_bRight() const {return this->m_bList[5];};

    // surface
    std::vector<std::vector<size_t>> get_sList() const {return this->m_sList;};
    std::vector<size_t> get_sList_joined() const {return this->m_sList_joined;};

    // obstacle
    std::vector<std::vector<size_t>> get_oList() const {return this->m_oList;};
    std::vector<size_t> get_oList_joined() const {return this->m_oList_joined;};
    std::vector<std::vector<size_t>> get_oFronts() const {return this->m_oFront;};
    std::vector<size_t> get_oFronts_joined() const {return this->m_oFront_joined;};
    std::vector<std::vector<size_t>> get_oBacks() const {return this->m_oBack;};
    std::vector<size_t> get_oBacks_joined() const {return this->m_oBack_joined;};
    std::vector<std::vector<size_t>> get_oTops() const {return this->m_oTop;};
    std::vector<size_t> get_oTops_joined() const {return this->m_oTop_joined;};
    std::vector<std::vector<size_t>> get_oBottoms() const {return this->m_oBottom;};
    std::vector<size_t> get_oBottoms_joined() const {return this->m_oBottom_joined;};
    std::vector<std::vector<size_t>> get_oLefts() const {return this->m_oLeft;};
    std::vector<size_t> get_oLefts_joined() const {return this->m_oLeft_joined;};
    std::vector<std::vector<size_t>> get_oRights() const {return this->m_oRight;};
    std::vector<size_t> get_oRights_joined() const {return this->m_oRight_joined;};
    std::vector<std::vector<size_t>> get_oInners() const {return this->m_oInner;};

    // inner and non-calculated
    std::vector<size_t> get_iList() const {return this->m_iList;};
	std::vector<size_t> get_iList_cg_red() const {return m_iList_cg_red;};
	std::vector<size_t> get_iList_cg_black() const {return m_iList_cg_black;};
	std::vector<size_t> get_nList() const {return this->m_nList;};

    // inner && boundary
    std::vector<size_t> get_ibList_joined() const{return this->m_ibList;};

    // get data pointer
    size_t* get_d_bList_joined() const {return this->m_d_bList_joined;};

    size_t* get_d_sList_joined() const {return this->m_d_sList_joined;};

    size_t* get_d_oList_joined() const {return this->m_d_oList_joined;};
	size_t* get_d_oFronts_joined() const {return this->m_d_oFront_joined;};
	size_t* get_d_oBacks_joined() const {return this->m_d_oBack_joined;};
	size_t* get_d_oTops_joined() const {return this->m_d_oTop_joined;};
	size_t* get_d_oBottoms_joined() const {return this->m_d_oBottom_joined;};
	size_t* get_d_oLefts_joined() const {return this->m_d_oLeft_joined;};
	size_t* get_d_oRights_joined() const {return this->m_d_oRight_joined;};

    size_t* get_d_iList() const {return this->m_d_iList;};
	size_t* get_d_iList_cg_red() const {return m_d_iList_cg_red;};
	size_t* get_d_iList_cg_black() const {return m_d_iList_cg_black;};

    // get lists/vectors for MG (only those necessary)
    std::vector<std::vector<size_t>> get_MG_bList_joined() const {return this->m_MG_bList_joined;};
    std::vector<std::vector<size_t>> get_MG_iList() const {return this->m_MG_iList;};
	size_t* get_d_MG_iList_level_joined_CG_red() const {return m_d_MG_iList_level_joined_CG_red;};
	size_t* get_d_MG_iList_level_joined_CG_black() const {return m_d_MG_iList_level_joined_CG_black;};

    // get data pointer
    size_t* get_d_MG_bList_level_joined() const {return this->m_d_MG_bList_level_joined;};
    size_t* get_d_MG_iList_level_joined() const {return this->m_d_MG_iList_level_joined;};
	std::vector<std::vector<size_t>> get_MG_iList_CG_red() const {return m_MG_iList_CG_red;};
	std::vector<std::vector<size_t>> get_MG_iList_CG_black() const {return m_MG_iList_CG_black;};

    // getter for vector sizes
    size_t get_bList_Size() const {return this->m_bList_size;};
    size_t get_oList_Size() const {return this->m_oList_size;};
    size_t get_iList_Size() const {return this->m_iList_size;};

    size_t get_MGbList_Size(int level) const {return m_MG_bList_joined[level].size();};
    size_t get_MGiList_Size(int level) const {return m_MG_iList[level].size();};
    size_t get_MGiList_Length() const {return m_MG_iList_level_joined.size();};


private:
    static Boundary* single; //Singleton
	void differentiateCGList(std::vector<size_t> list, std::vector<size_t>* red, std::vector<size_t>* black, size_t Nx, size_t Ny);

    BoundaryInfo* getBoundaryInfo(Types f, Patches p);
    SurfaceInfo* getSurfaceInfo(size_t id, Types f, Patches p);
    ObstacleInfo* getObstacleInfo(size_t id, Types f, Patches p);

    void calcLists(bool sync = true);
    void addMGLists(size_t nlevel, bool sync);		// for multigrid (dominant, obst not too thin else stop restriction)
    void GetDuplicates(size_t Nx, size_t Ny, size_t nx, size_t ny, size_t nz);
    void GetMGDuplicates(size_t level, size_t Nx, size_t Ny, size_t nx, size_t ny, size_t nz);
    void SendBoundaryListsToGPU();
    void SendMGBoundaryListsToGPU();

    void applyBCtoBoundary(real *d, size_t Nx, size_t Ny, size_t Nz, size_t nx, size_t ny, size_t nz, Types f, real dx, real dy, real dz, bool sync = true);
    void applyBCtoBoundary(real *d, size_t level, size_t Nx, size_t Ny, size_t Nz, size_t nx, size_t ny, size_t nz, Types f, real dx, real dy, real dz, bool sync = true);
    void applyBCtoBoundary(real *d, size_t Nx, size_t Ny, size_t Nz, size_t nx, size_t ny, size_t nz, Types f, real dx, real dy, real dz, real* val, bool sync = true);

    void applyBCtoSurface(real *d, size_t id, size_t Nx, size_t Ny, size_t Nz, size_t stride_x, size_t stride_y, size_t stride_z, Types f, real dx, real dy, real dz, bool sync = true);
    void applyBCtoSurface(real *d, size_t level, size_t id, size_t Nx, size_t Ny, size_t Nz, size_t stride_x, size_t stride_y, size_t stride_z, Types f, real dx, real dy, real dz, bool sync = true);
    void applyBCtoSurface(real *d, size_t id, size_t Nx, size_t Ny, size_t Nz, size_t stride_x, size_t stride_y, size_t stride_z, Types f, real dx, real dy, real dz, real* val, bool sync = true);

    void applyBCtoObstacle(real *d, size_t id, size_t Nx, size_t Ny, size_t Nz, size_t stride_x, size_t stride_y, size_t stride_z, Types f, real dx, real dy, real dz, bool sync = true);
	void applyBCtoObstacle(real * d, size_t level, size_t id, size_t Nx, size_t Ny, size_t Nz, size_t stride_x, size_t stride_y, size_t stride_z, Types f, real dx, real dy, real dz, bool sync = true);
	void applyBCtoObstacle(real *d, size_t id, size_t Nx, size_t Ny, size_t Nz, size_t stride_x, size_t stride_y, size_t stride_z, Types f, real dx, real dy, real dz, real* val, bool sync = true);

// MEMBER
	// vectors
    std::vector<BoundaryInfo*> m_boundaryInfos;
    std::vector<std::vector<SurfaceInfo*>> m_surfaceInfos;
    std::vector<std::vector<ObstacleInfo*>> m_obstacleInfos;

    std::vector<std::vector<size_t>> m_bList;	// boundary indices (Vector of 6 vectors <front, back, top, bottom, left, right> including corners and edges
    std::vector<size_t> m_bList_joined;
    std::vector<size_t> m_bFront, m_bBack, m_bTop, m_bBottom, m_bLeft, m_bRight;
    std::vector<size_t> m_bFrontDuplicates, m_bBackDuplicates, m_bTopDuplicates, m_bBottomDuplicates, m_bLeftDuplicates, m_bRightDuplicates;

    std::vector<std::vector<size_t>> m_sList;	// surface indices of multiple surfaces <s1,s2,s3,...> within boundary
    std::vector<size_t> m_sList_joined;

    std::vector<std::vector<size_t>> m_oList;	// obstacle indices of multiple obstacles <o1,o2,o3,...> within inner cells
	std::vector<size_t> m_oList_joined;
	std::vector<std::vector<size_t>> m_oFront;	// all fronts of all obstacles
	std::vector<size_t> m_oFront_joined;
	std::vector<std::vector<size_t>> m_oBack;	// all backs of all obstacles
	std::vector<size_t> m_oBack_joined;
	std::vector<std::vector<size_t>> m_oTop;	// all tops of all obstacles
	std::vector<size_t> m_oTop_joined;
	std::vector<std::vector<size_t>> m_oBottom;	// all bottoms of all obstacles
	std::vector<size_t> m_oBottom_joined;
	std::vector<std::vector<size_t>> m_oLeft;	// all lefts of all obstacles
	std::vector<size_t> m_oLeft_joined;
	std::vector<std::vector<size_t>> m_oRight;	// all rights of all obstacles
	std::vector<size_t> m_oRight_joined;
	std::vector<size_t> m_oFrontDuplicates_joined, m_oBackDuplicates_joined, m_oTopDuplicates_joined, m_oBottomDuplicates_joined, m_oLeftDuplicates_joined, m_oRightDuplicates_joined;
	std::vector<std::vector<size_t>> m_oInner;	// all inner of all obstacles

    std::vector<size_t> m_iList;				// inner cells indices
    std::vector<size_t> m_iList_cg_red; // red indices for cgs
    std::vector<size_t> m_iList_cg_black; // black indices for cgs
    std::vector<size_t> m_ibList;				// inner and boundary cells indices
    std::vector<size_t> m_nList;				// non-calculated indices

    // data pointer
    size_t* m_d_bList_joined;
	size_t* m_d_bFront;
	size_t* m_d_bBack;
	size_t* m_d_bTop;
	size_t* m_d_bBottom;
	size_t* m_d_bLeft;
	size_t* m_d_bRight;
	size_t* m_d_bFrontDuplicates, *m_d_bBackDuplicates, *m_d_bTopDuplicates, *m_d_bBottomDuplicates, *m_d_bLeftDuplicates, *m_d_bRightDuplicates;

    size_t* m_d_sList_joined;

    size_t* m_d_oList_joined;
    std::vector<size_t>* m_d_oFront;
    size_t* m_d_oFront_joined;
    size_t* m_d_oBack_joined;
    size_t* m_d_oTop_joined;
    size_t* m_d_oBottom_joined;
    size_t* m_d_oLeft_joined;
    size_t* m_d_oRight_joined;
    size_t* m_d_oFrontDuplicates_joined, *m_d_oBackDuplicates_joined, *m_d_oTopDuplicates_joined, *m_d_oBottomDuplicates_joined, *m_d_oLeftDuplicates_joined, *m_d_oRightDuplicates_joined;

    size_t* m_d_iList;
    size_t* m_d_iList_cg_red;
    size_t* m_d_iList_cg_black;

    // Lists for Multigrid
    std::vector<std::vector<std::vector<size_t>>> m_MG_bList; 	// boundary indices for each Multigrid level (level 0 = MG_bList[0]=bList)
	std::vector<std::vector<size_t>> m_MG_bList_joined; 		// all boundaries for level l
	std::vector<size_t> m_MG_bList_level_joined; 				// all boundaries for all levels
	std::vector<std::vector<size_t>> m_MG_bFront, m_MG_bBack, m_MG_bTop, m_MG_bBottom, m_MG_bLeft, m_MG_bRight;
	std::vector<size_t> m_MG_bFront_level_joined, m_MG_bBack_level_joined, m_MG_bTop_level_joined, m_MG_bBottom_level_joined, m_MG_bLeft_level_joined, m_MG_bRight_level_joined;
	std::vector<size_t> m_MG_bFrontDuplicates_joined, m_MG_bBackDuplicates_joined, m_MG_bTopDuplicates_joined, m_MG_bBottomDuplicates_joined, m_MG_bLeftDuplicates_joined, m_MG_bRightDuplicates_joined;

    std::vector<std::vector<std::vector<size_t>>> m_MG_sList; 	// surface indices for each Multigrid level (level 0 = MG_sList[0]=sList)
    std::vector<std::vector<size_t>> m_MG_sList_joined; 		// all surfaces for level l
    std::vector<size_t> m_MG_sList_level_joined; 				// all surfaces for all levels

    std::vector<std::vector<std::vector<size_t>>> m_MG_oList; 	// obstacle indices for each Multigrid level (level 0 = MGoList[0]=oList)
    std::vector<std::vector<size_t>> m_MG_oList_joined;
    std::vector<std::vector<std::vector<size_t>>> m_MG_oFront; 	// front obstacle indices for each Multigrid level (level 0 = MGoFront[0]=oFront)
    std::vector<std::vector<std::vector<size_t>>> m_MG_oBack; 	// back obstacle indices for each Multigrid level (level 0 = MGoBack[0]=oBack)
    std::vector<std::vector<std::vector<size_t>>> m_MG_oTop; 	// top obstacle indices for each Multigrid level (level 0 = MGoTop[0]=oTop)
    std::vector<std::vector<std::vector<size_t>>> m_MG_oBottom; // bottom obstacle indices for each Multigrid level (level 0 = MGoBottom[0]=oBottom)
    std::vector<std::vector<std::vector<size_t>>> m_MG_oLeft; 	// left obstacle indices for each Multigrid level (level 0 = MGoLeft[0]=oLeft)
    std::vector<std::vector<std::vector<size_t>>> m_MG_oRight; 	// right obstacle indices for each Multigrid level (level 0 = MGoRight[0]=oRight)
    std::vector<std::vector<std::vector<size_t>>> m_MG_oInner; 	// inner obstacle indices for each Multigrid level (level 0 = MGoInner[0]=oInner)
    std::vector<std::vector<size_t>> m_MG_oFront_joined, m_MG_oBack_joined, m_MG_oTop_joined, m_MG_oBottom_joined, m_MG_oLeft_joined, m_MG_oRight_joined;
    std::vector<size_t> m_MG_oFront_level_joined, m_MG_oBack_level_joined, m_MG_oTop_level_joined, m_MG_oBottom_level_joined, m_MG_oLeft_level_joined, m_MG_oRight_level_joined;
    std::vector<size_t> m_MG_oFrontDuplicates_joined, m_MG_oBackDuplicates_joined, m_MG_oTopDuplicates_joined, m_MG_oBottomDuplicates_joined, m_MG_oLeftDuplicates_joined, m_MG_oRightDuplicates_joined;

    std::vector<std::vector<size_t>> m_MG_iList;				// inner cells indices for each Multigrid level (level 0 = MG_iList[0]=iList)
    std::vector<size_t> m_MG_iList_level_joined; 				// all inner cells for all levels
    std::vector<std::vector<size_t>> m_MG_ibList;
    std::vector<std::vector<size_t>> m_MG_nList;	            // non-calculated indices for each Multigrid level (level 0 = MG_nList[0]=nList)
    std::vector<std::vector<size_t>> m_MG_iList_CG_red;
    std::vector<std::vector<size_t>> m_MG_iList_CG_black;
    std::vector<size_t> m_MG_iList_level_joined_CG_red;
    std::vector<size_t> m_MG_iList_level_joined_CG_black;

    // data pointer for MG lists
    size_t* m_d_MG_sList_level_joined;
    size_t* m_d_MG_bList_level_joined;
    size_t* m_d_MG_bFront_level_joined;
    size_t* m_d_MG_bBack_level_joined;
    size_t* m_d_MG_bTop_level_joined;
    size_t* m_d_MG_bBottom_level_joined;
    size_t* m_d_MG_bLeft_level_joined;
    size_t* m_d_MG_bRight_level_joined;
    size_t* m_d_MG_bFrontDuplicates_joined, *m_d_MG_bBackDuplicates_joined, *m_d_MG_bTopDuplicates_joined, *m_d_MG_bBottomDuplicates_joined, *m_d_MG_bLeftDuplicates_joined, *m_d_MG_bRightDuplicates_joined;
    size_t* m_d_MG_oFront_level_joined;
    size_t* m_d_MG_oBack_level_joined;
    size_t* m_d_MG_oTop_level_joined;
    size_t* m_d_MG_oBottom_level_joined;
    size_t* m_d_MG_oLeft_level_joined;
    size_t* m_d_MG_oRight_level_joined;
    size_t* m_d_MG_oFrontDuplicates_joined, *m_d_MG_oBackDuplicates_joined, *m_d_MG_oTopDuplicates_joined, *m_d_MG_oBottomDuplicates_joined, *m_d_MG_oLeftDuplicates_joined, *m_d_MG_oRightDuplicates_joined;
    size_t* m_d_MG_iList_level_joined;
    size_t* m_d_MG_iList_level_joined_CG_red;
    size_t* m_d_MG_iList_level_joined_CG_black;

    // sizes
    size_t m_bList_size = 0;
	size_t m_bFront_size= 0;
	size_t m_bBack_size= 0;
	size_t m_bTop_size= 0;
	size_t m_bBottom_size= 0;
	size_t m_bLeft_size= 0;
	size_t m_bRight_size= 0;

    size_t m_oList_size = 0;

    size_t m_iList_size = 0;
};

#endif /* SRC_BOUNDARY_H_ */
