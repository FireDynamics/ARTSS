/// \file       Boundary.h
/// \brief      Data class of boundary object
/// \date       Oct 01, 2019
/// \author     My Linh Würzburger
/// \copyright  <2015-2020> Forschungszentrum Juelich GmbH. All rights reserved.

#ifndef ARTSS_BOUNDARY_BOUNDARY_H_
#define ARTSS_BOUNDARY_BOUNDARY_H_


#include <vector>
#include "BoundaryData.h"
#include "Obstacle.h"

class Boundary {
public:
    ~Boundary();
    Boundary(Obstacle** obstacleList, size_t numberOfObstacles, size_t size_obstacles, size_t level=0);
    explicit Boundary(size_t level=0);
    void init(size_t size_obstacles);

    size_t* getBoundaryList() {return m_boundaryList;};
    size_t* getBoundaryFront() {return m_boundaryFront;};
    size_t* getBoundaryBack() {return m_boundaryBack;};
    size_t* getBoundaryTop() {return m_boundaryTop;};
    size_t* getBoundaryBottom() {return m_boundaryBottom;};
    size_t* getBoundaryLeft() {return m_boundaryLeft;};
    size_t* getBoundaryRight() { return m_boundaryRight;};

    size_t* getInnerList() { return m_innerList;};
    size_t getSize_innerList() { return m_size_innerList; };

    size_t getSize_boundaryList() {return m_size_boundaryList;};
    size_t getSize_boundaryFront() {return  m_size_boundaryFront;};
    size_t getSize_boundaryBack() {return   m_size_boundaryBack;};
    size_t getSize_boundaryTop() {return    m_size_boundaryTop;};
    size_t getSize_boundaryBottom() {return m_size_boundaryBottom;};
    size_t getSize_boundaryLeft() {return   m_size_boundaryLeft;};
    size_t getSize_boundaryRight() { return m_size_boundaryRight;};

    void updateLists(Obstacle** obstacleList, size_t numberOfObstacles, size_t size_obstacles);
    void updateLists();
    void control(size_t size_obstacles);
private:
    size_t m_level;

    size_t* m_boundaryList;
    size_t* m_boundaryFront;
    size_t* m_boundaryBack;
    size_t* m_boundaryTop;
    size_t* m_boundaryBottom;
    size_t* m_boundaryLeft;
    size_t* m_boundaryRight;

    size_t m_size_boundaryList;
    size_t m_size_boundaryFront;
    size_t m_size_boundaryBack;
    size_t m_size_boundaryTop;
    size_t m_size_boundaryBottom;
    size_t m_size_boundaryLeft;
    size_t m_size_boundaryRight;

    size_t* m_innerList;
    size_t m_size_innerList;

    void boundaryCells();
    void innerCells(Obstacle** obstacleList, size_t numberOfObstacles);
    void innerCells();
    void print(size_t size_obstacles);
    void clearLists();
};


#endif /* ARTSS_BOUNDARY_BOUNDARY_H_ */
