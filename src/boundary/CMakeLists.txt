target_sources(artss
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Boundary.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Boundary.h
        ${CMAKE_CURRENT_LIST_DIR}/BoundaryController.cpp
        ${CMAKE_CURRENT_LIST_DIR}/BoundaryController.h
        ${CMAKE_CURRENT_LIST_DIR}/BoundaryData.cpp
        ${CMAKE_CURRENT_LIST_DIR}/BoundaryData.h
        ${CMAKE_CURRENT_LIST_DIR}/BoundaryDataController.cpp
        ${CMAKE_CURRENT_LIST_DIR}/BoundaryDataController.h
        ${CMAKE_CURRENT_LIST_DIR}/Multigrid.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Multigrid.h
        ${CMAKE_CURRENT_LIST_DIR}/Obstacle.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Obstacle.h
        ${CMAKE_CURRENT_LIST_DIR}/Surface.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Surface.h
)

