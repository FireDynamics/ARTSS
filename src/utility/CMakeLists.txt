target_sources(artss
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/Algorithm.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Algorithm.h
        ${CMAKE_CURRENT_LIST_DIR}/Utility.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Utility.h
)

add_subdirectory(settings)
