target_sources(artss
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/CSVWriter.cpp
        ${CMAKE_CURRENT_LIST_DIR}/CSVWriter.h
        ${CMAKE_CURRENT_LIST_DIR}/FieldWriter.cpp
        ${CMAKE_CURRENT_LIST_DIR}/FieldWriter.h
        ${CMAKE_CURRENT_LIST_DIR}/Visual.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Visual.h
        ${CMAKE_CURRENT_LIST_DIR}/VTKWriter.cpp
        ${CMAKE_CURRENT_LIST_DIR}/VTKWriter.h
        ${CMAKE_CURRENT_LIST_DIR}/visit_writer.cpp
        ${CMAKE_CURRENT_LIST_DIR}/visit_writer.h
)

