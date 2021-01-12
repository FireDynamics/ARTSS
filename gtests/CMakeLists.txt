add_executable(artss_gtest ${ARTSS_SOURCE_FILES})

add_subdirectory(field)

target_include_directories(artss_gtest PRIVATE ${CMAKE_SOURCE_DIR})

target_link_libraries(artss_gtest gtest gtest_main)
