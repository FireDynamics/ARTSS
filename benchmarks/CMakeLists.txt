add_executable(artss_bench ${ARTSS_SOURCE_FILES})

add_subdirectory(field)

target_include_directories(artss_bench PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(artss_bench benchmark::benchmark)
target_link_libraries(artss_bench gtest)
target_link_libraries(artss_bench spdlog::spdlog)
target_link_libraries(artss_bench fmt::fmt)
target_link_libraries(artss_bench MPI::MPI_CXX)
target_link_libraries(artss_bench HighFive)

add_test(NAME artss_bench COMMAND artss_bench)
