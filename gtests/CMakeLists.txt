add_executable(artss_gtest ${ARTSS_SOURCE_FILES})

target_sources(artss_gtest
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/Environment.cpp
    )

add_subdirectory(advection)
add_subdirectory(domain)
add_subdirectory(field)
add_subdirectory(GPULists)
add_subdirectory(randomField)
add_subdirectory(utility)

target_include_directories(artss_gtest PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(artss_gtest gtest gtest_main)
target_link_libraries(artss_gtest spdlog::spdlog)
target_link_libraries(artss_gtest fmt::fmt)

add_test(NAME artss_gtest COMMAND artss_gtest)
