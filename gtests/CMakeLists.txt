add_executable(artss_gtest ${ARTSS_SOURCE_FILES})
add_executable(artss_gtest_gpu ${ARTSS_SOURCE_FILES})

target_sources(artss_gtest
    PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/Environment.cpp
    )

add_subdirectory(domain)
add_subdirectory(field)
add_subdirectory(GPULists)
add_subdirectory(randomField)
add_subdirectory(utility)

target_compile_options(artss_gtest_gpu PUBLIC -acc -Minfo=accel -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA})
target_link_options(artss_gtest_gpu PUBLIC -acc -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA})

target_include_directories(artss_gtest PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(artss_gtest gtest gtest_main)
target_link_libraries(artss_gtest spdlog::spdlog)
target_link_libraries(artss_gtest fmt::fmt)

target_include_directories(artss_gtest_gpu PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(artss_gtest_gpu gtest gtest_main)
target_link_libraries(artss_gtest_gpu spdlog::spdlog)
target_link_libraries(artss_gtest_gpu fmt::fmt)

add_test(NAME artss_gtest COMMAND artss_gtest)
