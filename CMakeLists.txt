# Documentation: Some useful options:
# -DCMAKE_C_COMPILER:STRING=pgcc -DCMAKE_CXX_COMPILER:STRING=pgc++ -DCMAKE_CXX_FLAGS:STRING="-acc -ta=tesla:managed -Minfo=all -fast -std=c++11"
# -DCMAKE_CXX_COMPILER:STRING=g++ -DCMAKE_CC_COMPILER:STRING=gcc
# -DCMAKE_CXX_COMPILER:STRING=clang++   -DCMAKE_C_COMPILER:STRING=clang
# -DCMAKE_BUILD_TYPE:STRING=Debug (default Release)
# -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON (default OFF)
# -DBUILD_TESTING=ON (default OFF)
# -DBUILD_CPPUNIT_TEST=ON (default OFF)
#--------------------------------------------------------------------------

cmake_minimum_required(VERSION 2.8  FATAL_ERROR)

project(ARTSS)
add_executable(artss "src/main.cpp")
set_target_properties(artss PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_subdirectory(src)
get_target_property(SOURCE_FILES artss SOURCES)

set(CMAKE_COLOR_MAKEFILE ON)

if(NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message( STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE} )

#------------------ set important directories --------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)
set(EXECUTABLE_OUTPUT_PATH bin)
set(LIBRARY_OUTPUT_PATH  bin/lib)

# Debug messages
message(STATUS "CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})
message(STATUS "CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR})
message(STATUS "CMAKE_EXE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH})

#----------------------------------------------------------------------
#---------------------------- compiler -------------------------------
if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(GPU_CC_TA "tesla:${GPU_CC}")
  message(STATUS "GPU Compute Compatiblity: ${GPU_CC}")

  if (NOT CUDA_VERSION)
    message(STATUS "Setting CUDA Version to default: 10.2")
    set(CUDA_VERSION "10.2")
  endif()

  if( ${CUDA_VERSION} STREQUAL "9" )
    set(CUDA_VERSION "9.0")
  elseif( ${CUDA_VERSION} STREQUAL "8" )
    set(CUDA_VERSION "8.0")
  elseif( ${CUDA_VERSION} STREQUAL "10" )
    set(CUDA_VERSION "10.2")
  endif()
  set(CUDA_VERSION_TA "cuda${CUDA_VERSION}" )
  message("-- CUDA Version: ${CUDA_VERSION}" )
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "PGI")
  set(CMAKE_CXX11_STANDARD_COMPILE_OPTION --c++11)
  set(CMAKE_CXX14_STANDARD_COMPILE_OPTION --c++14)
  set(CMAKE_CXX17_STANDARD_COMPILE_OPTION --c++17)
endif()

# Initialize CXXFLAGS.
target_compile_options(artss PUBLIC $<$<CONFIG:DEBUG>:-O0>)
target_compile_options(artss PUBLIC $<$<CONFIG:DEBUG>:-g>)
target_compile_options(artss PUBLIC $<$<CONFIG:DEBUG>:-Wall>)
target_compile_options(artss PUBLIC $<$<CONFIG:DEBUG>:-Wextra>)

target_compile_options(artss PUBLIC $<$<CONFIG:RELEASE>:-O3>)
target_compile_options(artss PUBLIC $<$<CONFIG:RELEASE>:-w>)
target_compile_options(artss PUBLIC $<$<CONFIG:RELEASE>:-DNDEBUG>)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  target_compile_options(artss PUBLIC $<$<CONFIG:DEBUG>:-Wno-cpp>)
  target_compile_options(artss PUBLIC $<$<CONFIG:DEBUG>:-Wno-unknown-pragmas>)
  target_compile_options(artss PUBLIC $<$<CONFIG:DEBUG>:-DTRACE_LOGGING>)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI" OR CMAKE_C_COMPILER_ID STREQUAL "PGI")
    target_compile_options(artss PUBLIC $<$<CONFIG:RELEASE>:-fast>)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
get_target_property(ARTSS_COMPILE_OPTIONS artss COMPILE_OPTIONS)


message("${COMPILER_SUPPORTS_CXX17}")

message(STATUS "C++ Compiler: "           ${CMAKE_CXX_COMPILER})
message(STATUS "Using ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_LDFLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
  message(STATUS "Debug flags: "  ${CMAKE_CXX_FLAGS_DEBUG} " " ${CMAKE_CXX_FLAGS})
else(CMAKE_BUILD_TYPE MATCHES Debug)
  message(STATUS "Release flags: " ${CMAKE_CXX_FLAGS_RELEASE} " " ${CMAKE_CXX_FLAGS})
endif(CMAKE_BUILD_TYPE MATCHES Debug)

# ========================== build platform ==========================
message(STATUS "")
message(STATUS "Platform: ")
message(STATUS "    Host: "             ${CMAKE_HOST_SYSTEM_NAME} ${CMAKE_HOST_SYSTEM_VERSION} ${CMAKE_HOST_SYSTEM_PROCESSOR})
if(CMAKE_CROSSCOMPILING)
  message(STATUS "    Target: "         ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION} ${CMAKE_SYSTEM_PROCESSOR})
endif()
message(STATUS "    CMake: "            ${CMAKE_VERSION})
message(STATUS "    CMake generator: "  ${CMAKE_GENERATOR})
message(STATUS "    CMake build tool: " ${CMAKE_BUILD_TOOL})
if(MSVC)
  message(STATUS "    MSVC: "           ${MSVC_VERSION})
endif()
if(CMAKE_GENERATOR MATCHES Xcode)
  message(STATUS "    Xcode: "          ${XCODE_VERSION})
endif()
if(NOT CMAKE_GENERATOR MATCHES "Xcode|Visual Studio")
  message(STATUS "    Configuration: "  ${CMAKE_BUILD_TYPE})
endif()
message(STATUS "")


# link cuda library for cuda features
link_directories($ENV{CUDA_LIB})
include_directories($ENV{CUDA_INC})

# find spdlog lib
if(NOT TARGET spdlog)
    # Stand-alone build
    find_package(spdlog QUIET)
endif()
# use local fallback
if(NOT TARGET spdlog)
    add_subdirectory(external/spdlog)
endif()

if(NOT TARGET fmt)
    # Stand-alone build
    find_package(fmt QUIET)
endif()
# use local fallback
if(NOT TARGET fmt)
    add_subdirectory(external/fmt)
endif()


#Targets
add_executable(artss_serial ${SOURCE_FILES} src/main.cpp)
target_compile_options(artss_serial PUBLIC ${ARTSS_COMPILE_OPTIONS})
target_link_libraries(artss_serial PRIVATE spdlog::spdlog)
target_link_libraries(artss_serial PRIVATE fmt::fmt)

add_executable(artss_serial_benchmark ${SOURCE_FILES})
target_compile_options(artss_serial_benchmark PUBLIC ${ARTSS_COMPILE_OPTIONS})
target_compile_options(artss_serial_benchmark PUBLIC -DBENCHMARKING -DUSE_NVTX)
target_link_libraries(artss_serial_benchmark)


add_executable(artss_multicore_cpu ${SOURCE_FILES})
target_compile_options(artss_multicore_cpu PUBLIC ${ARTSS_COMPILE_OPTIONS})
target_compile_options(artss_multicore_cpu PUBLIC -acc -ta=multicore)
target_link_options(artss_multicore_cpu PUBLIC -acc -ta=multicore)
target_link_libraries(artss_multicore_cpu PRIVATE spdlog::spdlog)

add_executable(artss_multicore_cpu_benchmark ${SOURCE_FILES})
target_compile_options(artss_multicore_cpu_benchmark PUBLIC ${ARTSS_COMPILE_OPTIONS})
target_compile_options(artss_multicore_cpu_benchmark PUBLIC -acc -ta=multicore -DBENCHMARKING -DUSE_NVTX)
target_link_options(artss_multicore_cpu_benchmark PUBLIC -acc -ta=multicore)
target_link_libraries(artss_multicore_cpu_benchmark)


add_executable(artss_gpu ${SOURCE_FILES})
target_compile_options(artss_gpu PUBLIC ${ARTSS_COMPILE_OPTIONS})
target_compile_options(artss_gpu PUBLIC -acc -Minfo=accel -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA})
target_link_options(artss_gpu PUBLIC -acc -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA})
target_link_libraries(artss_gpu PRIVATE spdlog::spdlog)

add_executable(artss_gpu_benchmark ${SOURCE_FILES})
target_compile_options(artss_gpu_benchmark PUBLIC ${ARTSS_COMPILE_OPTIONS})
target_compile_options(artss_gpu_benchmark PUBLIC -acc -Minfo=accel -DBENCHMARKING -DUSE_NVTX -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA})
target_link_options(artss_gpu_benchmark PUBLIC -acc -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA})
target_link_libraries(artss_gpu_benchmark)


enable_testing()
add_subdirectory(tests)
