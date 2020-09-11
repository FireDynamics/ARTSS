# Documentation: Some useful options:
# -DCMAKE_C_COMPILER:STRING=pgcc -DCMAKE_CXX_COMPILER:STRING=pgc++ -DCMAKE_CXX_FLAGS:STRING="-acc -ta=tesla:managed -Minfo=all -fast -std=c++11"
# -DCMAKE_CXX_COMPILER:STRING=g++ -DCMAKE_CC_COMPILER:STRING=gcc
# -DCMAKE_CXX_COMPILER:STRING=clang++   -DCMAKE_C_COMPILER:STRING=clang
# -DCMAKE_BUILD_TYPE:STRING=Debug (default Release)
# -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON (default OFF)
# -DBUILD_TESTING=ON (default OFF)
# -DBUILD_CPPUNIT_TEST=ON (default OFF)
#--------------------------------------------------------------------------

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(ARTSS)

set(CMAKE_COLOR_MAKEFILE ON)

option(ARTSS_BENCHMARK_FLAG "Flag if BENCHMARKING is enabled" OFF)
option(ARTSS_TESTS_FLAG "Flag if testing is enabled" OFF)
option(ARTSS_UNITTEST_FLAG "Flag if unit testing is enabled" OFF)
set(ARTSS_TARGET_ARCH "SERIAL" CACHE STRING "Target arch (SERIAL, MULTICORE, GPU")

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
# enable as many warnings as possible
set(warning_flags "-g -Wall -Wextra")
# issue all the warnings demanded by strict iso c and iso c++
# set(warning_flags "${warning_flags} -pedantic")
# treat warnings as errors
# set(warning_flags "${warning_flags} -Werror")
# definition shadows another
# set(warning_flags "${warning_flags} -Wshadow")
# do not issue warnings for system headers
# set(warning_flags "${warning_flags} -Wno-system-headers")
# overloaded virtual function has a different signature
# set(warning_flags "${warning_flags} -Woverloaded-virtual")
# make string constants const char*
# set(warning_flags "${warning_flags} -Wwrite-strings")


#----------------------------------------------------------------------
#---------------------------- compiler -------------------------------

# Set compiler
set(CMAKE_CXX_COMPILER		"${CMAKE_CXX_COMPILER}")
set(CMAKE_CC_COMPILER		"${CMAKE_CC_COMPILER}")
# Set GPU Flag
set(GPU_CC "${GPU_CC}")
# Set CUDA Flag
set(CUDA_VERSION "${CUDA_VERSION}")

if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(GPU_CC_TA "tesla:${GPU_CC}")
  message(STATUS "GPU Compute Compatiblity: ${GPU_CC}")

  if (NOT CUDA_VERSION)
    message(STATUS "Setting CUDA Version to default: 8.0")
    set(CUDA_VERSION "8.0")
  endif()

  if( ${CUDA_VERSION} STREQUAL "9" )
    set(CUDA_VERSION "9.0")
  elseif( ${CUDA_VERSION} STREQUAL "8" )
    set(CUDA_VERSION "8.0")
  elseif( ${CUDA_VERSION} STREQUAL "10" )
    set(CUDA_VERSION "10.0")
  endif()
  set(CUDA_VERSION_TA "cuda${CUDA_VERSION}" )
  message("-- CUDA Version: ${CUDA_VERSION}" )
endif()

# Initialize CXXFLAGS.
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  add_compile_options(-Wno-unused-parameter)
  add_compile_options(-Wno-unused-variable)
  add_compile_options(-Wno-unknown-pragmas)
  add_compile_options(-Wno-cpp)
  add_compile_options(-Wno-shadow)

  set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} ${warning_flags} -O0 -g -Wunused-variable -Wunused-parameter -Wunused-function -Wshadow -Wunused -DTRACE_LOGGING")
  set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -w -DNDEBUG")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI" OR CMAKE_C_COMPILER_ID STREQUAL "PGI")
  set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -std=c++17 -Mnollvm")
  set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
  set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fast -w -DNDEBUG")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "PGI")
  set(CMAKE_CXX11_STANDARD_COMPILE_OPTION  --c++11)
  set(CMAKE_CXX14_STANDARD_COMPILE_OPTION  --c++14)
  set(CMAKE_CXX17_STANDARD_COMPILE_OPTION  --c++17)
 endif()

#--------------------
if(NOT CMAKE_GENERATOR MATCHES "Xcode|Visual Studio")
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if( COMPILER_SUPPORTS_CXX17)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++17 support. Please use a different C++ compiler.")
  endif()
endif()


#---------------------
message(STATUS "C++ Compiler: "           ${CMAKE_CXX_COMPILER})
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  message(STATUS "Using Clang " ${CMAKE_CXX_COMPILER_VERSION})
endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
if(CMAKE_COMPILER_IS_GNUCXX)
  message( STATUS "Using g++ " ${CMAKE_CXX_COMPILER_VERSION})
endif(CMAKE_COMPILER_IS_GNUCXX)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
  message(STATUS "Using PGI " ${CMAKE_CXX_COMPILER_VERSION})
endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_LDFLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fprofile-arcs -ftest-coverage")
  message(STATUS "Debug flags: "  ${CMAKE_CXX_FLAGS} " " ${CMAKE_CXX_FLAGS_DEBUG} )
else(CMAKE_BUILD_TYPE MATCHES Debug)
  message(STATUS "Release flags: " ${CMAKE_CXX_FLAGS} " ,  " ${CMAKE_CXX_FLAGS_RELEASE} )
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
message(STATUS "    ARTSS_BENCHMARK_FLAG:   ${ARTSS_BENCHMARK_FLAG}")
message(STATUS "    ARTSS_TARGET_ARCH:      ${ARTSS_TARGET_ARCH}")
message(STATUS "    ARTSS_TESTS_FLAG:       ${ARTSS_TESTS_FLAG}")
message(STATUS "    ARTSS_UNITTEST_FLAG:    ${ARTSS_UNITTEST_FLAG}")
message(STATUS "")


# link cuda library for cuda features
link_directories( $ENV{CUDA_LIB} )
include_directories( $ENV{CUDA_INC} )


message(STATUS "")


#Targets
add_executable(artss "src/main.cpp")


if(NOT ARTSS_BENCHMARK_FLAG)
    # find spdlog lib
    if(NOT TARGET spdlog)
        add_subdirectory(external/spdlog)
    endif()

    target_link_libraries(artss PRIVATE spdlog::spdlog)
    spdlog_enable_warnings(artss)
endif()


set(OUTPUT_FILE "artss")

if(ARTSS_TARGET_ARCH STREQUAL "MULTICORE")
    message(STATUS "Enabled Multicore")
    set_target_properties(artss PROPERTIES
        COMPILE_FLAGS "-acc -ta=multicore"
        LINK_FLAGS "-acc -ta=multicore")
    set(OUTPUT_FILE "${OUTPUT_FILE}_multicore_cpu")
elseif(ARTSS_TARGET_ARCH STREQUAL "GPU")
    message(STATUS "Enabled GPU")
    set_target_properties(artss PROPERTIES
        COMPILE_FLAGS "-acc -Minfo=accel -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA}"
        LINK_FLAGS "-acc -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA}")
    set(OUTPUT_FILE "${OUTPUT_FILE}_gpu")
else()
    message(STATUS "Using serial version")
    set(OUTPUT_FILE "${OUTPUT_FILE}_serial")
endif()

if(ARTSS_BENCHMARK_FLAG)
    message(STATUS "Enabled Benchmarking")
    target_compile_options(artss PRIVATE "-DBENCHMARKING -DUSE_NVTX")
    set(OUTPUT_FILE "${OUTPUT_FILE}_benchmark")
endif()


string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE_LOWER)
set(OUTPUT_FILE "${OUTPUT_FILE}_${BUILD_TYPE_LOWER}")
add_custom_command(TARGET artss POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:artss> ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/${OUTPUT_FILE})


add_subdirectory(src)

if(ARTSS_TESTS_FLAG)
    enable_testing()
    add_subdirectory(tests)
endif()
