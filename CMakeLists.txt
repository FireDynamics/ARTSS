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
message( STATUS "CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR} )
message( STATUS "CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR} )
message( STATUS "CMAKE_EXE_OUTPUT_PATH: " ${EXECUTABLE_OUTPUT_PATH} )

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
message(STATUS "")

#statically link all gcc stuffs
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")

# add sources and headers
set(SOURCE_FILES
        "src/adaption/Adaption.cpp"
        "src/adaption/Layers.cpp"
        "src/adaption/Vortex.cpp"

        "src/advection/ExplicitAdvect.cpp"
        "src/advection/SLAdvect.cpp"

        "src/analysis/Analysis.cpp"
        "src/analysis/Solution.cpp"

        "src/boundary/Boundary.cpp"
        "src/boundary/BoundaryController.cpp"
        "src/boundary/BoundaryData.cpp"
        "src/boundary/BoundaryDataController.cpp"
        "src/boundary/Multigrid.cpp"
        "src/boundary/Obstacle.cpp"
        "src/boundary/Surface.cpp"

        "src/boundaryCondition/DomainBoundary.cpp"
        "src/boundaryCondition/ObstacleBoundary.cpp"

        "src/diffusion/ColoredGaussSeidelDiffuse.cpp"
        "src/diffusion/ExplicitDiffuse.cpp"
        "src/diffusion/JacobiDiffuse.cpp"

        "src/interfaces/PressureI.cpp"
        "src/interfaces/SolverI.cpp"
        "src/interfaces/SourceI.cpp"

        "src/pressure/VCycleMG.cpp"

        "src/solver/AdvectionSolver.cpp"
        "src/solver/AdvectionDiffusionSolver.cpp"
        "src/solver/DiffusionSolver.cpp"
        "src/solver/DiffusionTurbSolver.cpp"
        "src/solver/NSSolver.cpp"
        "src/solver/NSTurbSolver.cpp"
        "src/solver/NSTempSolver.cpp"
        "src/solver/NSTempConSolver.cpp"
        "src/solver/NSTempTurbSolver.cpp"
        "src/solver/NSTempTurbConSolver.cpp"
        "src/solver/PressureSolver.cpp"
        "src/solver/SolverSelection.cpp"

        "src/source/ExplicitEulerSource.cpp"

        "src/turbulence/ConstSmagorinsky.cpp"
        "src/turbulence/DynamicSmagorinsky.cpp"

        "src/utility/Parameters.cpp"
        "src/utility/tinyxml2.cpp"
        "src/utility/Utility.cpp"
        "src/utility/visit_writer.cpp"
        "src/utility/Visual.cpp"

        "src/Domain.cpp"
        "src/Field.cpp"
        "src/Functions.cpp"
        "src/TimeIntegration.cpp"
        )

set(HEADER_FILES
        "src/analysis/Analysis.h"
        "src/analysis/Solution.h"

        "src/adaption/Adaption.h"
        "src/adaption/Layers.h"
        "src/adaption/Vortex.h"

        "src/advection/ExplicitAdvect.h"
        "src/advection/SLAdvect.h"

        "src/boundary/Boundary.h"
        "src/boundary/BoundaryController.h"
        "src/boundary/BoundaryData.h"
        "src/boundary/BoundaryDataController.cpp"
        "src/boundary/Multigrid.h"
        "src/boundary/Obstacle.h"
        "src/boundary/Surface.h"

        "src/boundaryCondition/DomainBoundary.h"
        "src/boundaryCondition/ObstacleBoundary.h"

        "src/diffusion/ColoredGaussSeidelDiffuse.h"
        "src/diffusion/ExplicitDiffuse.h"
        "src/diffusion/JacobiDiffuse.h"

        "src/interfaces/AdaptionFunctionI.h"
        "src/interfaces/advectionI.h"
        "src/interfaces/diffusionI.h"
        "src/interfaces/PressureI.h"
        "src/interfaces/SolverI.h"
        "src/interfaces/SourceI.h"
        "src/interfaces/TurbulenceI.h"

        "src/pressure/VCycleMG.h"

        "src/solver/AdvectionSolver.h"
        "src/solver/AdvectionDiffusionSolver.h"
        "src/solver/DiffusionSolver.h"
        "src/solver/DiffusionTurbSolver.h"
        "src/solver/NSSolver.h"
        "src/solver/NSTurbSolver.h"
        "src/solver/NSTempSolver.h"
        "src/solver/NSTempConSolver.h"
        "src/solver/NSTempTurbSolver.h"
        "src/solver/NSTempTurbConSolver.h"
        "src/solver/PressureSolver.h"
        "src/solver/SolverSelection.h"

        "src/source/ExplicitEulerSource.h"

        "src/turbulence/ConstSmagorinsky.h"
        "src/turbulence/DynamicSmagorinsky.h"

        "src/utility/Parameters.h"
        "src/utility/tinyxml2.h"
        "src/utility/Tracer.h"
        "src/utility/Utility.h"
        "src/utility/visit_writer.h"
        "src/utility/Visual.h"

        "src/Domain.h"
        "src/Field.h"
        "src/Functions.h"
        "src/GlobalMacrosTyper.h"
        "src/TimeIntegration.h"
        )

# link cuda library for cuda features
link_directories( $ENV{CUDA_LIB} )
include_directories( $ENV{CUDA_INC} )

# find spdlog lib
if(NOT TARGET spdlog)
    # Stand-alone build
    find_package(spdlog)
endif()
if(NOT TARGET spdlog) # use local fallback
    add_subdirectory(external/spdlog)
endif()



#Targets
add_executable( artss_serial ${SOURCE_FILES} src/main.cpp)
target_link_libraries(artss_serial PRIVATE spdlog::spdlog)

add_executable(artss_serial_profile ${SOURCE_FILES} src/main.cpp)
set_target_properties(artss_serial_profile PROPERTIES COMPILE_FLAGS "-DPROFILING -DUSE_NVTX")
target_link_libraries(artss_serial_profile)
target_link_libraries(artss_serial_profile PRIVATE spdlog::spdlog)

add_executable(artss_multicore_cpu ${SOURCE_FILES} src/main.cpp)
set_target_properties(artss_multicore_cpu PROPERTIES COMPILE_FLAGS "-acc -ta=multicore" LINK_FLAGS "-acc -ta=multicore")
target_link_libraries(artss_multicore_cpu PRIVATE spdlog::spdlog)

add_executable(artss_multicore_cpu_profile ${SOURCE_FILES} src/main.cpp)
set_target_properties(artss_multicore_cpu_profile PROPERTIES COMPILE_FLAGS "-acc -ta=multicore -DPROFILING -DUSE_NVTX" LINK_FLAGS "-acc -ta=multicore")
target_link_libraries(artss_multicore_cpu_profile)
target_link_libraries(artss_multicore_cpu_profile PRIVATE spdlog::spdlog)

add_executable(artss_gpu ${SOURCE_FILES} src/main.cpp)
set_target_properties(artss_gpu PROPERTIES COMPILE_FLAGS "-acc -Minfo=accel -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA}" LINK_FLAGS "-acc -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA}")
target_link_libraries(artss_gpu PRIVATE spdlog::spdlog)

add_executable(artss_gpu_profile ${SOURCE_FILES} src/main.cpp)
set_target_properties(artss_gpu_profile PROPERTIES COMPILE_FLAGS "-acc -Minfo=accel -DPROFILING -DUSE_NVTX -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA}" LINK_FLAGS "-acc -ta=${GPU_CC_TA},lineinfo,${CUDA_VERSION_TA}")
target_link_libraries(artss_gpu_profile)
target_link_libraries(artss_gpu_profile PRIVATE spdlog::spdlog)

spdlog_enable_warnings(artss_serial)
spdlog_enable_warnings(artss_serial_profile)
spdlog_enable_warnings(artss_multicore_cpu)
spdlog_enable_warnings(artss_multicore_cpu_profile)
spdlog_enable_warnings(artss_gpu)
spdlog_enable_warnings(artss_gpu_profile)

enable_testing()
file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION ${CMAKE_BINARY_DIR} )

add_test( NAME AdvectionTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/advection )
add_test( NAME BurgersTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/burgers )
add_test( NAME DiffusionTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/diffusion )
add_test( NAME DiffusionHatTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/diffusion/hat )
add_test( NAME DiffusionTurbTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/diffusionTurb )
add_test( NAME DissipationTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/dissipation )
add_test( NAME PressureTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/pressure )
add_test( NAME NSMcDermottTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/mcDermott)
add_test( NAME NSVortexTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/vortex )
add_test( NAME NSCavityFlowTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/cavityFlow )
add_test( NAME NSChannelFlowTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/channelFlow )
add_test( NAME NSTurbVortexTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTurb/vortex )
add_test( NAME NSTurbMcDermottTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTurb/mcDermott )
add_test( NAME NSTempBuoyancyMMSTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTemp/mms )
add_test( NAME NSTempTurbBuoyancyMMSTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTempTurb/mms )
add_test( NAME NSTempTurbAdaptionTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTempTurb/dynamicBoundaries )

if(CMAKE_CXX_COMPILER_ID STREQUAL "PGI" OR CMAKE_C_COMPILER_ID STREQUAL "PGI")
  add_test( NAME AdvectionTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/advection )
  add_test( NAME BurgersTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/burgers )
  add_test( NAME DiffusionTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/diffusion )
  add_test( NAME DiffusionHatTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/diffusion/hat )
  add_test( NAME DiffusionTurbTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/diffusionTurb )
  add_test( NAME DissipationTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/dissipation )
  add_test( NAME PressureTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/pressure )
  add_test( NAME NSMcDermottTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/mcDermott)
  add_test( NAME NSVortexTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/vortex )
  add_test( NAME NSCavityFlowTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/cavityFlow )
  add_test( NAME NSChannelFlowTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/channelFlow )
  add_test( NAME NSTurbVortexTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTurb/vortex )
  add_test( NAME NSTurbMcDermottTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTurb/mcDermott )
  add_test( NAME NSTempBuoyancyMMSTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTemp/mms )
  add_test( NAME NSTempTurbBuoyancyMMSTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTempTurb/mms )
  add_test( NAME NSTempTurbAdaptionTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTempTurb/dynamicBoundaries )
  add_test( NAME AdvectionTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/advection )
  add_test( NAME BurgersTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/burgers )
  add_test( NAME DiffusionTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/diffusion )
  add_test( NAME DiffusionHatTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/diffusion/hat )
  add_test( NAME DiffusionTurbTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/diffusionTurb )
  add_test( NAME DissipationTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/dissipation )
  add_test( NAME PressureTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/pressure )
  add_test( NAME NSMcDermottTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/mcDermott)
  add_test( NAME NSVortexTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/vortex )
  add_test( NAME NSCavityFlowTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/cavityFlow )
  add_test( NAME NSChannelFlowTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokes/channelFlow )
  add_test( NAME NSTurbVortexTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTurb/vortex )
  add_test( NAME NSTurbMcDermottTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTurb/mcDermott )
  add_test( NAME NSTempBuoyancyMMSTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTemp/mms )
  add_test( NAME NSTempTurbBuoyancyMMSTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTempTurb/mms )
  add_test( NAME NSTempTurbAdaptionTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTempTurb/dynamicBoundaries )
endif()
