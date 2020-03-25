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
set(GPU_MODEL "${GPU_MODEL}")
# Set CUDA Flag
set(CUDA_VERSION "${CUDA_VERSION}")

if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if (NOT GPU_MODEL)
    message(STATUS "Setting GPU Model to default: P100")
    set(GPU_MODEL "P100")
  endif()

  if( ${GPU_MODEL} STREQUAL "K40" OR ${GPU_MODEL} STREQUAL "K80" )
    set(GPU_MODEL_TA "tesla:cc35")
  else( ${GPU_MODEL} STREQUAL "P100")
    set(GPU_MODEL_TA "tesla:cc60")
  endif()
  message(STATUS "GPU Model: ${GPU_MODEL}")

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
  set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} ${warning_flags} -O0 -g -Wunused-variable -Wunused-parameter -Wunused-function -Wshadow -Wunused -DTRACE_LOGGING")
  set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -w -DNDEBUG")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI" OR CMAKE_C_COMPILER_ID STREQUAL "PGI")
  set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -std=c++17 -Mnollvm")
  set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
  set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -fast -w -DNDEBUG")
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
        "Src/Adaption/Adaption.cpp"
        "Src/Adaption/Layers.cpp"
        "Src/Adaption/Vortex.cpp"

        "Src/Advection/ExplicitAdvect.cpp"
        "Src/Advection/SLAdvect.cpp"

        "Src/Analysis/Analysis.cpp"
        "Src/Analysis/Solution.cpp"

        "Src/Boundary/Boundary.cpp"
        "Src/Boundary/BoundaryController.cpp"
        "Src/Boundary/BoundaryData.cpp"
        "Src/Boundary/BoundaryDataController.cpp"
        "Src/Boundary/Multigrid.cpp"
        "Src/Boundary/Obstacle.cpp"
        "Src/Boundary/Surface.cpp"

        "Src/BoundaryCondition/DomainBoundary.cpp"
        "Src/BoundaryCondition/ObstacleBoundary.cpp"

        "Src/Diffusion/ColoredGaussSeidelDiffuse.cpp"
        "Src/Diffusion/ExplicitDiffuse.cpp"
        "Src/Diffusion/JacobiDiffuse.cpp"

        "Src/Interfaces/PressureI.cpp"
        "Src/Interfaces/SolverI.cpp"
        "Src/Interfaces/SourceI.cpp"

        "Src/Pressure/VCycleMG.cpp"

        "Src/Solver/AdvectionSolver.cpp"
        "Src/Solver/AdvectionDiffusionSolver.cpp"
        "Src/Solver/DiffusionSolver.cpp"
        "Src/Solver/DiffusionTurbSolver.cpp"
        "Src/Solver/NSSolver.cpp"
        "Src/Solver/NSTurbSolver.cpp"
        "Src/Solver/NSTempSolver.cpp"
        "Src/Solver/NSTempConSolver.cpp"
        "Src/Solver/NSTempTurbSolver.cpp"
        "Src/Solver/NSTempTurbConSolver.cpp"
        "Src/Solver/PressureSolver.cpp"
        "Src/Solver/SolverSelection.cpp"

        "Src/Source/ExplicitEulerSource.cpp"

        "Src/Turbulence/ConstSmagorinsky.cpp"
        "Src/Turbulence/DynamicSmagorinsky.cpp"

        "Src/Utility/Parameters.cpp"
        "Src/Utility/tinyxml2.cpp"
        "Src/Utility/Utility.cpp"
        "Src/Utility/visit_writer.cpp"
        "Src/Utility/Visual.cpp"

        "Src/Domain.cpp"
        "Src/Field.cpp"
        "Src/Functions.cpp"
        "Src/TimeIntegration.cpp"
        )

set(HEADER_FILES
        "Src/Analysis/Analysis.h"
        "Src/Analysis/Solution.h"

        "Src/Adaption/Adaption.h"
        "Src/Adaption/Layers.h"
        "Src/Adaption/Vortex.h"

        "Src/Advection/ExplicitAdvect.h"
        "Src/Advection/SLAdvect.h"

        "Src/Boundary/Boundary.h"
        "Src/Boundary/BoundaryController.h"
        "Src/Boundary/BoundaryData.h"
        "Src/Boundary/BoundaryDataController.cpp"
        "Src/Boundary/Multigrid.h"
        "Src/Boundary/Obstacle.h"
        "Src/Boundary/Surface.h"

        "Src/BoundaryCondition/DomainBoundary.h"
        "Src/BoundaryCondition/ObstacleBoundary.h"

        "Src/Diffusion/ColoredGaussSeidelDiffuse.h"
        "Src/Diffusion/ExplicitDiffuse.h"
        "Src/Diffusion/JacobiDiffuse.h"

        "Src/Interfaces/AdaptionFunctionI.h"
        "Src/Interfaces/AdvectionI.h"
        "Src/Interfaces/DiffusionI.h"
        "Src/Interfaces/PressureI.h"
        "Src/Interfaces/SolverI.h"
        "Src/Interfaces/SourceI.h"
        "Src/Interfaces/TurbulenceI.h"

        "Src/Pressure/VCycleMG.h"

        "Src/Solver/AdvectionSolver.h"
        "Src/Solver/AdvectionDiffusionSolver.h"
        "Src/Solver/DiffusionSolver.h"
        "Src/Solver/DiffusionTurbSolver.h"
        "Src/Solver/NSSolver.h"
        "Src/Solver/NSTurbSolver.h"
        "Src/Solver/NSTempSolver.h"
        "Src/Solver/NSTempConSolver.h"
        "Src/Solver/NSTempTurbSolver.h"
        "Src/Solver/NSTempTurbConSolver.h"
        "Src/Solver/PressureSolver.h"
        "Src/Solver/SolverSelection.h"

        "Src/Source/ExplicitEulerSource.h"

        "Src/Turbulence/ConstSmagorinsky.h"
        "Src/Turbulence/DynamicSmagorinsky.h"

        "Src/Utility/Parameters.h"
        "Src/Utility/tinyxml2.h"
        "Src/Utility/Tracer.h"
        "Src/Utility/Utility.h"
        "Src/Utility/visit_writer.h"
        "Src/Utility/Visual.h"

        "Src/Domain.h"
        "Src/Field.h"
        "Src/Functions.h"
        "Src/GlobalMacrosTyper.h"
        "Src/TimeIntegration.h"
        )

# link cuda library for cuda features
link_directories( $ENV{CUDA_LIB} )
include_directories( $ENV{CUDA_INC} )

#Targets
add_executable( artss_serial ${SOURCE_FILES} Src/main.cpp)

add_executable( artss_serial_profile ${SOURCE_FILES} Src/main.cpp  )
set_target_properties( artss_serial_profile PROPERTIES COMPILE_FLAGS "-DPROFILING -DUSE_NVTX" )
target_link_libraries( artss_serial_profile  )

add_executable( artss_multicore_cpu ${SOURCE_FILES} Src/main.cpp  )
set_target_properties( artss_multicore_cpu PROPERTIES COMPILE_FLAGS "-acc -ta=multicore" LINK_FLAGS "-acc -ta=multicore" )

add_executable( artss_multicore_cpu_profile ${SOURCE_FILES} Src/main.cpp  )
set_target_properties( artss_multicore_cpu_profile PROPERTIES COMPILE_FLAGS "-acc -ta=multicore -DPROFILING -DUSE_NVTX" LINK_FLAGS "-acc -ta=multicore" )
target_link_libraries( artss_multicore_cpu_profile  )

add_executable( artss_gpu ${SOURCE_FILES} Src/main.cpp  )
set_target_properties( artss_gpu PROPERTIES COMPILE_FLAGS "-acc -Minfo=accel -ta=${GPU_MODEL_TA},lineinfo,${CUDA_VERSION_TA}" LINK_FLAGS "-acc -ta=${GPU_MODEL_TA},lineinfo,${CUDA_VERSION_TA}" )

add_executable( artss_gpu_profile ${SOURCE_FILES} Src/main.cpp  )
set_target_properties( artss_gpu_profile PROPERTIES COMPILE_FLAGS "-acc -Minfo=accel -DPROFILING -DUSE_NVTX -ta=${GPU_MODEL_TA},lineinfo,${CUDA_VERSION_TA}" LINK_FLAGS "-acc -ta=${GPU_MODEL_TA},lineinfo,${CUDA_VERSION_TA}" )
target_link_libraries( artss_gpu_profile  )

enable_testing()
file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/Test DESTINATION ${CMAKE_BINARY_DIR} )

add_test( NAME AdvectionTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Advection )
add_test( NAME BurgersTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Burgers )
add_test( NAME DiffusionTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Diffusion )
add_test( NAME DiffusionHatTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Diffusion/Hat )
add_test( NAME DiffusionTurbTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/DiffusionTurb )
add_test( NAME DissipationTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Dissipation )
add_test( NAME PressureTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Pressure )
add_test( NAME NSMcDermottTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/McDermott)
add_test( NAME NSVortexTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/Vortex )
add_test( NAME NSCavityFlowTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/CavityFlow )
add_test( NAME NSChannelFlowTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/ChannelFlow )
add_test( NAME NSTurbVortexTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTurb/Vortex )
add_test( NAME NSTurbMcDermottTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTurb/McDermott )
add_test( NAME NSTempBuoyancyMMSTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTemp/MMS )
add_test( NAME NSTempTurbBuoyancyMMSTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTempTurb/MMS )

if(CMAKE_CXX_COMPILER_ID STREQUAL "PGI" OR CMAKE_C_COMPILER_ID STREQUAL "PGI")
  add_test( NAME AdvectionTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Advection )
  add_test( NAME BurgersTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Burgers )
  add_test( NAME DiffusionTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Diffusion )
  add_test( NAME DiffusionHatTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Diffusion/Hat )
  add_test( NAME DiffusionTurbTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/DiffusionTurb )
  add_test( NAME DissipationTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Dissipation )
  add_test( NAME PressureTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Pressure )
  add_test( NAME NSMcDermottTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/McDermott)
  add_test( NAME NSVortexTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/Vortex )
  add_test( NAME NSCavityFlowTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/CavityFlow )
  add_test( NAME NSChannelFlowTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/ChannelFlow )
  add_test( NAME NSTurbVortexTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTurb/Vortex )
  add_test( NAME NSTurbMcDermottTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTurb/McDermott )
  add_test( NAME NSTempBuoyancyMMSTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTemp/MMS )
  add_test( NAME NSTempTurbBuoyancyMMSTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTempTurb/MMS )
  add_test( NAME AdvectionTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Advection )
  add_test( NAME BurgersTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Burgers )
  add_test( NAME DiffusionTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Diffusion )
  add_test( NAME DiffusionHatTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Diffusion/Hat )
  add_test( NAME DiffusionTurbTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/DiffusionTurb )
  add_test( NAME DissipationTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Dissipation )
  add_test( NAME PressureTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/Pressure )
  add_test( NAME NSMcDermottTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/McDermott)
  add_test( NAME NSVortexTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/Vortex )
  add_test( NAME NSCavityFlowTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/CavityFlow )
  add_test( NAME NSChannelFlowTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokes/ChannelFlow )
  add_test( NAME NSTurbVortexTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTurb/Vortex )
  add_test( NAME NSTurbMcDermottTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTurb/McDermott )
  add_test( NAME NSTempBuoyancyMMSTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTemp/MMS )
  add_test( NAME NSTempTurbBuoyancyMMSTest_acc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/Test/NavierStokesTempTurb/MMS )
endif()
