
file( COPY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${CMAKE_BINARY_DIR} )

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
  #add_test( NAME NSTempTurbAdaptionTest COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_serial
  #  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTempTurb/dynamicBoundaries )

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
  #  add_test( NAME NSTempTurbAdaptionTest_mc COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_multicore_cpu
  #    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTempTurb/dynamicBoundaries )
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
  #  add_test( NAME NSTempTurbAdaptionTest_gpu COMMAND sh ./run.sh ${CMAKE_BINARY_DIR}/${EXECUTABLE_OUTPUT_PATH}/artss_gpu
  #    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/navierStokesTempTurb/dynamicBoundaries )
endif()
