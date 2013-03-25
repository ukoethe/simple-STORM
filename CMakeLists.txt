project(wienerstorm)
cmake_minimum_required(VERSION 2.8)
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

FIND_PACKAGE( Vigra 1.9.0 REQUIRED )
FIND_PACKAGE( FFTW REQUIRED )
FIND_PACKAGE( HDF5 )
FIND_PACKAGE( OpenMP )
FIND_PACKAGE( R REQUIRED)
FIND_PACKAGE( RudeConfig REQUIRED)

IF(HDF5_FOUND)
    ADD_DEFINITIONS(-DHDF5_FOUND)
ELSE(HDF5_FOUND)
    message(WARNING "Compiling without HDF5. No hdf5-input will be possible")
ENDIF(HDF5_FOUND)


IF(HDF5_FOUND)
    ADD_DEFINITIONS(-DHDF5_FOUND)
ELSE(HDF5_FOUND)
    message(WARNING "Compiling without HDF5. No hdf5-input will be possible")
ENDIF(HDF5_FOUND)


include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/storm
    ${Vigra_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIR}
    ${R_INCLUDE_DIRS}
    ${RUDECONFIG_INCLUDE_DIR}
)

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -std=c++11"
)

set(STORM_R_SCRIPT_DIR ${CMAKE_INSTALL_PREFIX}/share/simpleStorm/)
set(STORM_R_SCRIPT storm_fit.R)
add_definitions(-DSTORM_RSCRIPT=\"${STORM_R_SCRIPT}\" -DSTORM_RSCRIPT_DIR=\"${STORM_R_SCRIPT_DIR}\")

add_subdirectory(storm)
add_subdirectory(gui)
