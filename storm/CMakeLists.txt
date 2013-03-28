project(storm)
cmake_minimum_required(VERSION 2.8)

set(STORM_VERSION_MAJOR 0)
set(STORM_VERSION_MINOR 7)
set(STORM_VERSION_PATCH 0)

set(STORM_AUTHORS "© 2011-2013 Ullrich Koethe, Joachim Schleicher, Frank Herrmannsdörfer, and Ilia Kats")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/version.h @ONLY)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

include(ResourceCompiler)

set(STORM_R_SCRIPT storm_fit.R)

set(WIENERSTORM_SRCS
    stormparams.cpp
    dataparams.cpp
    wienerStorm.cpp
    util.cpp
)
set(WIENERSTORM_HDRS
    stormparams.h
    dataparams.h
    wienerStorm.hxx
    util.h
)

if (CMAKE_CXX_COMPILER_ID STREQUAL MSVC)
    set(WIENERSTORM_SRCS
        ${WIENERSTORM_SRCS}
        getoptMSVC.
    )
endif (CMAKE_CXX_COMPILER_ID STREQUAL MSVC)

if(WIN32)
    add_definitions(-D__WIN__)
endif(WIN32)

add_resource(${STORM_R_SCRIPT} storm_fit.h storm_fit WIENERSTORM_SRCS)

add_library(wienerstorm STATIC ${WIENERSTORM_SRCS})
target_link_libraries(wienerstorm vigraimpex ${FFTW_LIBRARIES} ${R_LIBRARIES} ${RUDECONFIG_LIBRARIES})

add_executable(storm storm.cpp)
target_link_libraries(storm wienerstorm vigraimpex)

if(WIN32)
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")
endif(WIN32)

IF(OPENMP_FOUND)
	SET_TARGET_PROPERTIES(wienerstorm PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
	SET_TARGET_PROPERTIES(wienerstorm PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	SET_TARGET_PROPERTIES(storm PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
	SET_TARGET_PROPERTIES(storm PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
    add_definitions(-DOPENMP_FOUND)
ENDIF(OPENMP_FOUND)

IF(HDF5_FOUND)
    TARGET_LINK_LIBRARIES(wienerstorm ${HDF5_LIBRARIES})
    INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )
ENDIF(HDF5_FOUND)

ADD_SUBDIRECTORY(test)
install(TARGETS storm DESTINATION ${BIN_INSTALL_DIR})
