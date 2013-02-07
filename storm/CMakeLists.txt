project(storm)
cmake_minimum_required(VERSION 2.8)
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${storm_SOURCE_DIR}/../cmake)

FIND_PACKAGE( Vigra 1.8.0 REQUIRED )
FIND_PACKAGE( FFTW REQUIRED )
FIND_PACKAGE( HDF5 )
FIND_PACKAGE( OpenMP )
include_directories(

/home/herrmannsdoerfer/master/opengm-2.0/src/external/MaxFlow-v3.02.src-patched
/home/herrmannsdoerfer/master/opengm-2.0/src/external/MRF-v2.1.src-patched
/home/herrmannsdoerfer/master/opengm-2.0/include
/usr/share/R/include
)

set(CMAKE_CXX_FLAGS
	${CMAKE_CXX_FLAGS}
	-std=c++0x)

set(CMAKE_EXE_LINKER_FLAGS
	"${CMAKE_EXE_LINKER_FLAGS} -L/usr/lib/R/lib")

IF(HDF5_FOUND)
	ADD_DEFINITIONS(-DHDF5_FOUND)
ELSE(HDF5_FOUND)
	message(WARNING "Compiling without HDF5. No hdf5-input will be possible")
ENDIF(HDF5_FOUND)

IF(CMAKE_COMPILER_IS_GNUCXX)
	ADD_EXECUTABLE(storm storm.cpp program_options_getopt.cpp myimportinfo.cpp util.cpp)
	ADD_EXECUTABLE(wienerfilter EXCLUDE_FROM_ALL wienerfilter.cpp program_options_getopt.cpp myimportinfo.cpp util.cpp)
	target_link_libraries(storm R)
	target_link_libraries(wienerfilter R)
ELSE(CMAKE_COMPILER_IS_GNUCXX)
	ADD_DEFINITIONS(-DEMULATE_GETOPT)
	ADD_EXECUTABLE(storm storm.cpp getoptMSVC.c program_options_getopt.cpp myimportinfo.cpp util.cpp)
	ADD_EXECUTABLE(wienerfilter EXCLUDE_FROM_ALL wienerfilter.cpp getoptMSVC.c program_options_getopt.cpp myimportinfo.cpp util.cpp)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

IF(OPENMP_FOUND)
	SET_TARGET_PROPERTIES(storm PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
	SET_TARGET_PROPERTIES(storm PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
	SET_TARGET_PROPERTIES(wienerfilter PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
	SET_TARGET_PROPERTIES(wienerfilter PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
ENDIF(OPENMP_FOUND)

IF(HDF5_FOUND)
    TARGET_LINK_LIBRARIES(storm ${HDF5_LIBRARIES}
/home/herrmannsdoerfer/master/opengm-2.0/build/src/external/libexternal-library-maxflow.a
)
    TARGET_LINK_LIBRARIES(wienerfilter ${HDF5_LIBRARIES} 
/home/herrmannsdoerfer/master/opengm-2.0/build/src/external/libexternal-library-maxflow.a)
    INCLUDE_DIRECTORIES( ${HDF5_INCLUDE_DIRS} )
ENDIF(HDF5_FOUND)

TARGET_LINK_LIBRARIES(storm vigraimpex ${FFTW_LIBRARIES})
TARGET_LINK_LIBRARIES(wienerfilter vigraimpex ${FFTW_LIBRARIES})
include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${Vigra_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIRS}
    )

ADD_SUBDIRECTORY(test)
set(BIN_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS storm DESTINATION ${BIN_INSTALL_DIR})
