project(storm-util)
cmake_minimum_required(VERSION 2.8)
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${storm_SOURCE_DIR}/../cmake)

FIND_PACKAGE( Vigra 1.8.0 REQUIRED )
FIND_PACKAGE( HDF5 )


ADD_EXECUTABLE(sif2hdf5 EXCLUDE_FROM_ALL sif2hdf5.cpp)
ADD_EXECUTABLE(sif2tiff  sif2tiff.cpp)
ADD_EXECUTABLE(hdf52tiff EXCLUDE_FROM_ALL hdf52tiff.cpp )
ADD_EXECUTABLE(tiff2hdf5 EXCLUDE_FROM_ALL tiff2hdf5.cpp)
#ADD_DEPENDENCIES(experiments sif2hdf5 sif2tiff hdf52tiff tiff2hdf5)
TARGET_LINK_LIBRARIES(sif2hdf5 vigraimpex)
TARGET_LINK_LIBRARIES(sif2tiff vigraimpex)
TARGET_LINK_LIBRARIES(hdf52tiff vigraimpex)
TARGET_LINK_LIBRARIES(tiff2hdf5 vigraimpex)

include_directories(
    ${Vigra_INCLUDE_DIRS}
    )



#ADD_EXECUTABLE(conv_3d convert.cpp ../storm/myimportinfo.cpp )
#TARGET_LINK_LIBRARIES(conv_3d vigraimpex)

set(BIN_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS sif2hdf5 DESTINATION ${BIN_INSTALL_DIR})