project(storm-gui)
cmake_minimum_required(VERSION 2.8)
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${storm-gui_SOURCE_DIR}/../cmake)

find_package(Qt4 4.6 COMPONENTS QtCore QtGui REQUIRED)
FIND_PACKAGE(Vigra 1.9.0 REQUIRED)
FIND_PACKAGE(FFTW REQUIRED)
FIND_PACKAGE( HDF5 )
FIND_PACKAGE( OpenMP )
FIND_PACKAGE( R REQUIRED)
FIND_PACKAGE( RudeConfig REQUIRED)

include(${QT_USE_FILE} )
add_definitions(${QT_DEFINITIONS})

IF(HDF5_FOUND)
ADD_DEFINITIONS(-DHDF5_FOUND)
ELSE(HDF5_FOUND)
message(WARNING "Compiling without HDF5. No hdf5-input will be possible")
ENDIF(HDF5_FOUND)

set(CMAKE_CXX_FLAGS
"${CMAKE_CXX_FLAGS} -std=c++11"
)

set(STORM_R_SCRIPT_DIR ${CMAKE_INSTALL_PREFIX}/share/simpleStorm/)
set(STORM_R_SCRIPT storm_fit.R)
add_definitions(-DSTORM_RSCRIPT=\"${STORM_R_SCRIPT}\" -DSTORM_RSCRIPT_DIR=\"${STORM_R_SCRIPT_DIR}\")

if(WIN32)
    add_definitions(-D__WIN__)
endif(WIN32)

# to be able to run the binaries without installing
file(COPY ../storm/${STORM_R_SCRIPT} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(SRCS
    analysisworker.cpp
    main.cpp
    mainwindow.cpp
    inputwidget.cpp
    resultwidget.cpp
    previewimage.cpp
    previewscrollarea.cpp
    filenamelineedit.cpp
    globals.cpp
    guiparams.cpp
    ../storm/stormparams.cpp
    ../storm/dataparams.cpp
    ../storm/wienerStorm.cpp
    ../storm/util.cpp
)

SET(UIS
    inputwidget.ui
    backgroundlevelgroupbox.ui
    advancedsettingsgroupbox.ui
    resultwidget.ui
)
include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../storm
    ${Vigra_INCLUDE_DIRS}
    ${FFTW_INCLUDE_DIRS}
    ${R_INCLUDE_DIRS}
    ${RUDECONFIG_INCLUDE_DIR}
    )
qt4_wrap_ui(UI_HEADERS ${UIS})
add_executable(storm-gui ${SRCS} ${UI_HEADERS})
set_target_properties(storm-gui PROPERTIES AUTOMOC TRUE)

target_link_libraries(storm-gui ${QT_LIBRARIES} vigraimpex ${FFTW_LIBRARIES} ${R_LIBRARIES} ${RUDECONFIG_LIBRARIES})
IF(HDF5_FOUND)
    target_link_libraries(storm-gui ${HDF5_LIBRARIES})
ENDIF(HDF5_FOUND)

set(BIN_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)
set(DATA_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
install(TARGETS storm-gui DESTINATION ${BIN_INSTALL_DIR})
install(FILES ${STORM_R_SCRIPT} DESTINATION ${STORM_R_SCRIPT_DIR})
