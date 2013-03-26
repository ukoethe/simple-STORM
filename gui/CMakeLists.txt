project(storm-gui)
cmake_minimum_required(VERSION 2.8)

set(STORMGUI_VERSION_MAJOR 0)
set(STORMGUI_VERSION_MINOR 3)
set(STORMGUI_VERSION_PATCH 0)

set(STORMGUI_AUTHORS "© 2013 Ilia Kats")

SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${storm-gui_SOURCE_DIR}/../cmake)

find_package(Qt4 4.6 COMPONENTS QtCore QtGui REQUIRED)

include(${QT_USE_FILE} )
add_definitions(${QT_DEFINITIONS})

if(WIN32)
    add_definitions(-D__WIN__)
endif(WIN32)

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
)

SET(UIS
    inputwidget.ui
    backgroundlevelgroupbox.ui
    advancedsettingsgroupbox.ui
    resultwidget.ui
)

if(WIN32 OR APPLE)
    qt4_add_resources(ICONS_RCC icons/icons.qrc)
endif(WIN32 OR APPLE)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/version.h @ONLY)

include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
)
get_property(INC DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
message(${INC})
qt4_wrap_ui(UI_HEADERS ${UIS})
add_executable(storm-gui ${SRCS} ${UI_HEADERS} ${ICONS_RCC})
set_target_properties(storm-gui PROPERTIES AUTOMOC TRUE)

target_link_libraries(storm-gui ${QT_LIBRARIES} wienerstorm vigraimpex)

IF(OPENMP_FOUND)
    set_target_properties(storm-gui PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
    set_target_properties(storm-gui PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
    add_definitions(-DOPENMP_FOUND)
ENDIF(OPENMP_FOUND)

set(BIN_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/bin)
set(DATA_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
install(TARGETS storm-gui DESTINATION ${BIN_INSTALL_DIR})
