project(stormgui)
cmake_minimum_required(VERSION 2.8)

set(STORMGUI_VERSION_MAJOR 0)
set(STORMGUI_VERSION_MINOR 3)
set(STORMGUI_VERSION_PATCH 0)

set(STORMGUI_AUTHORS "© 2013 Ilia Kats Frank Herrmannsdoerfer und Ullrich Koethe")

set(GUI_EXECUTABLE storm-gui)
set(GUI_NAME simpleSTORM)

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

qt4_wrap_ui(UI_HEADERS ${UIS})
add_executable(${GUI_EXECUTABLE} WIN32 ${SRCS} ${UI_HEADERS} ${ICONS_RCC})
set_target_properties(${GUI_EXECUTABLE} PROPERTIES AUTOMOC TRUE)

target_link_libraries(${GUI_EXECUTABLE} ${QT_LIBRARIES} wienerstorm vigraimpex)

IF(OPENMP_FOUND)
    set_target_properties(${GUI_EXECUTABLE} PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
    set_target_properties(${GUI_EXECUTABLE} PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
    add_definitions(-DOPENMP_FOUND)
ENDIF(OPENMP_FOUND)

install(TARGETS ${GUI_EXECUTABLE} RUNTIME DESTINATION ${BIN_INSTALL_DIR} BUNDLE DESTINATION ${BIN_INSTALL_DIR})

set(CPACK_PACKAGE_VERSION_MAJOR ${STORMGUI_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${STORMGUI_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${STORMGUI_VERSION_PATCH})

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${QT_LIBRARY_DIR} ${QT_BINARY_DIR})
string(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_LOWER)
if (NOT "${BUILD_TYPE_LOWER}" STREQUAL "debug")
    set(CPACK_STRIP_FILES ON)
endif (NOT "${BUILD_TYPE_LOWER}" STREQUAL "debug")

option(BUILD_INSTALLER_PACKAGE "Build an installer." ON)
if (WIN32)
    option(BUILD_ZIP_PACKAGE "Build a zip package." ON)
endif (WIN32)


if (WIN32 OR APPLE)
    set(plugin_dest_dir ${BIN_INSTALL_DIR})
    string(LENGTH ${CMAKE_SHARED_LIBRARY_SUFFIX} SHARED_LIBRARY_SUFFIX_LENGTH)
    file(GLOB_RECURSE QT_PLUGINS_EXIST
    "${QT_PLUGINS_DIR}/imageformats/*${CMAKE_SHARED_LIBRARY_SUFFIX}")
    foreach(QT_PLUGIN ${QT_PLUGINS_EXIST})
        string(LENGTH ${QT_PLUGIN} QT_PLUGIN_LENGTH)
        math(EXPR QT_PLUGIN_SUBSTR_END ${QT_PLUGIN_LENGTH}-${SHARED_LIBRARY_SUFFIX_LENGTH}-1)
        string(SUBSTRING ${QT_PLUGIN} 0 ${QT_PLUGIN_SUBSTR_END} QT_PLUGIN_BASE)
        set(QT_PLUGIN_DEBUG "${QT_PLUGIN_BASE}d4${CMAKE_SHARED_LIBRARY_SUFFIX}")
        if (NOT ${QT_PLUGIN_DEBUG_FOUND} STREQUAL -1)
            if ("${BUILD_TYPE_LOWER}" STREQUAL "debug")
                set(QT_PLUGINS ${QT_PLUGINS} ${QT_PLUGIN_DEBUG})
            else ("${BUILD_TYPE_LOWER}" STREQUAL "debug")
                set(QT_PLUGINS ${QT_PLUGINS} ${QT_PLUGIN})
            endif ("${BUILD_TYPE_LOWER}" STREQUAL "debug")
        else(NOT ${QT_PLUGIN_DEBUG_FOUND} STREQUAL -1)
			set(QT_PLUGINS ${QT_PLUGINS} ${QT_PLUGIN})
        endif(NOT ${QT_PLUGIN_DEBUG_FOUND} STREQUAL -1)
    endforeach(QT_PLUGIN ${QT_PLUGINS_EXIST})
    install(CODE "
if (\${CMAKE_INSTALL_PREFIX} MATCHES \".*_CPack_Packages.*\")
file(INSTALL ${QT_PLUGINS} DESTINATION \${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/imageformats)
file(GLOB_RECURSE QTPLUGINS
\"\${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/imageformats/*${CMAKE_SHARED_LIBRARY_SUFFIX}\")
include(BundleUtilities)
fixup_bundle(\"\${CMAKE_INSTALL_PREFIX}/${BIN_INSTALL_DIR}/${GUI_EXECUTABLE}${CMAKE_EXECUTABLE_SUFFIX}\" \"\${QTPLUGINS}\" \"${CMAKE_LIBRARY_PATH};${CMAKE_SYSTEM_LIBRARY_PATH}\")
if (APPLE)
    file(WRITE \"\${CMAKE_INSTALL_PREFIX}/${BIN_INSTALL_DIR}/qt.conf\"
\"[Paths]
Plugins = .
\")
else (APPLE)
# on MacOS, the library is linked with a path relative to the executable
file(GLOB R_LIBS \"\${CMAKE_INSTALL_PREFIX}/\${BIN_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}R*${CMAKE_SHARED_LIBRARY_SUFFIX}\")
foreach(R_LIB_PATH \${R_LIBS})
    file(REMOVE \"\${R_LIB_PATH}\")
endforeach()
endif (APPLE)
endif (\${CMAKE_INSTALL_PREFIX} MATCHES \".*_CPack_Packages.*\")
")
endif(WIN32 OR APPLE)

set(CPACK_PACKAGE_INSTALL_DIRECTORY simpleSTORM)
set(CPACK_PACKAGE_EXECUTABLES
        ${GUI_EXECUTABLE} ${GUI_NAME}
)

if (BUILD_ZIP_PACKAGE)
	set(CPACK_GENERATOR
		${CPACK_GENERATOR}
		ZIP
	)
endif (BUILD_ZIP_PACKAGE)

if (WIN32 AND BUILD_INSTALLER_PACKAGE)
    set(CPACK_GENERATOR
        ${CPACK_GENERATOR}
        NSIS
    )
    set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON)
    set(CPACK_NSIS_COMPRESSOR lzma)
    set(CPACK_NSIS_MUI_FINISHPAGE_RUN "${GUI_EXECUTABLE}${CMAKE_EXECUTABLE_SUFFIX}")
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(CPACK_NSIS_PROGRAMFILES PROGRAMFILES64)
    else (CMAKE_SIZEOF_VOID_P EQUAL 8)
        set(CPACK_NSIS_PROGRAMFILES PROGRAMFILES)
    endif (CMAKE_SIZEOF_VOID_P EQUAL 8)
endif (WIN32 AND BUILD_INSTALLER_PACKAGE)

if (APPLE AND BUILD_INSTALLER_PACKAGE)
    set(CPACK_GENERATOR
        ${CPACK_GENERATOR}
        Bundle
    )
    set(CPACK_BUNDLE_NAME ${GUI_NAME})
    set(CPACK_BUNDLE_ICON ${CMAKE_CURRENT_SOURCE_DIR}/icons/StormIcon.icns)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Info.plist.in ${CMAKE_CURRENT_BINARY_DIR}/Info.plist @ONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/MacBundle.sh.in ${CMAKE_CURRENT_BINARY_DIR}/MacBundle.sh @ONLY)
    set(CPACK_BUNDLE_PLIST ${CMAKE_CURRENT_BINARY_DIR}/Info.plist)
    set(CPACK_BUNDLE_STARTUP_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/MacBundle.sh)
endif (APPLE AND BUILD_INSTALLER_PACKAGE)

include(CPack)
