# - Try to find R
# Once done this will define
#
#  R_FOUND - system has OpenConnect
#  R_INCLUDE_DIRS - the OpenConnect include directories
#  R_LIBRARIES - the libraries needed to use OpenConnect
#  R_CFLAGS - Compiler switches required for using OpenConnect
#  R_VERSION - version number of OpenConnect



IF (R_INCLUDE_DIRS)
    # in cache already
    SET(R_FIND_QUIETLY TRUE)
ENDIF (R_INCLUDE_DIRS)

# use pkg-config to get the directories and then use these values
# in the FIND_PATH() and FIND_LIBRARY() calls
find_package(PkgConfig)
pkg_search_module(R libR)

IF (R_FOUND)
    IF (NOT R_FIND_QUIETLY)
        MESSAGE(STATUS "Found R ${R_VERSION}: ${R_LIBRARIES}")
    ENDIF (NOT R_FIND_QUIETLY)
ELSE (R_FOUND)
    IF (R_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could NOT find R, check FindPkgConfig output above!")
    ENDIF (R_FIND_REQUIRED)
ENDIF (R_FOUND)

MARK_AS_ADVANCED(R_INCLUDE_DIRS R_LIBRARIES R_STATIC_LIBRARIES)
