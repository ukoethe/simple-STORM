MACRO(COPY_TEST_DATA)
    FOREACH(test_data ${ARGN})
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${test_data}
                       ${CMAKE_CURRENT_BINARY_DIR}/${test_data}
                       COPYONLY)
    ENDFOREACH(test_data)
ENDMACRO(COPY_TEST_DATA)

project(test)

set(WIENERSTORM_SRCS
    stormparams.cpp
    dataparams.cpp
    wienerStorm.cpp
    util.cpp
)

set(FILES_TO_COPY
    test.hdf5
    test.tif
)

cmake_minimum_required(VERSION 2.8)
set(STORM_R_SCRIPT storm_fit.R)
if (NOT TARGET wienerstorm)
    add_resource(${STORM_R_SCRIPT} storm_fit.h storm_fit WIENERSTORM_SRCS)
    add_library(wienerstorm STATIC ${WIENERSTORM_SRCS})
    target_link_libraries(wienerstorm vigraimpex ${FFTW_LIBRARIES} ${R_LIBRARIES} ${RUDECONFIG_LIBRARIES})
endif() 



add_executable(test test.cpp)



# enable_testing()
# add_test(test test.cpp)
target_link_libraries(test wienerstorm)
install (TARGETS test DESTINATION ${BIN_INSTALL_DIR})

FOREACH(file ${FILES_TO_COPY})
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${file}
	    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
ENDFOREACH(file)