######################################################################
# CMakeLists.txt
#
# Description: configuration file for libraries
# 
# Generated with: mmk 
# Author: ibrahim
#   
#  
# Date: Mon Nov 30 15:01:18 2015
# Changelog: 
#
######################################################################

######################################################################
# Library 
######################################################################
set(LIB_NAME voronoi)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

FILE(GLOB SRC_HEADERS "*.h*")

FILE(GLOB SRC_SOURCES "*.cpp")

IF(BUILD_FOR_RELEASE)
    add_library(${LIB_NAME} STATIC ${SRC_SOURCES})
    link_libraries(${LIB_NAME} )
    install(TARGETS ${LIB_NAME} DESTINATION lib )
ELSE(BUILD_FOR_RELEASE)
    IF(STATIC)
        IF(SHARED)
            add_library(${LIB_NAME}_static STATIC ${SRC_SOURCES})
            link_libraries(${LIB_NAME}_static )
            install(TARGETS ${LIB_NAME}_static DESTINATION lib )
        ELSE(SHARED)
            add_library(${LIB_NAME} STATIC ${SRC_SOURCES})
            link_libraries(${LIB_NAME} )
            install(TARGETS ${LIB_NAME} DESTINATION lib )
        ENDIF(SHARED)
    ENDIF(STATIC)

    IF(SHARED)
        add_library(${LIB_NAME} SHARED ${SRC_SOURCES})
        link_libraries(${LIB_NAME} )
        install(TARGETS ${LIB_NAME} DESTINATION lib )
    ENDIF(SHARED)
ENDIF(BUILD_FOR_RELEASE)

####################################################################
# Installation 
####################################################################
install(FILES ${SRC_HEADERS} DESTINATION include/${PROJECT_NAME})
install(FILES ${SRC_MODULES} DESTINATION modules)
install(TARGETS ${LIB_NAME}
   RUNTIME DESTINATION bin
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib)

