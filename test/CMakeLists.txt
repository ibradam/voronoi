######################################################################
# CMakeLists.txt
#
# Description: configuration file for tests
# 
# Generated with: mmk 
# Author: ibrahim
#  
#   
# Date: Mon Nov 30 15:05:45 2015
# Changelog: 
#
######################################################################
include_directories (${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src)
link_libraries(${${PKG_NAME}_LIBRARIES})
######################################################################
## Build tests 
######################################################################
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test)

FILE(GLOB cppFiles RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

FOREACH(file ${cppFiles})
    STRING(REGEX REPLACE ".cpp\$" "" outfileName "${file}")
    add_executable(${PROJECT_NAME}_${outfileName} ${CMAKE_CURRENT_SOURCE_DIR}/${file} )
    add_test(${PROJECT_NAME}_${outfileName} ${CMAKE_BINARY_DIR}/test/${PROJECT_NAME}_${outfileName} )
ENDFOREACH(file)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY)
