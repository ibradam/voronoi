######################################################################
# CMakeLists.txt
#
# Description: configuration file for libraries
# 
# Generated with: mmk 
# Author: mourrain
#   
#  
# Date: Fri Dec  4 08:44:43 2015
# Changelog: 
#
######################################################################

FILE(GLOB INCLUDE_HEADERS "*.hpp*")

add_custom_target (${PROJECT_NAME}_include_headers SOURCES ${INCLUDE_HEADERS})

