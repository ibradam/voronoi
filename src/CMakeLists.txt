######################################################################
# CMakeLists.txt
#
# Description: configuration file for libraries
# 
# Generated with: mmk 
# Author: ibrahim
#   
#  
# Date: Mon Nov 30 14:47:28 2015
# Changelog: 
#
######################################################################
include_directories (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)


######################################################################
# Subdirectories
######################################################################
set (SUBDIR)
# <subdirectories> #
set (SUBDIR ${SUBDIR} voronoi)
# </subdirectories> #

foreach(OPT ${SUBDIR})
   add_subdirectory (${OPT})
endforeach(OPT)

