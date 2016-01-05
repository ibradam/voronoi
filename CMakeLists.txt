######################################################################
# CMakeLists.txt
#
# Description: configuration file for package voronoi
# 
# Generated with: mmk 
# Author: ibrahim
#  
#   
# Date: Mon Nov 30 14:47:28 2015
# Changelog: 
#
######################################################################
cmake_minimum_required (VERSION 2.8)

if(PROJECT_NAME)
  set(PKG_USED ON)
endif(PROJECT_NAME)

project (voronoi)

set(${PROJECT_NAME}_VERSION_MAJOR 0)
set(${PROJECT_NAME}_VERSION_MINOR 0)
set(${PROJECT_NAME}_VERSION_PATCH 2)
set(${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_PATCH}.${${PROJECT_NAME}_VERSION_PATCH})

string (TOUPPER ${PROJECT_NAME} PKG_NAME)
if(PKG_USED)
 set(PKG "${PKG_NAME}_")
endif(PKG_USED)

set (PROJECT_VERSION    ${${PROJECT_NAME}_VERSION})
set (PROJECT_DATA_DIR ${PROJECT_SOURCE_DIR}/data)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_BINARY_DIR}/cmake/)

message ("== Configuring ${PROJECT_NAME} ${${PROJECT_NAME}_VERSION}")
######################################################################
# Compilation
######################################################################
set (SHARED OFF) 
set (STATIC ON)
set (INSTALL ON)

if (NOT PKG_USED)
  set (CMAKE_INSTALL_PREFIX "/usr/local/share/voronoi")
endif (NOT PKG_USED)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "")

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

if(APPLE)
   set(CMAKE_MACOSX_RPATH 1)
endif (APPLE)

option (MDEBUG "enable mdebug message" OFF)

if(MDEBUG)
   add_definitions ("-DMDEBUG_COUT")
endif()

enable_testing()

if (EXISTS ${PROJECT_SOURCE_DIR}/cmake/voronoi-config.h.in)
   configure_file("${PROJECT_SOURCE_DIR}/cmake/voronoi-config.h.in"
			  "${CMAKE_BINARY_DIR}/include/voronoi-config.h"
               		  @ONLY IMMEDIATE)
endif (EXISTS ${PROJECT_SOURCE_DIR}/cmake/voronoi-config.h.in)

include_directories ( ${CMAKE_BINARY_DIR}/include ${PROJECT_SOURCE_DIR}/src)

######################################################################
# Subdirectories
######################################################################
set (${PKG}DOC_DEFAULT OFF)
set (SUBDIR) 
# <subdirectories> #
set (SUBDIR ${SUBDIR} include)
set (SUBDIR ${SUBDIR} src)
set (SUBDIR ${SUBDIR} test)
set (SUBDIR ${SUBDIR} app)
# </subdirectories> #
 
foreach(dir ${SUBDIR})

string (TOUPPER ${dir} OPT)

if (NOT DEFINED ${PKG}${OPT})
  if (DEFINED ${PKG}${OPT}_DEFAULT)
     option (${PKG}${OPT} "add/skip subdirectory ${OPT}" ${${PKG}${OPT}_DEFAULT})
  else()
     option (${PKG}${OPT} "add/skip subdirectory ${OPT}" ON)
  endif()
endif()

if (${${PKG}${OPT}})
     add_subdirectory (${dir})
     message(" |-- ${dir} [${${PKG}${OPT}}]" )
else (${${PKG}${OPT}})
     message (" ... ${dir} [${OPT}=OFF]")
endif (${${PKG}${OPT}})

endforeach(dir)

######################################################################
# Update
######################################################################
add_custom_target(${PKG}update
                   COMMAND git pull origin master
                   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                   COMMENT "update the packages")

######################################################################
# Packaging
######################################################################
set(CPACK_PACKAGE_VENDOR "MyCompany")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_NAME})
set(CPACK_PACKAGE_DESCRIPTION ${PROJECT_NAME})
#set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/COPYING.txt)

set (CPACK_PACKAGE_VERSION ${${PROJECT_NAME}_VERSION})
set (CPACK_SOURCE_PACKAGE_FILE_NAME "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}-src")

  set (CPACK_PACKAGE_CONTACT "ibrahim@mycompany.com")
  set (CPACK_PACKAGE_SUMMARY ${PROJECT_NAME})
  set (CPACK_PACKAGE_DEFAULT_LOCATION "${PROJECT_SOURCE_DIR}")
  set (CPACK_SOURCE_INSTALLED_DIRECTORIES "${PROJECT_SOURCE_DIR};.")
  set (CPACK_BINARY_STGZ "ON")
  set (CPACK_BINARY_TBZ2 "OFF")
  set (CPACK_BINARY_TGZ "ON")
  set (CPACK_BINARY_TZ "ON")
  set (CPACK_GENERATOR "")
  set (CPACK_SOURCE_GENERATOR "TGZ")

include (CPack)

####################################################################
## Installation
####################################################################
if (BUILD_FOR_RELEASE) #--------------------------------------------

set(MMK_SCRIPT ${CMAKE_BINARY_DIR}/for_installer/bin/mmk)
file(WRITE "${MMK_SCRIPT}" "#!/bin/sh\n")
file(APPEND "${MMK_SCRIPT}" "SRC=[INSTALLER_TARGET_DIR]/share/voronoi\n")
file(APPEND "${MMK_SCRIPT}" "$SRC/bin/mmk.pl $*\n")

install(PROGRAMS ${CMAKE_BINARY_DIR}/for_installer/bin/mmk DESTINATION bin)

endif (BUILD_FOR_RELEASE) #-----------------------------------------

install(FILES CMakeLists.txt DESTINATION share/${PROJECT_NAME})

foreach (dir ${SUBDIR})
   install(DIRECTORY ${dir} DESTINATION share/${PROJECT_NAME})
endforeach (dir)

install(PROGRAMS bin/mmk.pl DESTINATION share/${PROJECT_NAME}/bin)
