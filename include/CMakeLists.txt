####################################################################
## Subdirs
####################################################################
set (SUBDIR) 
# <subdirectories> #
set (SUBDIR ${SUBDIR} voronoi)
# </subdirectories> #

foreach(OPT ${SUBDIR})
     add_subdirectory (${OPT})
endforeach(OPT)

FILE(GLOB INCLUDE_HEADERS "*.h*")

add_custom_target (include_headers SOURCES ${INCLUDE_HEADERS})  
####################################################################
## Installation
####################################################################
install (FILES ${INCLUDE_HEADERS} DESTINATION include/${PROJECT_NAME})
