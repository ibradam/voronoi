
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

file (GLOB CPP_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(cppfile ${CPP_FILES})
    string (REGEX REPLACE ".cpp\$" "" appname "${cppfile}")
    add_executable (${appname} ${cppfile})
    install (TARGETS ${appname} RUNTIME DESTINATION bin BUNDLE DESTINATION bin)
endforeach(cppfile)
