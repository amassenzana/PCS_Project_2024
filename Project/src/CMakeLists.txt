list(APPEND DFN_sources "${CMAKE_CURRENT_SOURCE_DIR}/DFNLibrary.cpp")
list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/DFNLibrary.hpp")
# list(APPEND DFN_include ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/PolygonalMesh.hpp")


set (DFN_sources ${DFN_sources} PARENT_SCOPE)
set (DFN_headers ${DFN_headers} PARENT_SCOPE)
set (DFN_include ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

