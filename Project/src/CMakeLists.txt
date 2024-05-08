list(APPEND DFN_sources "${CMAKE_CURRENT_SOURCE_DIR}/DFNLibrary.cpp")
list(APPEND DFN_headers "${CMAKE_CURRENT_SOURCE_DIR}/DFNLibrary.hpp")
list(APPEND DFN_include ${CMAKE_CURRENT_SOURCE_DIR})

set (DFN_sources ${DFN_sources} PARENT_SCOPE)
set (DFN_headers ${DFN_headers} PARENT_SCOPE)
set (DFN_include ${DFN_include} PARENT_SCOPE)
