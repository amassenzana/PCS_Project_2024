list(APPEND DFN_TEST_sources "${CMAKE_CURRENT_SOURCE_DIR}/testDFN.cpp")
list(APPEND DFN_TEST_headers "${CMAKE_CURRENT_SOURCE_DIR}/testDFN.hpp")
list(APPEND DFN_TEST_include "${CMAKE_CURRENT_SOURCE_DIR}")

set (DFN_TEST_sources ${DFN_TEST_sources} PARENT_SCOPE)
set (DFN_TEST_headers ${DFN_TEST_headers} PARENT_SCOPE)
set (DFN_TEST_include ${DFN_TEST_include} PARENT_SCOPE)
