cmake_minimum_required (VERSION 2.8.10...3.20)

project (Libpll)

add_subdirectory(src)

set(PLL_EXPORTED_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/src PARENT_SCOPE)
