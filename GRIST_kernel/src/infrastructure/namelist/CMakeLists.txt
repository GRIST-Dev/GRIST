# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)
enable_language(Fortran)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/utils)
include_directories(${CMAKE_BINARY_DIR}/grist_lib/src)

set (DIR_SRCS grist_nml_module.f90)

add_library(namelist_lib ${DIR_SRCS})
TARGET_LINK_LIBRARIES(namelist_lib utils_lib grist_lib)
set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
