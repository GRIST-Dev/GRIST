# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)

enable_language(Fortran)

LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/utils)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/setting)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/operators)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/io)
include_directories(${CMAKE_BINARY_DIR}/grist_lib/src)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/namelist)

set (DIR_SRCS grist_config_mesh.f90 grist_grid_file_read.f90  grist_grid_file_vars.f90)
add_library(grid_lib ${DIR_SRCS})
TARGET_LINK_LIBRARIES(grid_lib io_lib)
set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
