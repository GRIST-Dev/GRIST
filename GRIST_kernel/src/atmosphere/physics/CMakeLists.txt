# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)

enable_language(Fortran)

include_directories(${CMAKE_BINARY_DIR}/grist_lib/src)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/operators)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/hpe)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/grid)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/utils)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/io)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/setting)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/namelist)

set (DIR_SRCS grist_phys_interface.F90)

add_library(phys_lib ${DIR_SRCS})
target_link_libraries(phys_lib io_lib  hpe_lib)
set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
