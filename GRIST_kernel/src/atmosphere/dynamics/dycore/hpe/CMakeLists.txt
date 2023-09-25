# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)
enable_language(Fortran)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/utils)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/setting)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/operators)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/io)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/namelist)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)

set (DIR_SRCS grist_hpe_constants.f90   grist_hpe_hydro_pgf.f90
				grist_hpe_continuity.f90  grist_hpe_vertial_advection.f90)
add_library(hpe_lib ${DIR_SRCS})
TARGET_LINK_LIBRARIES(hpe_lib io_lib)
set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
