# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)
enable_language(Fortran)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
include_directories(${CMAKE_BINARY_DIR}/grist_lib/src)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/grid)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/hpe)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/utils)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/setting)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/operators)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/io)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/namelist)
include_directories(${CMAKE_BINARY_DIR}/src/drivers/phys)

set (DIR_SRCS grist_config_partition.f90 grist_diffusion_module.f90         grist_nh_implicit_module.f90
				grist_nh_driver_module.f90         grist_tridiagonal_solver.f90
				grist_nh_explicit_tend_module.f90 ) 

add_library(nh_dynamics_lib ${DIR_SRCS})
TARGET_LINK_LIBRARIES(nh_dynamics_lib io_lib  hpe_lib grid_lib phys_lib)
set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
