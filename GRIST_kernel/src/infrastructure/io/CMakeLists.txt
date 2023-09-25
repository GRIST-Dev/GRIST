# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)
enable_language(Fortran)

include_directories(${CMAKE_BINARY_DIR}/grist_lib/src)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/operators)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/utils)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/setting)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/namelist)

set (DIR_SRCS grist_fileio_list_1d_module_par.f90  grist_wrap_nf.f90
				grist_fileio_0d_module.f90       grist_fileio_list_2d_module.f90      grist_wrap_pf.f90
				grist_fileio_0d_module_gcm.f90   grist_fileio_list_2d_module_par.f90  README
				grist_fileio_list_1d_module.f90  grist_list_array.f90)

add_library(io_lib ${DIR_SRCS})
target_link_libraries(io_lib operators_lib utils_lib setting_lib namelist_lib pnetcdf netcdf netcdff)
set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
