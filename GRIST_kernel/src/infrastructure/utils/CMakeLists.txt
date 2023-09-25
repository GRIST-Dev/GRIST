# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)

enable_language(Fortran)

set (DIR_SRCS grist_constants.F90 grist_handle_error.f90 grist_list.f90 grist_svd_module.f90 grist_time_manager.f90)

add_library(utils_lib ${DIR_SRCS})
target_link_libraries(utils_lib lapack blas)
set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
