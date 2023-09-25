# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)

enable_language(Fortran)

set (DIR_SRCS grist_data_types.f90 grist_domain_types.f90 grist_element_types.f90)

add_library(setting_lib ${DIR_SRCS})
link_directories(${CMAKE_BINARY_DIR}/lib)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/core/utils)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/grist_lib/src)
target_link_libraries(setting_lib utils_lib grist_lib)
set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
