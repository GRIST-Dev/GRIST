# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)
enable_language(Fortran)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/utils)
include_directories(${CMAKE_BINARY_DIR}/grist_lib/src)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/namelist)
include_directories(${CMAKE_BINARY_DIR}/src/infrastructure/core/setting)
set (DIR_SRCS grist_advection_module.f90  grist_hori_swe_module.f90
				grist_base_flux_module.f90  grist_limiter_module.f90
				grist_destruct_module.f90   grist_math_module.f90
				grist_ffsl_flux_module.f90  grist_mesh_weight_new.f90
				grist_flux_operators.f90    grist_recon_module.f90
				grist_gcd_module.f90        grist_util_module.f90
				grist_geometric_info.f90)
add_library(operators_lib ${DIR_SRCS})
TARGET_LINK_LIBRARIES(operators_lib utils_lib  setting_lib)
set (LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
