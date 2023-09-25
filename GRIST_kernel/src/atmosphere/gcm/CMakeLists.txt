# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)

enable_language(Fortran)

set (DIR_SRCS gcm_diagnose_module.f90  gcm_dyn_vars_module.f90  gcm_inout_module.f90      gcm_time_integration_rkfb.f90
				gcm_config_partition.f90  gcm_driver.f90           gcm_init_module.f90      gcm_phys_vars_module.f90  grist.F90)

link_directories(${CMAKE_BINARY_DIR}/lib)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/grist_lib/src)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/namelist)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/grid)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/hpe)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/core/io)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/core/setting)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/core/utils)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/core/operators)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/drivers/phys)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/src/infrastructure/nh-dynamics)

ADD_EXECUTABLE(gcrm.exe ${DIR_SRCS})
SET_TARGET_PROPERTIES(gcrm.exe PROPERTIES LINKER_LANGUAGE Fortran)
TARGET_LINK_LIBRARIES(gcrm.exe grist_lib namelist_lib grid_lib hpe_lib io_lib setting_lib phys_lib utils_lib nh_dynamics_lib operators_lib)
set (EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/run)
