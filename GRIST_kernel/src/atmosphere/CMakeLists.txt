# CMake 最低版本号要求
cmake_minimum_required (VERSION 2.8)

#项目名称
project (grist3)
enable_language(Fortran)

#子目录
add_subdirectory(gcrm)
add_subdirectory(phys)
add_subdirectory(swe)
