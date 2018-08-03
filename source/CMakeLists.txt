cmake_minimum_required(VERSION 3.4)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE on)

# Only set project name if built as stand-alone library
if("${PROJECT_NAME}" STREQUAL "")
    project("vkm_main")
endif()


# OpenMesh -------------------------------
include_directories(${VKM_OPENMESH_SOURCE_DIR}/src)
link_directories(${VKM_OPENMESH_BINARY_DIR}/Build/lib)


# VXL ------------------------------------
set(VKM_VXL_SOURCE_DIR ${VKM_EXTERNAL_DIR}/vxl)
set(VKM_VXL_BINARY_DIR ${CMAKE_BINARY_DIR}/vxl-build)

set(VXL_USE_DCMTK "FALSE" CACHE BOOL "USE DCMTK")

add_subdirectory(${VKM_VXL_SOURCE_DIR} ${VKM_VXL_BINARY_DIR} EXCLUDE_FROM_ALL)

include_directories(${VKM_VXL_SOURCE_DIR})
include_directories(${VKM_VXL_SOURCE_DIR}/core)
include_directories(${VKM_VXL_SOURCE_DIR}/vcl)
include_directories(${VKM_VXL_SOURCE_DIR}/contrib/brl)
include_directories(${VKM_VXL_SOURCE_DIR}/contrib/brl/bseg)
include_directories(${VKM_VXL_SOURCE_DIR}/contrib/brl/bbas)
include_directories(${VKM_VXL_SOURCE_DIR}/contrib/gel)

include_directories(${VKM_VXL_BINARY_DIR}/core)
include_directories(${VKM_VXL_BINARY_DIR}/vcl)
include_directories(${VKM_VXL_BINARY_DIR}/contrib)


# Source ---------------------------------
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/vkm ${CMAKE_BINARY_DIR}/vkm-build)
include_directories(${CMAKE_CURRENT_LIST_DIR}/vkm)


# PYBIND ---------------------------------
find_package(PythonLibs 3 REQUIRED)
add_subdirectory(${VKM_EXTERNAL_DIR}/pybind11 ${CMAKE_BINARY_DIR}/pybind11-build)
add_subdirectory(${VKM_EXTERNAL_DIR}/pyvxl ${CMAKE_BINARY_DIR}/pyvxl-build)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/pyvkm ${CMAKE_BINARY_DIR}/pyvkm-build)