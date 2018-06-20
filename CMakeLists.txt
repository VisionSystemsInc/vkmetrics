
cmake_minimum_required(VERSION 3.4)
include(ExternalProject)
project("vkm_super")
set(VERBOSE ON)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)

# print cmake VERSION
message(STATUS "CMAKE_VERSION = ${CMAKE_VERSION}**********************************")

# directories of interest
set(VKM_TOP_LEVEL_DIR ${CMAKE_CURRENT_LIST_DIR})
set(VKM_EXTERNAL_DIR ${VKM_TOP_LEVEL_DIR}/external)

set(VKM_OPENMESH_SOURCE_DIR ${VKM_EXTERNAL_DIR}/openmesh)
set(VKM_OPENMESH_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/openmesh-build)
set(VKM_OPENMESH_STAMP_DIR ${CMAKE_CURRENT_BINARY_DIR}/openmesh-stamp)

set(VKM_VKM_SOURCE_DIR ${VKM_TOP_LEVEL_DIR}/source)
set(VKM_VKM_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/source-build)
set(VKM_VKM_STAMP_DIR ${CMAKE_CURRENT_BINARY_DIR}/source-stamp)

# add python site package directory
set(PYTHON_SITE /install/lib/python3.5/site-packages CACHE STRING "Python installation directory")

# default cmake arguments
list(APPEND CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
                       "-DCMAKE_CXX_FLAGS=-fdiagnostics-color=always"
                       "-DCMAKE_C_FLAGS=-fdiagnostics-color=always"
                       "-DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}"
                       #-DCMAKE_DEBUG_POSTFIX="DAN"
    )

# Building openmesh as an external project
ExternalProject_Add(openmesh
        SOURCE_DIR ${VKM_OPENMESH_SOURCE_DIR}
        CMAKE_ARGS "${CMAKE_ARGS}"
        CMAKE_ARGS -DOPENMESH_BUILD_PYTHON_BINDINGS:BOOL=0
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        BINARY_DIR ${VKM_OPENMESH_BINARY_DIR}
        STAMP_DIR ${VKM_OPENMESH_STAMP_DIR}
        BUILD_ALWAYS 1
        USES_TERMINAL_DOWNLOAD 1
        USES_TERMINAL_UPDATE 1
        USES_TERMINAL_CONFIGURE 1
        USES_TERMINAL_BUILD 1
        USES_TERMINAL_TEST 1
        USES_TERMINAL_INSTALL 1
)


# Build everything else as external project depending on OpenMesh & VXL
list(APPEND CMAKE_ARGS "-DVKM_EXTERNAL_DIR=${VKM_EXTERNAL_DIR}"
                       "-DVKM_OPENMESH_SOURCE_DIR=${VKM_OPENMESH_SOURCE_DIR}"
                       "-DVKM_OPENMESH_BINARY_DIR=${VKM_OPENMESH_BINARY_DIR}"
                       "-DPYTHON_SITE=${PYTHON_SITE}"
    )

ExternalProject_add(vkm_main
        DEPENDS openmesh
        SOURCE_DIR ${VKM_VKM_SOURCE_DIR}
        CMAKE_ARGS "${CMAKE_ARGS}"
        CMAKE_ARGS -DVXL_USE_DCMTK=0
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}
        BINARY_DIR ${VKM_VKM_BINARY_DIR}
        STAMP_DIR ${VKM_VKM_STAMP_DIR}
        USES_TERMINAL_DOWNLOAD 1
        USES_TERMINAL_UPDATE 1
        USES_TERMINAL_CONFIGURE 1
        USES_TERMINAL_BUILD 1
        USES_TERMINAL_TEST 1
        USES_TERMINAL_INSTALL 1
)

