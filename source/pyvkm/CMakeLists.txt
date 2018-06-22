project("pyvkm")

include_directories(${PYTHON_INCLUDE_PATH})

set(PYVKM_SOURCES
  pyvkm.cxx
  pyvkm_ground_truth.h  pyvkm_ground_truth.cxx
  pyvkm_metrics.h  pyvkm_metrics.cxx
  pyvkm_obj_io.h  pyvkm_obj_io.cxx
)

pybind11_add_module(pyvkm ${PYVKM_SOURCES})
target_link_libraries(pyvkm PRIVATE vkm OpenMeshCore OpenMeshTools)
set_target_properties(pyvkm PROPERTIES OUTPUT_NAME "vkm")

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
    from distutils import sysconfig as sc
    print(sc.get_python_lib(plat_specific=True))"
  OUTPUT_VARIABLE PYTHON_SITE_DEFAULT
  OUTPUT_STRIP_TRAILING_WHITESPACE)
set(PYTHON_SITE ${PYTHON_SITE_DEFAULT} CACHE STRING "Python installation directory")

install(TARGETS pyvkm DESTINATION ${PYTHON_SITE})
