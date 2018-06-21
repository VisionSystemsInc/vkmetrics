project("vkm")

# ----------------------------------------
# OpenMesh requirement
if (WIN32)
  add_definitions(
      -D_USE_MATH_DEFINES -DNOMINMAX
      -D_CRT_SECURE_NO_WARNINGS
  )
endif ()
# ----------------------------------------

set(VKM_SOURCES
  vkm_histogram.h
  vkm_containment_tree.h  vkm_containment_tree.cxx
  vkm_binary_regions.h  vkm_binary_regions.cxx
  vkm_obj_io.h
  vkm_ground_truth.h  vkm_ground_truth.cxx
  vkm_metrics.h  vkm_metrics.cxx
)

add_library(vkm ${VKM_SOURCES} )
target_link_libraries(vkm OpenMeshCore OpenMeshTools vgl vgl_algo vnl vil vil_algo vcl netlib)

add_subdirectory(tests)

