add_executable( vkm_test_all
  test_driver.cxx
  test_binary_regions.cxx
  test_containment_tree.cxx
  test_obj_io.cxx
  test_ground_truth.cxx
  test_metrics.cxx
)

target_link_libraries( vkm_test_all vkm OpenMeshCore OpenMeshTools vgl vgl_algo vnl vil vil_algo vcl netlib testlib )
add_test( vkm_test_binary_regions COMMAND $<TARGET_FILE:vkm_test_all> test_binary_regions)
add_test( vkm_test_containment_tree COMMAND $<TARGET_FILE:vkm_test_all> test_containment_tree)
add_test( vkm_test_obj_io COMMAND $<TARGET_FILE:vkm_test_all> test_obj_io)
#add_test( vkm_test_ground_truth COMMAND $<TARGET_FILE:vkm_test_all> test_ground_truth)
#add_test( vkm_test_metrics COMMAND $<TARGET_FILE:vkm_test_all> test_metrics)

install(TARGETS vkm_test_all DESTINATION bin)
