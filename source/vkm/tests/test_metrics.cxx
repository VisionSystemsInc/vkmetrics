#include <iostream>
#include <testlib/testlib_test.h>
#include <limits>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>

// OpenMesh -------------------------------
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
// ----------------------------------------

#include "../vkm_metrics.h"

static void test_metrics()
{
  //std::string site = "yasu_house";
  //std::string site = "firehouse";
  //std::string site = "directorate";
  //std::string site = "square_complex";
  //std::string site = "small_square_complex";
  //std::string site = "multi_storey";
  std::string site = "flat_roof_complex";

  typedef OpenMesh::PolyMesh_ArrayKernelT<> PolyMesh;
  std::string gt_dir ="D:/core3d_experiments/phase1a_D2/GTruth_data_v2/";
  std::string model_dir ="D:/core3d_experiments/phase1a_D2/";
  std::map<std::string, std::string> test_model_path, registered_model_path;
  test_model_path["yasu_house"] = model_dir + "yasu_house/alpha_meshes/alpha_hull_planar_meshes.obj";
  test_model_path["firehouse"] = model_dir + "firehouse/alpha_meshes/alpha_hull_planar_meshes.obj";
  test_model_path["directorate"] = model_dir + "directorate/alpha_meshes/alpha_hull_planar_meshes.obj";
  test_model_path["square_complex"] = model_dir + "square_complex/alpha_meshes/alpha_hull_planar_meshes.obj";
  test_model_path["small_square_complex"] = model_dir + "small_square_complex/alpha_meshes/alpha_hull_planar_meshes.obj";
  test_model_path["multi_storey"] = model_dir + "multi_storey/alpha_meshes/alpha_hull_planar_meshes.obj";
  test_model_path["flat_roof_complex"] = model_dir + "flat_roof_complex/alpha_meshes/alpha_hull_planar_meshes.obj";

  registered_model_path["yasu_house"] = model_dir + "yasu_house/registered_model.obj";
  registered_model_path["firehouse"] = model_dir + "firehouse/registered_model.obj";
  registered_model_path["directorate"] = model_dir + "directorate/registered_model.obj";
  registered_model_path["square_complex"] = model_dir + "square_complex/registered_model.obj";
  registered_model_path["small_square_complex"] = model_dir + "small_square_complex/registered_model.obj";
  registered_model_path["multi_storey"] = model_dir + "multi_storey/registered_model.obj";
  registered_model_path["flat_roof_complex"] = model_dir + "flat_roof_complex/registered_model.obj";

  std::string gt_img_region_path = gt_dir + "CORE3D-Phase1a.kw18.regions";
  std::string gt_pc_region_path = gt_dir + "CORE3D-Phase1a.kw18.regions.xyz";
  std::string gt_dem_path = gt_dir + "LiDAR.tif";
  std::string gt_types_path = gt_dir + "CORE3D-Phase1a.kw18.types";
  std::string gt_surf_path = gt_dir + "CORE3D-Phase1a.kw18_surfaces.obj";
  std::string gt_region_path = gt_dir + "CORE3D-Phase1a.kw18_all_regions.obj";
  std::string ground_plane_path = gt_dir + "site_ground_planes.txt";
  // translation found by vkm_binary_regions, align binary images
  // see tests/test_binary_regions
  std::map<std::string, std::pair<double, double> > tr;
  // more robust that intersection/union of footprint polygons
  tr["square_complex"] = std::pair<double, double>(2.0, -3.0);
  tr["yasu_house"] = std::pair<double, double>(2.0, -4.0);
  tr["firehouse"] = std::pair<double, double>(3.0, 1.0);
  tr["small_square_complex"] = std::pair<double, double>(2.0, 0.0);
  tr["multi_storey"] = std::pair<double, double>(3.0, 1.0);

  // offset due to lvcs
  double z_off =   232.111831665;

  // should be replaced by modl_ground_plane
  double z_gnd_elev = 17.9843;

  vkm_metrics met(z_off, z_gnd_elev);
  met.set_translation(tr[site].first, tr[site].second);
  met.load_ground_truth_model(gt_surf_path);
  met.load_ground_planes(ground_plane_path);
  met.load_simply_connected_test_model(test_model_path[site]);
  met.construct_xy_regions();
  met.translate_test_model_xy();
  met.match_xy_regions();
  // === not a particularly useful display===
  //met.compute_match_array();
  //std::string match_array_path = model_dir + "yasu_house/yasu_house";
  //met.save_match_array(match_array_path);
  //=-================
  met.compute_best_match_2d_score_stats();
  met.find_z_offset();
  met.save_transformed_regions_as_meshes(site, registered_model_path[site]);
  met.compute_best_match_3d_score_stats();
  met.print_score_array();
  //std::cout << "======JSON ENCODING======\n"<< met.json() << std::endl;
}
TESTMAIN(test_metrics);
