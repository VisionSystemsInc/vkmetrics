#include <iostream>
#include <testlib/testlib_test.h>
#include <limits>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>

// VXL ------------------------------------
#include <vgl/vgl_box_2d.h>
// ----------------------------------------

// OpenMesh -------------------------------
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
// ----------------------------------------

#include "../vkm_ground_truth.h"


static void test_ground_truth()
{
  typedef OpenMesh::PolyMesh_ArrayKernelT<> PolyMesh;
  std::string site_dir ="D:/core3d_experiments/phase1a_D2/GTruth_data_v2/";
  std::string ground_truth_img_region_path = site_dir + "CORE3D-Phase1a.kw18.regions";
  std::string ground_truth_pc_region_path = site_dir + "CORE3D-Phase1a.kw18.regions.xyz";
  std::string ground_truth_dem_path = site_dir + "LiDAR.tif";
  std::string ground_truth_types_path = site_dir + "CORE3D-Phase1a.kw18.types";
  std::string ground_truth_surf_path = site_dir + "CORE3D-Phase1a.kw18_surfaces.obj";
  std::string ground_truth_region_path = site_dir + "CORE3D-Phase1a.kw18_all_regions.obj";

  std::map<std::string, std::string> perimeter_paths;
  perimeter_paths["yasu_house"] = site_dir + "yasu_house_perimeter.txt";
  perimeter_paths["firehouse"] = site_dir + "firehouse_perimeter.txt";
  perimeter_paths["square_complex"] = site_dir + "square_complex_perimeter.txt";
  perimeter_paths["small_square_complex"] = site_dir + "small_square_complex_perimeter.txt";
  perimeter_paths["directorate"] = site_dir + "directorate_perimeter.txt";
  perimeter_paths["domed_bldg"] = site_dir + "domed_bldg_perimeter.txt";
  perimeter_paths["multi_storey"] = site_dir + "multi_storey_perimeter.txt";
  perimeter_paths["flat_roof_complex"] = site_dir + "flat_roof_complex_perimeter.txt";

  std::map<std::string, std::string> xy_poly_paths;
  xy_poly_paths["yasu_house"] = site_dir + "yasu_house_xy_polys.txt";
  xy_poly_paths["firehouse"] = site_dir + "firehouse_xy_polys.txt";
  xy_poly_paths["square_complex"] = site_dir + "square_complex_xy_polys.txt";
  xy_poly_paths["small_square_complex"] = site_dir + "small_square_complex_xy_polys.txt";
  xy_poly_paths["directorate"] = site_dir + "directorate_xy_polys.txt";
  xy_poly_paths["domed_bldg"] = site_dir + "domed_bldg_xy_polys.txt";
  xy_poly_paths["multi_storey"] = site_dir + "multi_storey_xy_polys.txt";
  xy_poly_paths["flat_roof_complex"] = site_dir + "flat_roof_complex_xy_polys.txt";

  std::map<std::string, std::string> ifc_paths;
  ifc_paths["yasu_house"] = site_dir + "yasu_house.ifc";
  ifc_paths["firehouse"] = site_dir + "firehouse.ifc";
  ifc_paths["square_complex"] = site_dir + "square_complex.ifc";
  ifc_paths["small_square_complex"] = site_dir + "small_square_complex.ifc";
  ifc_paths["directorate"] = site_dir + "directorate.ifc";
  ifc_paths["domed_bldg"] = site_dir + "domed_bldg.ifc";
  ifc_paths["multi_storey"] = site_dir + "multi_storey.ifc";
  ifc_paths["flat_roof_complex"] = site_dir + "flat_roof_complex.ifc";

  std::map<std::string, std::string> extruded_obj_paths;
  extruded_obj_paths["yasu_house"] = site_dir + "extruded_yasu_house.obj";
  extruded_obj_paths["firehouse"] = site_dir + "extruded_firehouse.obj";
  extruded_obj_paths["square_complex"] = site_dir + "extruded_square_complex.obj";
  extruded_obj_paths["small_square_complex"] = site_dir + "extruded_small_square_complex.obj";
  extruded_obj_paths["directorate"] = site_dir + "extruded_directorate.obj";
  extruded_obj_paths["domed_bldg"] = site_dir + "extruded_domed_bldg.obj";
  extruded_obj_paths["multi_storey"] = site_dir + "extruded_multi_storey.obj";
  extruded_obj_paths["flat_roof_complex"] = site_dir + "extruded_flat_roof_complex.obj";

  double z_off = 232.111831665;
  std::string ground_plane_path = site_dir + "site_ground_planes.txt";

  vkm_ground_truth gt(z_off);
  bool good = gt.load_ground_truth_img_regions(ground_truth_img_region_path);
   good = good && gt.load_ground_truth_pc_regions(ground_truth_pc_region_path);
   good = good && gt.compute_img_to_xy_trans();
   good = good && gt.load_dem_image(ground_truth_dem_path);
   good = good && gt.load_surface_types(ground_truth_types_path);
   good = good && gt.load_site_perimeter("square_complex", perimeter_paths["square_complex"]);
   good = good && gt.load_site_perimeter("small_square_complex", perimeter_paths["small_square_complex"]);
   good = good && gt.load_site_perimeter("yasu_house", perimeter_paths["yasu_house"]);
   good = good && gt.load_site_perimeter("firehouse", perimeter_paths["firehouse"]);
   good = good && gt.load_site_perimeter("directorate", perimeter_paths["directorate"]);
   good = good && gt.load_site_perimeter("multi_storey", perimeter_paths["multi_storey"]);
   good = good && gt.load_site_perimeter("flat_roof_complex", perimeter_paths["flat_roof_complex"]);
   good = good && gt.load_site_perimeter("domed_bldg", perimeter_paths["domed_bldg"]);
   if (good) {
	   gt.snap_image_region_vertices();
	   gt.convert_img_regions_to_meshes();
	   gt.process_region_containment();
	   gt.fit_region_planes();
	   gt.construct_polygon_soup();
	   gt.convert_to_meshes();
	   gt.write_xy_polys("square_complex", xy_poly_paths["square_complex"]);
	   gt.write_xy_polys("small_square_complex", xy_poly_paths["small_square_complex"]);
	   gt.write_xy_polys("yasu_house", xy_poly_paths["yasu_house"]);
	   gt.write_xy_polys("firehouse", xy_poly_paths["firehouse"]);
	   gt.write_xy_polys("directorate", xy_poly_paths["directorate"]);
	   gt.write_xy_polys("multi_storey", xy_poly_paths["multi_storey"]);
	   gt.write_xy_polys("flat_roof_complex", xy_poly_paths["flat_roof_complex"]);
	   gt.write_xy_polys("domed_bldg", xy_poly_paths["domed_bldg"]);
	   if (!gt.write_processed_ground_truth(ground_truth_surf_path)) {
		   std::cout << "can't write - " << ground_truth_surf_path << std::endl;
	   }
	   gt.write_ground_planes(ground_plane_path);
	   gt.construct_extruded_gt_model("square_complex");
	   gt.construct_extruded_gt_model("small_square_complex");
	   gt.construct_extruded_gt_model("yasu_house");
	   gt.construct_extruded_gt_model("firehouse");
	   gt.construct_extruded_gt_model("directorate");
	   gt.construct_extruded_gt_model("domed_bldg");
	   gt.construct_extruded_gt_model("multi_storey");
	   gt.construct_extruded_gt_model("flat_roof_complex");
	   gt.write_extruded_gt_model("square_complex", extruded_obj_paths["square_complex"]);
	   gt.write_extruded_gt_model("small_square_complex", extruded_obj_paths["small_square_complex"]);
	   gt.write_extruded_gt_model("yasu_house", extruded_obj_paths["yasu_house"]);
	   gt.write_extruded_gt_model("firehouse", extruded_obj_paths["firehouse"]);
	   gt.write_extruded_gt_model("directorate", extruded_obj_paths["directorate"]);
	   gt.write_extruded_gt_model("domed_bldg", extruded_obj_paths["domed_bldg"]);
	   gt.write_extruded_gt_model("flat_roof_complex", extruded_obj_paths["flat_roof_complex"]);
	   gt.write_extruded_gt_model("multi_storey", extruded_obj_paths["multi_storey"]);
   }
  }
TESTMAIN(test_ground_truth);
