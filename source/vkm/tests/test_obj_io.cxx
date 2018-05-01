#include <iostream>
#include <testlib/testlib_test.h>
#include <limits>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>

// VXL ------------------------------------
#include <vgl/vgl_point_3d.h>
#include <vgl/vgl_pointset_3d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_box_3d.h>
// ----------------------------------------

//// OpenMesh -------------------------------
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
// ----------------------------------------

#include "../vkm_obj_io.h"

static void test_obj_io()
{
  using namespace OpenMesh;
  typedef OpenMesh::TriMesh_ArrayKernelT<> TriMesh;
  typedef OpenMesh::PolyMesh_ArrayKernelT<> PolyMesh;

  std::string test_dir ="C:/VisionSystems/CORE3D/CityEngine/siteplan/hypo_tests/";
  std::string scene_path = test_dir + "roof_footprint.obj";
  std::string junk_path = test_dir + "junk.obj";
  std::map<size_t, std::map<size_t,PolyMesh> > scene;
  vkm_obj_io::read_composite_obj_file(scene_path, scene);
  PolyMesh& mesh = scene[1][3];
  for(PolyMesh::VertexIter vit = mesh.vertices_begin();
	  vit != mesh.vertices_end(); ++vit)
	  std::cout << *vit << ' ' << mesh.point(*vit) << std::endl;
 // vkm_obj_io::write_composite_obj_file(junk_path,scene);  /// to do make real test
 // vkm_obj_io::write_composite_obj_file(junk_path, scene[1]);
}
TESTMAIN(test_obj_io);
