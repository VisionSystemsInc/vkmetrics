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
  //********************
  // TODO: unit test
  //********************
}
TESTMAIN(test_obj_io);
