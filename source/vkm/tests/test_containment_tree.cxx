#include <iostream>
#include <testlib/testlib_test.h>
#include <limits>
#include <string>
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>

// VXL ------------------------------------
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_vector_3d.h>
#include <vgl/vgl_pointset_3d.h>
#include <vgl/vgl_fit_oriented_box_2d.h>
#include <vgl/vgl_polygon.h>
#include <vgl/vgl_oriented_box_2d.h>
#include <vgl/vgl_clip.h>
// ----------------------------------------

#include "../vkm_containment_tree.h"


// cntx
static void test_containment_tree()
{
  vgl_point_2d<double> p00(0.0, 0.0), p01(10.0, 0.0), p02(10.0, 10.0), p03(0.0, 10.0);
  vgl_point_2d<double> p10(2.0, 2.0), p11(8.0, 2.0), p12(8.0, 8.0), p13(2.0, 8.0);
  vgl_point_2d<double> p20(5.0, 4.0), p21(6.0, 4.0), p22(6.0, 6.0), p23(5.0, 6.0);
  vgl_point_2d<double> p30(3.0, 3.0), p31(4.0, 3.0), p32(4.0, 3.5), p33(3.0, 3.5);
  vgl_point_2d<double> p40(15.0, 0.0), p41(20.0, 0.0), p42(20.0, 10.0), p43(15.0, 10.0);
  std::vector<vgl_point_2d<double> > verts0, verts1, verts2, verts3, verts4;
  verts0.push_back(p00);   verts0.push_back(p01);   verts0.push_back(p02);   verts0.push_back(p03);
  verts1.push_back(p10);   verts1.push_back(p11);   verts1.push_back(p12);   verts1.push_back(p13);
  verts2.push_back(p20);   verts2.push_back(p21);   verts2.push_back(p22);   verts2.push_back(p23);
  verts3.push_back(p30);   verts3.push_back(p31);   verts3.push_back(p32);   verts3.push_back(p33);
  verts4.push_back(p40);   verts4.push_back(p41);   verts4.push_back(p42);   verts4.push_back(p43);
  vgl_polygon<double> poly0(verts0), poly1(verts1), poly2(verts2), poly3(verts3), poly4(verts4);
  std::map<size_t, vgl_polygon<double> > poly_map;
  poly_map[0]=poly0; poly_map[1]=poly1; poly_map[2]=poly2; poly_map[3]=poly3; poly_map[4]=poly4;
  vkm_containment_tree cont_tree(poly_map);
  cont_tree.process();
  std::map < size_t,  mc_region_2d>& mcrs = cont_tree.mc_regions();
 bool good =  mcrs[0].holes_.size() == 1;
  good = good && mcrs[1].holes_.size() == 2;
  good = good && mcrs[2].holes_.size() == 0;
  good = good && mcrs[3].holes_.size() == 0;
  good = good && mcrs[4].holes_.size() == 0;
  TEST("test region containment", good, true);
}
TESTMAIN(test_containment_tree);
