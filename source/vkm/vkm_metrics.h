#ifndef vkm_metrics_h_
#define vkm_metrics_h_
//:
// \file
// \brief Compute metrics using the ground truth model
// \author J.L. Mundy
// \date  February 8, 2018
//
// \verbatim
//  Modifications
//   <none yet>
// \endverbatim

#include <iostream>
#include <iosfwd>
#include <string>
#include <vector>
#include <map>

// VXL ------------------------------------
#include <vgl/vgl_box_2d.h>
#include <vgl/vgl_plane_3d.h>
#include <vgl/vgl_polygon.h>
#include <vil/vil_image_view.h>
//-----------------------------------------

// -------------------- OpenMesh--------
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh> //need to include even if no IO
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
//--------------------------------------

#include "vkm_obj_io.h"


// a useful data structure for managing 2-d and 3-d metric computations
// manages both the 2-d region boundary in world x-y coordinates
// and a plane that is fit to the vertices of the incoming mesh
struct xy_region{
  typedef OpenMesh::PolyMesh_ArrayKernelT<> PolyMesh;
  xy_region(){}
  // construct from a mesh - a plane is fit to the vertices
  xy_region(const PolyMesh& mesh);

  // set of 3d vertices for the region sheet 0 only (no holes)
  std::vector<vgl_point_3d<double> > verts_3d() const;

  // convert back to a mesh for display purposes
  PolyMesh mesh();

  // print the vertices of the polygon for display purposes
  static void print_poly(vgl_polygon<double> const& poly);
  //
  void print_region() const{
    print_poly(region_);
  }
  //: translate the region to align the test model with the ground truth model
  xy_region translate(double tx, double ty, double tz = 0.0) const;

  //: the centroid of the region
  vgl_point_2d<double> xy_centroid() const;

  //: z translation from plane 0 to plane 1 at x-y point p
  static double z_off(vgl_plane_3d<double> const& pl0, vgl_plane_3d<double> const& pl1, vgl_point_2d<double> const& p);

  vgl_polygon<double> region_;
  vgl_box_2d<double> bounding_box_;
  vgl_plane_3d<double> plane_;
};


//: manage the metric values (scores) for matching test and ground truth regions
struct score{
  score(): gt_id_(0), gt_name_(""), test_id_(0), test_name_(""),
           gt_area_(0.0), test_area_(0.0), comp_(0.0), corr_(0.0), iou_(0.0),
           normal_ang_diff_(0.0), z_error_(0.0) {}

  //: encode json string for a score instance
  std::string json() const;

  //: members
  size_t gt_id_;            // ground truth region id
  std::string gt_name_;     // ground truth region name
  size_t test_id_;          // test region id
  std::string test_name_;   // test region name
  double gt_area_;          // ground truth region area
  double test_area_;        // test region area
  double comp_;             // completeness of the matching test region
  double corr_;             // correctness of the matching test region
  double iou_;              // intersection over union of the gt and test regions
  double normal_ang_diff_;  // angle between gt plane normal and test plane normal
  double z_error_;          // signed elevation difference at the centroid of the region intersection
};

// various operators for maintaining scores
score operator +(score const& sa, score const& sb);
score operator /=(score& a, double v);
std::ostream& operator<<(std::ostream& str, score const& s);


// parameters required in computing metrics
struct metrics_params{
  metrics_params(): min_iou_(0.05),tol_(0.5), search_radius_(5.0) {}
  double min_iou_;                  // the minimum IoU such that a gt region and test region are matched
  double tol_;                      // the spatial tolerance for geometric tests
  double search_radius_;            // radius for x y translation search
};


// main metrics class
class vkm_metrics
{
public:
  typedef OpenMesh::TriMesh_ArrayKernelT<> TriMesh;
  typedef OpenMesh::PolyMesh_ArrayKernelT<> PolyMesh;

  vkm_metrics(): z_off_(0.0), tx_(0.0), ty_(0.0), verbose_(false){}

  void set_verbose(bool verbose){verbose_ = verbose;}
  void set_translation(double tx, double ty){tx_ = tx; ty_=ty;}

  //: file input for ground truth models (covers entire area, e.g. D2)
  bool load_ground_truth_model(std::string const& path);

  //: load ground plane information
  bool load_ground_planes(std::string const& path);

  //: read test model in .OBJ format, "site" is typically a single building
  // TODO - handle multiply connected test models (requires IFC or advanced City GML format)
  bool load_simply_connected_test_model(std::string const& path, std::string site_name = "");

  // delete isolated vertices from ground truth & test model
  void delete_isolated_vertices();

  //: convert meshes to xy regions - the projection of 3-d regions onto the x-y plane
  void construct_xy_regions();

  //: find corresponding test regions with gt regions based on IoU scores
  void match_xy_regions();

  //: translate the test model using current value of tx_ and ty_
  void translate_test_model_xy();

  //: find the z translation that aligns the test region 3-d planes with gt planes
  void find_z_offset();

  //: compute 2-d scores for region pairs with maximum IoU
  void compute_best_match_2d_score_stats();

  //: save the registered test regions for dispay purposes
  bool save_transformed_regions_as_meshes(std::string const& path, std::string site_name) const;

  //: compute 3-d scores for region pairs with maximum IoU
  void compute_best_match_3d_score_stats();

  //:: output json formatted score array
  std::string json() const;

  //: print the scores in tablular form
  void print_score_array();

private:

  bool verbose_;// print debug output
  metrics_params params_;
  std::string site_name_;
  double z_off_; // shift to align test with gt in z
  double tx_;         // alignment x translation
  double ty_;         // alignment y translation
  score avg_best_match_score_; // average score for the highest IoU match to a given gt region
  score footprint_score_;      // the score for the registered footprints

  //: meshes for the ground truth region
  std::map<std::string, PolyMesh> gt_model_;

  //: local ground plane for site derived from the DEM using the perimeter
  //   used for extruding test models
  std::map<std::string, vgl_plane_3d<double> > local_gnd_planes_;

  //: multiply connected regions as polygons in x-y world coordinates for the ground truth
  std::map<size_t, std::string> gt_names_;
  std::map<size_t, xy_region> gt_regions_;

  //: meshes for test model
  std::map<std::string, PolyMesh> test_model_;

  //: multiply connected regions as polygons in x-y world coordinates for the test model
  //  currently test model region are simply connected (no holes)
  std::map<size_t, std::string> test_names_;
  std::map<size_t, xy_region> test_regions_;
  std::map<size_t, xy_region> registered_test_regions_;

  //: the matches for a gt region to a set of test regions
  //    gt region id    test region id
  //         |               |
  std::map<size_t, std::map<size_t, score> >  gt_to_test_;

  //: the matches for a test region to a set of gt_regions
  //    test region id    gt region id
  //         |               |
  std::map<size_t, std::map<size_t, score> >  test_to_gt_;

  //: the best matches with compute scores
  //      gt region id     test region id
  //         |                 |
  std::map<size_t, std::pair<size_t, score> > max_scores_;

};

#endif
