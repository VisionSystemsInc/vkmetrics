#ifndef vkm_ground_truth_h_
#define vkm_ground_truth_h_

#include <iostream>
#include <iosfwd>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <utility>

// VXL ------------------------------------
#include <vcl_compiler.h>
#include <vgl/vgl_box_2d.h>
#include <vgl/vgl_plane_3d.h>
#include <vgl/vgl_polygon.h>
#include <vgl/algo/vgl_h_matrix_2d.h>
#include <vil/vil_image_view.h>
// ----------------------------------------

// OpenMesh -------------------------------
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh> //need to include even if no IO
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
// ----------------------------------------

#include "vkm_containment_tree.h"


// parameters that control plane fitting
struct ground_truth_fit_params{
  ground_truth_fit_params() : tol_(0.25), min_consensus_(0.40), n_hypo_factor_(5.0) {}
  double tol_;  // surface fitting tolerance
  double min_consensus_;  // minimum fraction of points for a successful plane fit
  double n_hypo_factor_;  // multiplier on n^2 to define number of trials for plane hypotheses
};

// a set of algorithms for fitting planes
// depending on the type of surface, i.e. "flat" or "sloped"
class fit_plane_3d_by_type{

public:

  fit_plane_3d_by_type(std::vector<vgl_point_3d<double> > const& pts): pts_(pts), verbose_(false){}
  bool set_verbose(bool verbose){ verbose_ = verbose;}
  bool fit(std::string surface_type, vgl_plane_3d<double>& plane);
  bool fit_flat(vgl_plane_3d<double>& plane);
  bool fit_sloped(vgl_plane_3d<double>& plane);//linear least squares algorithm
  bool fit_sloped2(vgl_plane_3d<double>& plane);//simple ransac algorthm (better)

private:
  bool verbose_;
  ground_truth_fit_params params_;
  std::vector<vgl_point_3d<double> > pts_;

};

class vkm_ground_truth{

public:

  typedef OpenMesh::TriMesh_ArrayKernelT<> TriMesh;
  typedef OpenMesh::PolyMesh_ArrayKernelT<> PolyMesh;

  // surface region labels defined by Kitware
  enum surface_t {FLAT, SLOPED, ARCHED, DOMED, MISC};

  // z_off the local vertical CS elevation - tangent plane to the Earth relative to the dsm elevation
  vkm_ground_truth(): z_off_(0.0), verbose_(false) {H_.set_identity();}

  vkm_ground_truth(double z_off):z_off_(z_off), verbose_(false) {H_.set_identity();}

  void set_verbose(bool verbose){verbose_ = verbose;}

  //: file input
  // original annotation polygons in ground truth dsm image coordinates
  bool load_ground_truth_img_regions(std::string const& path);

  // region ids with assoicated type, e.g. "flat", "sloped" ...
  bool load_surface_types(std::string const& path);

  // the ground truth DSM
  bool load_dsm_image(std::string const& path);

  // a set of boundary points for a site in DSM image coordinates
  // (used to delineate a structure for metric analysis and display)
  bool load_site_perimeter(std::string const& site_name,  std::string const& bwm_ptset_path);

  //: the set of 2-d regions formed by projecting the 3-d model regions onto the x-y plane
  //  requried by the image-based footprint methods in vkm_binary_regions
  bool write_xy_polys(std::string const& site_name, std::string const& path);

  //: IO for the ground truth model
  bool write_processed_ground_truth(std::string const& path);
  static bool load_processed_ground_truth(std::string const& path,
                                          std::map<size_t, std::map<size_t, PolyMesh> >& region_meshes,
                                          std::map<size_t, surface_t>& surface_types);

  //: directly set conversion from geotif dsm image coordinates to local CS
  void set_img_to_xy_trans(vnl_matrix_fixed<double,3,3> const& M);

  //: compute conversion from geotif dsm image coordinates to local CS
  bool compute_img_to_xy_trans(std::string const& path);

  //: snap vertices to the same average location if within tol
  void snap_image_region_vertices(double tol = 2.5);

  //: discover the containmenent relations between ground truth regions
  void process_region_containment(){
    cont_tree_.process();
  }

  //: given a mc region extract the interior 3-d points and fit a 3-d plane
  void fit_region_planes();

  //: correct the topology of the multiply-connected region due to annotation choices
  // occurs mainly when hole and outer cycle boundaries are coincident but not enclosing
  bool ensure_consistent_topology(size_t outer_index, mc_region_2d& mcr);

  //: construct 3-d polygons corresponding to the multiply connected regions
  void construct_polygon_soup();

  //: convert 3-d mc regions to meshes
  void convert_to_meshes();
  const std::map<size_t, std::map<size_t, PolyMesh> >& region_meshes() const{return region_meshes_;}

  //: extrude 3-d regions to enclosing region plane otherwise the ground plane
  //  outer one cycle boundary only - useable by all mesh viewers
  bool construct_extruded_gt_model(std::string const& name);

  //: write extruded gt model in obj format
  bool write_extruded_gt_model(std::string const& name, std::string const& path);

  //: to convert 3-d regions to extruded meshes for better viewing context. Returns the face handle of the base face
  static PolyMesh::FaceHandle extrude_base_pts(const std::vector<vgl_point_3d<double> >& base_points,
                                               vgl_plane_3d<double> ground_plane, PolyMesh& mesh);

  //debug
  // convert the raw input regions to meshes for debug purposes (verts will be wildly wrong, just dem z)
  void convert_img_regions_to_meshes();
  const std::map<size_t, std::map<size_t, PolyMesh> >& img_meshes() const{return img_meshes_;}

  //: for use by metrics code to extrude test regions
  bool write_ground_planes(std::string const& gnd_plane_path) const;

private:

  //: internal functions
  // determine if reg is entirely inside xy_poly
  bool region_contained(mc_region_3d const& reg, vgl_polygon<double> const& xy_poly);

  // surface type
  //: convert string to enum
  static surface_t string_to_type(std::string const& stype);
  //: convert enum to string
  static std::string type_to_string(surface_t stype);

  // print debug output
  bool verbose_;

  //: the tranformation from dsm image coordinates to 3-d LVCS coordinates
  vgl_point_3d<double> dsm_to_world(unsigned u_dsm, unsigned v_dsm) const;

  //: the input regions defined by annotator
  std::map<size_t, std::pair<vgl_polygon<double>, vgl_box_2d<double> >  > img_regions_;

  //: a vkm_containment_tree is class for computing the containment relation between model regions
  // e.g.  the top surface of an air conditioner on a horizontal roof of a storage
  // enclosure on a lower horizonal roof. Each containing region is multiply connected
  vkm_containment_tree cont_tree_;

  //: the surface types defined by the annotator
  std::map<size_t, surface_t> surface_types_;

  //: meshes for the raw input regions
  std::map<size_t, std::map<size_t, PolyMesh> > img_meshes_;

  //: the digital elevation model
  vil_image_view<float> dsm_;

  //: the offset between dsm elevation and LVCS elevation
  double z_off_;

  //: the mapping between dsm image coordinates and LVCS X-Y coordinates
  vgl_h_matrix_2d<double> H_;

  //: the 3d fitted plane for planar regions ("flat" , "sloped")
  std::map<size_t, vgl_plane_3d<double> > region_planes_;

  //: the multiply connected regions in 3-d LVCS coordinates
  std::map<size_t, mc_region_3d> regions_3d_;

  //: meshes for the fitted regions
  std::map<size_t, std::map<size_t, PolyMesh> > region_meshes_;

  //: site perimeter typically surrounding a single building (in local LVCS coordinates )
  std::map<std::string, vgl_polygon<double> > perimeters_;

  //: local ground plane for site derived from the DSM using the perimeter
  std::map<std::string, vgl_plane_3d<double> > local_gnd_planes_;

  //: extruded 3d regions as a PolyMesh only outer boundaries
  std::map<std::string, std::map<size_t, PolyMesh> > extruded_regions_;

};

#endif
