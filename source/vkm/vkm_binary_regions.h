#ifndef vkm_binary_regions_h
#define vkm_binary_regions_h

//:
// \file
// \brief convert a set of vgl_polygons, potentially with holes, to a binary image
// \author J.L. Mundy
// \date  March 18, 2018
//
// \verbatim
//  Modifications
//   <none yet>
// \endverbatim

#include <map>
#include <vector>
#include <string>


// VXL ------------------------------------
#include <vil/vil_image_view.h>
#include <vil/algo/vil_blob.h>
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_vector_2d.h>
#include <vgl/vgl_box_2d.h>
// ----------------------------------------


class vkm_binary_regions{

public:

  //default constructor
  vkm_binary_regions():verbose_(false){}
  bool create_binary_image(std::string const& name, std::vector<vgl_polygon<double> > const& xy_polys);
  vil_image_view<unsigned char> bin_image(std::string const& name)  {return bin_images_[name];}

  //: not used in current implmentation but can provide a speed up if needed
  // region boundaries are shown in red
  void  extract_region_boundaries(std::string const& name);
  vil_image_view<unsigned char> region_image(std::string const& name) {return region_images_[name];}

  //: compute the intersection over union if two binary images,
  // potentially with a translation of b with respect to a
  void region_metrics(std::string const& name_a, std::string const& name_b, double& int_over_union, double& completness, double& correctness,
                      double b_rel_to_a_x = 0.0, double b_rel_to_a_y = 0.0);

  //: given a translation of b with respect to a compute a color image showing the
  //  overlap of a and b (a in red channel, b in green channel)
  vil_image_view<unsigned char> generate_overlay_image(std::string const& name_a, std::string const& name_b,
                                                       double b_rel_to_a_x = 0.0, double b_rel_to_a_y = 0.0);

private:

  // internal methods
  //: map Cartesian coordinates to image coordinates for a given region indicated by "name"
  void img_coords(std::string const& name, vgl_point_2d<double> const& p2d, size_t& u, size_t& v) {
    size_t ni = bin_images_[name].ni(), nj = bin_images_[name].nj();
    const vgl_vector_2d<double>& origin = origins_[name];
    vgl_point_2d<double> dp = p2d -origin ;
    double dx = dp.x(); if(dx<0.0) dx = 0.0;
    double dy = nj-(dp.y()+1);
    if(dy<0) dy = 0.0;
    u = static_cast<size_t>(dx); v = static_cast<size_t>(dy);
    if(u>=ni) u = ni-1; if(v>=nj) v = nj-1;
  }

  //: map image coordinates to Cartesian coordinates for a given region indicated by "name"
  void xy_coords(std::string const& name, size_t u, size_t v, vgl_point_2d<double>& p2d) {
    const vgl_vector_2d<double>& origin = origins_[name];
    double ud = static_cast<double>(u), vd = static_cast<double>(bin_images_[name].nj()-(v+1));
    p2d.set(origin.x() + ud, origin.y() + vd);
  }

  bool verbose_;
  std::map<std::string, vgl_vector_2d<double> > origins_;              // origin of a given region
  std::map<std::string, vgl_box_2d<double> > bounding_boxes_;          // bounding box for the region
  std::map<std::string, vil_image_view<unsigned char> > bin_images_;   // binary image generated for the region
  std::map<std::string, vil_image_view<unsigned char> > region_images_;// color image showing region boundaries
  std::map<std::string, std::vector<vil_blob_region> > regions_;       // vil regions for the binary image

};

#endif
