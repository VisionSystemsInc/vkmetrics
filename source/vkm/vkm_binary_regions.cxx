// VXL ------------------------------------
#include <vgl/vgl_polygon_scan_iterator.h>
#include <vgl/vgl_intersection.h>
#include <vil/algo/vil_blob.h>
#include <vil/algo/vil_find_4con_boundary.h>
// ----------------------------------------

#include "vkm_binary_regions.h"


bool vkm_binary_regions::create_binary_image(std::string const& name, std::vector<vgl_polygon<double> > const& xy_polys){
  size_t n = xy_polys.size();
  if( n == 0){
    std::cout << "No xy_polys to create binary image" << std::endl;
    return false;
  }
  // compute bounding box
  vgl_box_2d<double> bb;
  for(size_t i = 0; i<n; ++i){
    for(std::vector<vgl_point_2d<double> >::const_iterator vit = xy_polys[i][0].begin();
        vit != xy_polys[i][0].end(); ++vit)
      bb.add(*vit);
  }
  bb.scale_about_centroid(1.1);
  bounding_boxes_[name] = bb;
  vgl_vector_2d<double> origin(bb.min_x(), bb.min_y());
  origins_[name] = origin;
  double w = bb.width(), h = bb.height();
  size_t ni = static_cast<size_t>(w), nj = static_cast<size_t>(h);
  ni++; nj++; // border to account for round off
  vil_image_view<unsigned char> bin_image(ni, nj);
  bin_image.fill((unsigned char) 0);
  bin_images_[name] = bin_image;
  std::vector<vgl_polygon<double> > img_polys;
  for(size_t i = 0; i<n; ++i){
    size_t ns = xy_polys[i].num_sheets();
    if(ns == 0)
      continue;
    vgl_polygon<double> img_poly;
    std::vector<vgl_point_2d<double> > img_outer;
    for(std::vector<vgl_point_2d<double> >::const_iterator vit = xy_polys[i][0].begin();
        vit != xy_polys[i][0].end(); ++vit){
      size_t u, v;
      img_coords(name,*vit, u, v);
      img_outer.push_back(vgl_point_2d<double>(u, v));
    }
    img_poly.push_back(img_outer);
    for(size_t j = 1; j<ns; ++j){
      std::vector<vgl_point_2d<double> > hole_verts;
      size_t nsi = xy_polys[i][j].size();
      for (size_t k = 0; k < nsi; ++k) {
        size_t u, v;
        img_coords(name, xy_polys[i][j][k], u, v);
        hole_verts.push_back(vgl_point_2d<double>(u, v));
      }
      img_poly.push_back(hole_verts);
    }
    img_polys.push_back(img_poly);
  }
  // iterate through the image polygons an set the binary image pixel values
  for(size_t i = 0; i<n; ++i){
    const vgl_polygon<double>& img_poly = img_polys[i];
    vgl_polygon_scan_iterator<double> pscan(img_poly);
    for(pscan.reset(); pscan.next();){
      int v  = pscan.scany();
      for(int u = pscan.startx(); u <= pscan.endx(); ++u){
        bin_images_[name](u, v) = (unsigned char) 255;
      }
    }
  }
  return true;
}

void vkm_binary_regions::extract_region_boundaries(std::string const& name){
  size_t ni = bin_images_[name].ni(), nj = bin_images_[name].nj();
  vil_image_view<unsigned char> region_image(ni, nj, 3);
  region_image.fill((unsigned char) 0);
  vil_image_view<bool> bin(ni, nj);
  bin.fill(false);
  for(size_t j = 0; j<nj; ++j)
    for(size_t i = 0; i<ni; ++i)
      if(bin_images_[name](i,j)>0){
        bin(i,j) = true;
        region_image(i, j, 0) = (unsigned char) 255;
        region_image(i, j, 1) = (unsigned char) 255;
        region_image(i, j, 2) = (unsigned char) 255;
      }

  vil_image_view<unsigned> labels(ni, nj);
  std::vector<vil_blob_region> regions;
  vil_blob_labels(bin, vil_blob_4_conn, labels);
  vil_blob_labels_to_regions(labels, regions);
  regions_[name]=regions;
  size_t nr = regions.size();
  for (size_t ir = 0; ir < nr; ++ir) {
    size_t nc = regions[ir].size();
    for(size_t k = 0; k<nc; ++k){
      size_t v = regions[ir][k].j;
      size_t u_start = regions[ir][k].ilo, u_end = regions[ir][k].ihi;
      region_image(u_start, v, 0) = (unsigned char) 255;
      region_image(u_start, v, 1) = (unsigned char) 0;
      region_image(u_start, v, 2) = (unsigned char) 0;

      region_image(u_end, v, 0) = (unsigned char) 255;
      region_image(u_end, v, 1) = (unsigned char) 0;
      region_image(u_end, v, 2) = (unsigned char) 0;

    }
  }
  region_images_[name] = region_image;
}
void vkm_binary_regions::region_metrics(std::string const& name_a, std::string const& name_b, double& int_over_union, double& int_over_area_a, double& int_over_area_b,
                    double b_rel_to_a_x, double b_rel_to_a_y ){
  const vgl_box_2d<double>& bb_a = bounding_boxes_[name_a];
  const vgl_box_2d<double>& bb_b = bounding_boxes_[name_b];
  // intersect the bounding boxes
  vgl_box_2d<double> int_bb = vgl_intersection(bb_a, bb_b);
  double w = int_bb.width(), h = int_bb.height();
  double xs = int_bb.min_x(), ys = int_bb.min_y();
  double xe = xs + w, ye = ys +h;
  double intersection_area = 0.0, union_area = 0.0;
  double area_a = 0.0, area_b = 0.0;
  for(double y = ys; y<=ye; ++y)
    for(double x = xs; x<xe; ++x){
      vgl_point_2d<double> pa(x, y);
      vgl_point_2d<double> pb(x + b_rel_to_a_x, y + b_rel_to_a_y);
      size_t u_a, v_a, u_b, v_b;
      img_coords(name_a, pa, u_a, v_a);
      img_coords(name_b, pb, u_b, v_b);
      double val_a = static_cast<double>(bin_images_[name_a](u_a, v_a))/255.0;
      double val_b = static_cast<double>(bin_images_[name_b](u_b, v_b))/255.0;
      area_a += val_a; area_b += val_b;
      double inter = val_a*val_b;
      double un = 0.0;
      if(val_a >0 || val_b > 0)
        un = 1.0;
      intersection_area += inter;
      union_area += un;
    }
  if(union_area == 0.0)
    int_over_union = 0.0;
  else
    int_over_union = intersection_area/union_area;

  if(area_a == 0.0)
    int_over_area_a = 0.0;
  else
    int_over_area_a = intersection_area/area_a;

  if(area_b == 0.0)
    int_over_area_b = 0.0;
  else
    int_over_area_b = intersection_area/area_b;
}

vil_image_view<unsigned char> vkm_binary_regions::generate_overlay_image(std::string const& name_a, std::string const& name_b,
                                                                     double b_rel_to_a_x , double b_rel_to_a_y){
  const vgl_box_2d<double>& bb_a = bounding_boxes_[name_a];
  const vgl_box_2d<double>& bb_b = bounding_boxes_[name_b];
  // intersect the bounding boxes
  vgl_box_2d<double> int_bb = vgl_intersection(bb_a, bb_b);
  double w = int_bb.width(), h = int_bb.height();
  double xs = int_bb.min_x(), ys = int_bb.min_y();
  double xe = xs + w, ye = ys +h;
  size_t ni = w, nj = h;
  ni++; nj++;
  vil_image_view<unsigned char> overlay(ni, nj, 3);
  size_t u = 0, v = 0;
  for (double y = ys; y <= ye; ++y, ++v) {
	  u = 0;
	  for (double x = xs; x < xe; ++x, ++u) {
		  vgl_point_2d<double> pa(x, y);
		  vgl_point_2d<double> pb(x + b_rel_to_a_x, y + b_rel_to_a_y);
		  size_t u_a, v_a, u_b, v_b;
		  img_coords(name_a, pa, u_a, v_a);
		  img_coords(name_b, pb, u_b, v_b);
		  unsigned char val_a = bin_images_[name_a](u_a, v_a);
		  unsigned char val_b = bin_images_[name_b](u_b, v_b);
		  overlay(u, v, 0) = val_a;
		  overlay(u, v, 1) = val_b;
		  overlay(u, v, 2) = (unsigned char) 0;
	  }
  }
  return overlay;
}
