#include <fstream>
#include <math.h>

// VXL ------------------------------------
#include <vgl/vgl_point_2d.h>
#include <vgl/vgl_homg_point_2d.h>
#include <vgl/vgl_homg_point_3d.h>
#include <vgl/vgl_area.h>
#include <vgl/vgl_distance.h>
#include <vgl/vgl_intersection.h>
#include <vgl/vgl_clip.h>
#include <vgl/algo/vgl_h_matrix_2d_compute_linear.h>
#include <vgl/algo/vgl_fit_plane_3d.h>
#include <vgl/vgl_polygon_scan_iterator.h>
#include <vgl/vgl_closest_point.h>
#include <vnl/vnl_random.h>
#include <vil/vil_load.h>
// ----------------------------------------

#include "vkm_histogram.h"
#include "vkm_obj_io.h"
#include "vkm_ground_truth.h"


#define test_region_indx 4


bool fit_plane_3d_by_type::fit(
    std::string surface_type, vgl_plane_3d<double>& plane)
{
  if(surface_type == "Flat")
    return fit_flat(plane);
  else if(surface_type == "Sloped")
    return fit_sloped2(plane);
  std::cout << "surface type: " << surface_type << " unknown" << std::endl;
  return false;
}


bool fit_plane_3d_by_type::fit_flat(vgl_plane_3d<double>& plane)
{
  // find range of elevations
  size_t n = pts_.size();
  if (n == 0) {
    std::cout << "no points to fit plane" << std::endl;
    return false;
  }
  double avg_z = 0.0;
  double z_min = std::numeric_limits<double>::max(), z_max = -z_min;
  for(size_t i = 0; i<n; ++i){
    double z = pts_[i].z();
    if(z < z_min) z_min = z;
    if(z > z_max) z_max = z;
  }
  if ((z_max - z_min) <= params_.tol_) { // all elevations within tolerance
    avg_z = 0.5*(z_max + z_min);
    plane.set(0.0, 0.0, 1.0, -avg_z);
    return true;
  }else{// filter outliers
    size_t nbins = static_cast<size_t>((z_max - z_min) / params_.tol_);
    vkm_histogram<double> h(nbins, z_min, z_max);
    for (size_t i = 0; i < n; ++i)
      h.upcount(pts_[i].z(), 1.0);
    //find most probable z
    double most_prob_z = h.most_probable_value();
    double avg_z = 0.0, cnt = 0.0;
    for (size_t i = 0; i < n; ++i) {
      double z = pts_[i].z();
      if (fabs(z - most_prob_z) <= params_.tol_) {
        avg_z += z;
        cnt += 1.0;
      }
    }
    if (cnt > 0.0) {
      avg_z /= cnt;
      plane.set(0.0, 0.0, 1.0, -avg_z);
      return true;
    }
  }
  std::cout << "no points left after filtering" << std::endl;
  return false;
}


bool fit_plane_3d_by_type::fit_sloped(vgl_plane_3d<double>& plane)
{
  size_t n = pts_.size();
  if (n == 0) {
	  std::cout << "no points to fit plane" << std::endl;
	  return false;
  }
  std::vector<vgl_homg_point_3d<double> > plane_pts;
  for (size_t i = 0; i < n; ++i){
    vgl_homg_point_3d<double> hp3d(pts_[i]);
    plane_pts.push_back(hp3d);
  }
  // now fit the plane including outliers
  vgl_fit_plane_3d<double> fitter(plane_pts);
  double large_tol = 1000.0;
  if(!fitter.fit(large_tol, &std::cout)){
    std::cout << "initial plane fit failed" << std::endl;
    return false;
  }
  vgl_homg_plane_3d<double> plf = fitter.get_plane();
  plf.normalize();
  double avg_dist = 0.0, nd = 0.0;
  std::vector<double> dists;
  for (size_t i = 0; i < plane_pts.size(); ++i, nd+=1.0) {
    double alg_dist = plane_pts[i].x()*plf.a() + plane_pts[i].y()*plf.b() + plane_pts[i].z()*plf.c() + plane_pts[i].w()*plf.d();
    dists.push_back( fabs(alg_dist));
    avg_dist += fabs(alg_dist);
  }
  avg_dist /= nd;
  std::cout << "Avg error " << avg_dist << std::endl;
  std::vector<vgl_homg_point_3d<double> > temp;
  for (size_t i = 0; i < plane_pts.size(); ++i) {
    if (dists[i] > avg_dist)
      continue;
    temp.push_back(plane_pts[i]);
  }
  if (temp.size() == 0) {
    std::cout << "plane_fit failed  - no points after filter " << std::endl;
    return false;
  }
  // now fit the plane without outliers
  vgl_fit_plane_3d<double> fitter2(temp);
  double small_tol = 4.0*avg_dist;
  if (!fitter2.fit(small_tol, &std::cout)) {
    std::cout << "plane_fit failed after removal of outliers " << std::endl;
    return false;
  }
  vgl_homg_plane_3d<double> plf2 = fitter2.get_plane();
  plf2.normalize();
  plane.set(plf2.a(), plf2.b(), plf2.c(), plf2.d());
  return true;
}


bool fit_plane_3d_by_type::fit_sloped2(vgl_plane_3d<double>& plane)
{
  size_t n = pts_.size();
  if (n <3 ) {
    std::cout << "insufficient points to fit plane" << std::endl;
    return false;
  }
  double max_support = 0.0;
  vgl_plane_3d<double> max_pl;
  vnl_random rand;
  // generate k x n^2 random triples
  size_t n_hypos = params_.n_hypo_factor_*n*n;
  for(size_t i = 0; i<n_hypos; ++i){
  size_t i0 = rand(n);
  size_t i1 = rand(n);
  while(i1 == i0)
    i1 = rand(n);
  size_t i2 = rand(n);
  while(i2 == i1 || i2 == i0)
    i2 = rand(n);
  const vgl_point_3d<double>& p0 = pts_[i0];
  const vgl_point_3d<double>& p1 = pts_[i1];
  const vgl_point_3d<double>& p2 = pts_[i2];
  vgl_vector_3d<double> v01 = p1-p0;
  vgl_vector_3d<double> v02 = p2-p0;
  // want nearby points
  if (v01.length() > 1.0 || v02.length() > 1.0)
    continue;
  v01/=v01.length(); v02/=v02.length();
  vgl_vector_3d<double> normal = cross_product(v01, v02);
  // don't want collinear points
  if(normal.length() < 0.025)
    continue;
  normal/=(normal.length());
  vgl_plane_3d<double> pl(normal, p0);
  double sup = 0.0;
  for(size_t k = 0; k<n; ++k)
    if(vgl_distance(pl,pts_[k])<params_.tol_)
      sup += 1.0;
  if(sup>max_support){
    max_support = sup;
    max_pl = pl;
  }
  }
  if(verbose_) std::cout << "npts " << n << " % sup " << 100*max_support/n << std::endl;
  if((max_support/n) < params_.min_consensus_){
    std::cout << "insufficient support for plane fit " << 100*max_support/n << " %" << std::endl;
    return false;
  }
  plane.set(max_pl.a(), max_pl.b(), max_pl.c(), max_pl.d());
  return true;
}


bool vkm_ground_truth::load_ground_truth_img_regions(std::string const& path)
{
  std::ifstream istr(path.c_str());
  if(!istr){
    std::cout << "Can't open " << path << " to read ground truth regions" << std::endl;
    return false;
  }

  size_t indx, key0, key1, nverts;
  while(istr >> indx >> key0 >> key1 >> nverts){
    if(nverts == 0)
      continue;
    double x= 0.0, y = 0.0;
    std::vector<vgl_point_2d<double> > verts;
    vgl_box_2d<double> bb;
    for(size_t i = 0; i<nverts; ++i){
      istr >> x >> y;
      vgl_point_2d<double> p(x, y);
      bb.add(p);
      verts.push_back(p);
    }
    istr >> std::ws;
    if(vgl_area(bb) == 0){
      std::cout << "zero area region[" << indx << "] - fatal " << std::endl;
      return false;
    }
    std::pair<vgl_polygon<double>, vgl_box_2d<double> > pr(vgl_polygon<double>(verts), bb);
    img_regions_[indx] = pr;
  }
  cont_tree_ = vkm_containment_tree(img_regions_);
  return true;
}


// find close vertex pairs and replace with average
void vkm_ground_truth::snap_image_region_vertices(double tol)
{
  for(std::map<size_t, std::pair<vgl_polygon<double>, vgl_box_2d<double> >  >::iterator riti = img_regions_.begin();
      riti != img_regions_.end();++ riti){
    size_t indx_i = riti->first;
    std::vector<vgl_point_2d<double> >& verts_i = riti->second.first[0];
    size_t ni = verts_i.size();
    for(std::map<size_t, std::pair<vgl_polygon<double>, vgl_box_2d<double> >  >::iterator ritj = img_regions_.begin();
      ritj != img_regions_.end();++ ritj){
      size_t indx_j = ritj->first;
      if(indx_j <= indx_i)
        continue;
      std::vector<vgl_point_2d<double> >& verts_j = ritj->second.first[0];
      size_t nj = verts_j.size();
      for(size_t i = 0; i<ni; ++i){
        vgl_point_2d<double>& pi = verts_i[i];
        for(size_t j =0; j<nj; ++j){
          vgl_point_2d<double>& pj = verts_j[j];
          double d = vgl_distance(pi, pj);
          if(d<tol){ // average vertex position
            vgl_point_2d<double> vavg(0.5*(pi.x()+pj.x()),0.5*(pi.y()+pj.y()));
            pi.set(vavg.x(), vavg.y());
            pj.set(vavg.x(), vavg.y());
          }
        }
      }
    }
  }
}


vkm_ground_truth::surface_t vkm_ground_truth::string_to_type(std::string const& stype)
{
  std::map<std::string, surface_t> type_map;
  type_map["Flat"]   = FLAT;
  type_map["Sloped"] = SLOPED;
  type_map["Arched"] = ARCHED;
  type_map["Domed"]  = DOMED;
  type_map["Misc"]   = MISC;
  return type_map[stype];
}


std::string vkm_ground_truth::type_to_string(vkm_ground_truth::surface_t stype)
{
  switch(stype){
  case FLAT:
    return "Flat";
  case SLOPED:
    return "Sloped";
  case ARCHED:
    return "Arched";
  case DOMED:
    return "Domed";
  case MISC:
    return "Misc";
  default:{
    std::cout << "unknown type " << std::endl;
    return "";}
  }
}


bool vkm_ground_truth::load_surface_types(std::string const& path)
{
  std::ifstream istr(path.c_str());
  if(!istr){
    std::cout << "Can't open " << path << " to read ground truth surface types" << std::endl;
    return false;
  }
  size_t index = 0;
  std::string stype;
  while(istr >> index >> stype){
    surface_types_[index] = string_to_type(stype);
    istr >> std::ws;
  }
  return true;
}


bool vkm_ground_truth::compute_img_to_xy_trans(std::string const& path)
{
  //: region vertex positions in LVCS coordinates
  std::map<size_t, std::vector<vgl_point_2d<double > > > pc_regions;

  // locate file
  std::ifstream istr(path.c_str());
  if(!istr){
    std::cout << "Can't open " << path << " to read ground truth point cloud regions" << std::endl;
    return false;
  }

  // read file
  double x=0.0, y=0.0, dindx=0.0;
  while(istr >> x >> y >> dindx){
    size_t indx = static_cast<size_t>(dindx);
    vgl_point_2d<double> vert(x, y);
    pc_regions[indx].push_back(vert);
    istr >> std::ws;
  }

  // get corresponding points
  std::vector<vgl_homg_point_2d<double> >img_pts;
  std::vector<vgl_homg_point_2d<double> > xy_pts;
  for (std::map<size_t, std::pair<vgl_polygon<double>, vgl_box_2d<double> >  >::iterator rit = img_regions_.begin();
       rit != img_regions_.end(); ++rit) {
    size_t indx = rit->first;
    const std::vector<vgl_point_2d<double> >& iverts = (rit->second.first)[0];
    const std::vector<vgl_point_2d<double> >& xyverts = pc_regions[indx];
    size_t ni = iverts.size(), nxy = xyverts.size();
    if(ni == 0 || nxy == 0){
      std::cout << "no points available - fatal" << std::endl;
      return false;
    }
    if(ni != nxy){
      std::cout << "mismatching verts in compute trans-- fatal" << std::endl;
      return false;
    }
    for(size_t i = 0; i< ni; ++i){
      img_pts.push_back(vgl_homg_point_2d<double>(iverts[i]));
      xy_pts.push_back(vgl_homg_point_2d<double>(xyverts[i]));
    }
  }

  // compute homography
  vgl_h_matrix_2d_compute_linear h_comp;
  if(!h_comp.compute(img_pts, xy_pts, H_)){
    std::cout << " compute homography failed - fatal" << std::endl;
    return false;
  }

  // report
  std::cout << "Computed homography \n" << H_ << std::endl;
  return true;
}


void vkm_ground_truth::set_img_to_xy_trans(vnl_matrix_fixed<double,3,3> const& M)
{
  H_.set(M);
  std::cout << "Manual homography \n" << H_ << std::endl;
}


bool vkm_ground_truth::load_dem_image(std::string const& path)
{
  dem_ = vil_load(path.c_str());
  if(dem_.ni()==0){
    std::cout << "load dem failed for path " << path << std::endl;
    return false;
  }
  return true;
}


vgl_point_3d<double> vkm_ground_truth::dem_to_world(unsigned u_dem, unsigned v_dem) const
{
  double z = dem_(u_dem, v_dem)-z_off_;
  // transform to x-y
  double du = static_cast<double>(u_dem), dv = static_cast<double>(v_dem);
  vgl_homg_point_2d<double> img_pt_h(du, dv), xy_pt_h;
  xy_pt_h = H_*img_pt_h;
  vgl_point_2d<double> xy_pt(xy_pt_h);
  return vgl_point_3d<double>(xy_pt.x(), xy_pt.y(), z);
}


void vkm_ground_truth::fit_region_planes()
{
  //size_t debug_indx = 110;
  for (std::map<size_t, mc_region_2d >::const_iterator mit = cont_tree_.mc_regions().begin();
       mit != cont_tree_.mc_regions().end(); ++mit) {
    size_t indx = mit->first;
    surface_t stype = surface_types_[indx];
    if( stype == ARCHED || stype == DOMED || stype == MISC)
      continue;

    std::vector<vgl_point_3d<double> > ppts;
    //if(indx == debug_indx)
    //std::cout << "pts for index " << indx << std::endl;

    //         =========      don't include poly boundary =========>V<
    vgl_polygon_scan_iterator<double> scan_it(mit->second.poly(), false);
    for(scan_it.reset();scan_it.next();){
      int v = scan_it.scany();
      for( int u = (scan_it.startx()); u <= (scan_it.endx()); ++u){
        vgl_point_3d<double> p3d = dem_to_world(u, v);
        ppts.push_back(p3d);
        //     if(indx == debug_indx)
        //std::cout << p3d.x() << ' ' << p3d.y() << ' ' << p3d.z() << std::endl;
      }
    }
    vgl_plane_3d<double> pl;
    fit_plane_3d_by_type fpt(ppts);
    if(!fpt.fit(type_to_string(stype),  pl)){
      std::cout << "plane_fit failed for index " << indx << std::endl;
      continue;
    }
    std::cout << " region id " << indx << std::endl;
    region_planes_[indx]= pl;
  }
}


static bool proj_vertical(
    vgl_point_3d<double> const& p, vgl_plane_3d<double> const& plane,
    vgl_point_3d<double>& proj_pt)
{
  if(fabs(plane.c())<0.1){
    std::cout << "nearly vertical plane - can't project point vertically" << std::endl;
    return false;
  }
  double z = -(plane.a()*p.x() + plane.b()*p.y() + plane.d())/plane.c();
  proj_pt.set(p.x(), p.y(), z);
  return true;
}


static bool matched(
    vgl_line_segment_2d<double> const& seg,
    std::vector<vgl_line_segment_2d<double> > const& segs,
    std::vector<std::pair<size_t, size_t> > const& edge_prs,
    std::pair<size_t, size_t>& segs_match, double tol)
{
  double ang_tol = 0.01; //tol in radians
   size_t n = segs.size();
  if(n == 0)
    return false;
  // the segment, seg,  is matched if it lies on or in a segment from segs for example,
  // 0--------------0 segs[i]
  //   o------------o seg
  //

  vgl_vector_2d<double> seg_dir = seg.direction();
  const vgl_point_2d<double>& p1 = seg.point1(), p2 = seg.point2();
  for(size_t i = 0; i<n; ++i){
    // 1) are directions the same ?
    vgl_vector_2d<double> dir = segs[i].direction();
    double dp = dot_product(seg_dir, dir);
    if((1.0-fabs(dp)) > ang_tol)
      continue;
    // yes. 2) do the endpoints of seg lie on segs[i]
    double d1 = vgl_distance(p1, segs[i]);
    double d2 = vgl_distance(p2, segs[i]);
    if(d1>tol)
      continue;
    if(d2>tol)
      continue;
    segs_match = edge_prs[i];
    // yes found match
    return true;
  }
  return false;
}


static bool repair_topology_caseI(
    mc_region_2d& mcr_2d,
    std::vector<std::pair<size_t, size_t> > const& matched_edges,
    vgl_line_segment_2d<double> const& unmatched_seg,
    size_t hole_index, double tol)
{
  size_t n = matched_edges.size();
  if (n == 0) {
    std::cout << "shouldn't happen - no matched edges" << std::endl;
    return false;
  }

  // the outer boundary will be edited to include hole
  std::vector<vgl_point_2d<double> > out_verts = mcr_2d.outer_cycle_;

  //if outer boundary and hole boundary are exactly coincient can cause drawing issues
  vgl_vector_2d<double> hole_translation = unmatched_seg.normal();

  // must translate at least one dem pixel distance (discrete dem image coordinates)
  double max_component = fabs(hole_translation.x());
  if (fabs(hole_translation.y()) > max_component)
	  max_component = fabs(hole_translation.y());
  double scale = 1.1 / max_component;
  hole_translation *= scale;

  // find the sense of the perpendicular translation
  std::vector<vgl_point_2d<double> > hole_verts = mcr_2d.holes_[hole_index];
  size_t nh = hole_verts.size();
  if(!nh){
    std::cout << "hole " << hole_index << " null boundary" << std::endl;
    return false;
  }
  // compute hole centroid
  double cx = 0.0, cy = 0.0;
  for(size_t i = 0; i<nh; ++i){
    cx += hole_verts[i].x();
    cy += hole_verts[i].y();
  }
  vgl_vector_2d<double> c(cx/nh, cy/nh);

  //mid point of seg
  vgl_vector_2d<double> m(0.5*(unmatched_seg.point1().x() + unmatched_seg.point2().x()),
                          0.5*(unmatched_seg.point1().y() + unmatched_seg.point2().y()));

  // translation from seg midpoint to hole centroid defines the sense of the hole translation
  vgl_vector_2d<double> t = c-m;
  double sign = dot_product(t, hole_translation);
  if(sign < 0)
    hole_translation *= -1.0;

  // actually translate the hole
  for(size_t i = 0; i<n; ++i)
    hole_verts[i].set(hole_verts[i].x() + hole_translation.x(),
                      hole_verts[i].y() + hole_translation.y());
  mcr_2d.holes_[hole_index] = hole_verts;

  // handle 4 cases of seg intersection with the outer boundary

  // convert edge pairs to unique vert ids
  std::set<size_t> unique_vert_indices;
  for (std::vector<std::pair<size_t, size_t> >::const_iterator iit = matched_edges.begin();
       iit != matched_edges.end(); ++iit) {
    unique_vert_indices.insert(iit->first);
    unique_vert_indices.insert(iit->second);
  }

  //search for p1
  bool p1_matched = false;
  const vgl_point_2d<double>& p1 = unmatched_seg.point1();
  size_t p1_match_index = -1;
  for(std::set<size_t>::iterator vit = unique_vert_indices.begin();
      vit != unique_vert_indices.end(); ++vit){
    const vgl_point_2d<double>& p = out_verts[*vit];
    if((p1-p).length() < tol){
      p1_matched = true;
      p1_match_index = *vit;
    }
  }
  // if not matched p1 must split a matched edge
  std::pair<size_t, size_t> p1_split_edge;
  bool p1_edge_found = false;
  if(!p1_matched){
    for(std::vector<std::pair<size_t, size_t> >::const_iterator iit = matched_edges.begin();
        iit != matched_edges.end()&&!p1_edge_found; ++iit){
      const vgl_point_2d<double>& m1 = out_verts[iit->first];
      const vgl_point_2d<double>& m2 = out_verts[iit->second];
      vgl_line_segment_2d<double> s12(m1, m2);
      double d = vgl_distance(p1, s12);
      if(d<tol){
        p1_split_edge = *iit;
        p1_edge_found = true;
      }
    }
  }
  //search for p2
  bool p2_matched = false;
  const vgl_point_2d<double>& p2 = unmatched_seg.point2();
  size_t p2_match_index = -1;
  for (std::set<size_t>::iterator vit = unique_vert_indices.begin();
       vit != unique_vert_indices.end(); ++vit) {
    const vgl_point_2d<double>& p = out_verts[*vit];
    if ((p2 - p).length() < tol) {
      p2_matched = true;
      p2_match_index = *vit;
    }
  }
  // if not matched p2 must split a matched edge
  std::pair<size_t, size_t> p2_split_edge;
  bool p2_edge_found = false;
  if(!p2_matched){
    for(std::vector<std::pair<size_t, size_t> >::const_iterator iit = matched_edges.begin();
        iit != matched_edges.end()&&!p2_edge_found; ++iit){
      const vgl_point_2d<double>& m1 = out_verts[iit->first];
      const vgl_point_2d<double>& m2 = out_verts[iit->second];
      vgl_line_segment_2d<double> s12(m1, m2);
      double d = vgl_distance(p2, s12);
      if(d<tol){
        p2_split_edge = *iit;
        p2_edge_found = true;
      }
    }
  }
  std::cout << "p1 matched " << p1_matched << " p2 matched " << p2_matched << std::endl;
  if(p1_matched){
    std::cout << " p1 " << p1 ;
    if(p2_edge_found)
      std::cout << " p2 edge split (" << p2_split_edge.first << ' ' << p2_split_edge.second << ")" << std::endl;
  }
  if(p2_matched){
    std::cout << " p2 " << p2 ;
    if(p1_edge_found)
      std::cout << " p1 edge split " << p1_split_edge.first << ' ' << p1_split_edge.second << ")" << std::endl;
  }
  std::vector<vgl_point_2d<double> >  to_erase; // vertices to remove from outer cycle
  // case I p1 and p2 matched - just remove matched edges
  if(p1_matched && p2_matched){
    for (std::set<size_t>::iterator vit = unique_vert_indices.begin();
         vit != unique_vert_indices.end(); ++vit){
      if(*vit == p1_match_index || *vit == p2_match_index)
        continue;
      to_erase.push_back(out_verts[*vit]);
    }
    for(std::vector<vgl_point_2d<double> >::iterator rit = to_erase.begin();
        rit != to_erase.end(); ++rit){
      std::vector<vgl_point_2d<double> >::iterator eit = std::find(out_verts.begin(), out_verts.end(), *rit);
      out_verts.erase(eit);
    }
    mcr_2d.outer_cycle_ = out_verts;
    return true;
  }

  // case II p1 only matched  - break the edge that p2 intersects. remove matching partial edge
  if (p1_matched && !p2_matched) {
    size_t vidxf = p2_split_edge.first;
    size_t vidxs = p2_split_edge.second;

    // insert p2 vertex between first and second vertices of split edge
    size_t low_vert_idx = vidxf;
    if (vidxs < low_vert_idx)
      low_vert_idx = vidxs;
    const vgl_point_2d<double> lv = out_verts[low_vert_idx];

    // vit is the location to insert the split vertex p2
    std::vector<vgl_point_2d<double> >::iterator vit = std::find(out_verts.begin(), out_verts.end(), lv);
    if (vit == out_verts.end()) {
      std::cout << "insert p2 vert location not found -- fatal" << std::endl;
      return false;
    }

    // find which split edge vert (first or second) is inside the hole boundary
    bool first_inside = false;
    for (std::vector<std::pair<size_t, size_t> >::const_iterator iit = matched_edges.begin();
         iit != matched_edges.end() && !first_inside; ++iit) {
      if (*iit == p2_split_edge)
        continue;
      if (vidxf == iit->first || vidxf == iit->second)
        first_inside = true;
    }
    // set the vertices of the outer boundary to erase
    for (std::set<size_t>::iterator uit = unique_vert_indices.begin();
         uit != unique_vert_indices.end(); ++uit) {
      if ((first_inside && *uit == vidxs) || (!first_inside && *uit == vidxf))
        continue;
      to_erase.push_back(out_verts[*uit]);
    }
    // insert the extra p2 vert after finding verts to erase since insert changes indices
    out_verts.insert(vit, p2);

    for (std::vector<vgl_point_2d<double> >::iterator rit = to_erase.begin();
         rit != to_erase.end(); ++rit) {
      std::vector<vgl_point_2d<double> >::iterator eit = std::find(out_verts.begin(), out_verts.end(), *rit);
      out_verts.erase(eit);
    }
    mcr_2d.outer_cycle_ = out_verts;
    return true;
  }

  // case III p2 only matched  - break the edge that p1 intersects. remove matching partial edge
  if(!p1_matched && p2_matched){
    size_t vidxf = p1_split_edge.first;
    size_t vidxs = p1_split_edge.second;
    // insert p1 vertex between first and second vertices of the split edge
    size_t low_vert_idx = vidxf;
    if (vidxs < low_vert_idx)
      low_vert_idx = vidxs;
    const vgl_point_2d<double> lv = out_verts[low_vert_idx];

    // vit is the location to insert the split vertex p1
    std::vector<vgl_point_2d<double> >::iterator vit = std::find(out_verts.begin(), out_verts.end(), lv);
    if (vit == out_verts.end()) {
      std::cout << "insert p1 vert location not found -- fatal" << std::endl;
      return false;
    }

    // find which vert (first or second) is outside the hole boundary
    bool first_inside = false;
    for (std::vector<std::pair<size_t, size_t> >::const_iterator iit = matched_edges.begin();
         iit != matched_edges.end() && !first_inside; ++iit) {
      if (*iit == p1_split_edge)
        continue;
      if (vidxf == iit->first || vidxf == iit->second)
        first_inside = true;
    }
    // set the vertices of the outer boundary to erase
    for (std::set<size_t>::iterator uit = unique_vert_indices.begin();
         uit != unique_vert_indices.end(); ++uit) {
      if ((first_inside && *uit == vidxs) || (!first_inside && *uit == vidxf))
        continue;
      to_erase.push_back(out_verts[*uit]);
    }

    // insert the extra p1 vert after finding verts to erase since insert changes indices
    out_verts.insert(vit, p1);
    for (std::vector<vgl_point_2d<double> >::iterator rit = to_erase.begin();
         rit != to_erase.end(); ++rit) {
      std::vector<vgl_point_2d<double> >::iterator eit = std::find(out_verts.begin(), out_verts.end(), *rit);
      out_verts.erase(eit);
    }
    mcr_2d.outer_cycle_ = out_verts;
    return true;

  }

  // cast IV  no outer cycle vertex match - split both edges
  if(!p1_matched && !p2_matched){
    //insert p1 and p2 into outer boundary
    size_t vidxf = p1_split_edge.first;
    size_t vidxs = p1_split_edge.second;

    // insert p1 vertex between first and second
    size_t low_vert_idx = vidxf;
    if (vidxs < low_vert_idx)
      low_vert_idx = vidxs;
    const vgl_point_2d<double> lv = out_verts[low_vert_idx];

    // location to insert p1
    std::vector<vgl_point_2d<double> >::iterator vit = std::find(out_verts.begin(), out_verts.end(), lv);
    if (vit == out_verts.end()) {
      std::cout << "insert p1 vert location not found -- fatal" << std::endl;
      return false;
    }
    size_t vidxf2 = p2_split_edge.first;
    size_t vidxs2 = p2_split_edge.second;

    // insert p2 vertex between first and second
    size_t low_vert_idx2 = vidxf2;
    if (vidxs2 < low_vert_idx2)
      low_vert_idx2 = vidxs2;
    const vgl_point_2d<double> lv2 = out_verts[low_vert_idx2];

    // vit2 is the location to insert p2
    std::vector<vgl_point_2d<double> >::iterator vit2 = std::find(out_verts.begin(), out_verts.end(), lv2);
    if (vit2 == out_verts.end()) {
      std::cout << "insert p2 vert location not found -- fatal" << std::endl;
      return false;
    }
    bool first1_inside = false, first2_inside = false;
    for (std::vector<std::pair<size_t, size_t> >::const_iterator iit = matched_edges.begin();
         iit != matched_edges.end(); ++iit) {
      if (*iit == p1_split_edge|| *iit == p2_split_edge)
        continue;
      if (vidxf == iit->first || vidxf == iit->second)
        first1_inside = true;
      if (vidxf2 == iit->first || vidxf2 == iit->second)
        first2_inside = true;
    }
    for (std::set<size_t>::iterator uit = unique_vert_indices.begin();
         uit != unique_vert_indices.end(); ++uit) {
      if ((first1_inside && *uit == vidxs) || (!first1_inside && *uit == vidxf)||
          (first2_inside && *uit == vidxs2) || (!first2_inside && *uit == vidxf2))
        continue;
      to_erase.push_back(out_verts[*uit]);
    }

    out_verts.insert(vit, p1); // insert after finding verts to erase
    out_verts.insert(vit2, p2);

    for (std::vector<vgl_point_2d<double> >::iterator rit = to_erase.begin();
         rit != to_erase.end(); ++rit) {
      std::vector<vgl_point_2d<double> >::iterator eit = std::find(out_verts.begin(), out_verts.end(), *rit);
      out_verts.erase(eit);
    }
    mcr_2d.outer_cycle_ = out_verts;
    return true;
  }
  return true;
}


bool vkm_ground_truth::ensure_consistent_topology(
    size_t outer_index, mc_region_2d& mcr_2d)
{
  //
  // It is anticipated that a number of annotation issues will be repaired by this function
  //
  //
  double tol = 0.15;
  std::vector<vgl_line_segment_2d<double> > outer_edge_segs;
  std::vector<std::pair<size_t, size_t> > outer_edge_vert_ids;
  size_t nout = mcr_2d.outer_cycle_.size();
  for (size_t i = 0; i < nout; ++i) {
    size_t ip1 = (i + 1) % nout;
    std::pair<size_t, size_t> pr(i, ip1);
    outer_edge_vert_ids.push_back(pr);
    vgl_line_segment_2d<double> seg(mcr_2d.outer_cycle_[i], mcr_2d.outer_cycle_[ip1]);
    outer_edge_segs.push_back(seg);
  }

  for (std::map<size_t, std::vector<vgl_point_2d<double> > >::const_iterator hit = mcr_2d.holes_.begin();
       hit != mcr_2d.holes_.end(); ++hit) {
    size_t hole_idx = hit->first;
    // extract hole edges
    const std::vector<vgl_point_2d<double> >& hole_verts = hit->second;
    std::vector<std::pair<size_t, size_t> > matched_outer_edge_vert_ids;
    std::vector<vgl_line_segment_2d<double> > unmatched_hole_edge_segs;
    std::vector<std::pair<size_t, size_t> > unmatched_hole_edge_vert_ids;
    size_t nh = hole_verts.size();
    for (size_t j = 0; j < nh; ++j) {
      size_t jp1 = (j + 1) % nh;
      vgl_line_segment_2d<double> seg(hole_verts[j], hole_verts[jp1]);
      std::pair<size_t, size_t> segs_match;
      if (!matched(seg, outer_edge_segs, outer_edge_vert_ids, segs_match, tol)) {
        std::pair<size_t, size_t> pr(j, jp1);
        unmatched_hole_edge_segs.push_back(seg);
        unmatched_hole_edge_vert_ids.push_back(pr);
      }
      else {
        matched_outer_edge_vert_ids.push_back(segs_match);
      }
    }
    size_t nm = unmatched_hole_edge_segs.size();
    if (nm != 1)
      continue;
    else{
      // Case I An annotator draws an incomplete boundary on the base surface for a planar
      // region above, e.g. a stair well right at the edge of a roof.
      // The outer cycle is the base region, and the hole is generated by the region above the base.
      //    ----------------
      //    |  outer cycle |
      //    |  p1          |
      //    --0----------  |
      //      * + + + + +| |
      //      * +  hole +| |
      //      * + + + + +| |
      //    --0----------  |
      //    |  p2          |
      //    |              |
      //    ---------------
      // a hole region coincident with the boundary of the parent but not actually enclosed
      // in this case the outer cycle should be modified.
      // a) a new edge or edges added to the outer one chain to enclose the hole (the edge with stars)
      // b) remove the parts of the outer boundary that are now inside the added boundary
      //
      // Note that region labeled "hole" is actually outside the outer cycle but all its vertices touch
      // the outer cycle boundary and thus are considered "inside", which is a correct result
      std::cout << "In face " << outer_index << " Hole " << hole_idx << " a partial match of " << nm << " not matched and " << nh - nm << "  matched out of " << nh << " hole edges" << std::endl;
      repair_topology_caseI(mcr_2d, matched_outer_edge_vert_ids, unmatched_hole_edge_segs[0], hole_idx, tol);
    }//end CaseI
  }
  return true;
}


void vkm_ground_truth::construct_polygon_soup()
{
  for(std::map<size_t, vgl_plane_3d<double> >::iterator rit = region_planes_.begin();
      rit != region_planes_.end(); ++rit){
    size_t indx = rit->first;

    // copy region since it might be repaired
    mc_region_2d mcr_2d = cont_tree_.mc_regions()[indx];

    // potentially alter the topology of the region to repair inconsistencies
    if(!ensure_consistent_topology(indx, mcr_2d))
      continue;

    //project the region vertices onto the fitted plane
    const vgl_plane_3d<double>& pl = rit->second;
    mc_region_3d mcr_3d;
    mcr_3d.enclosing_region_ = mcr_2d.enclosing_region_;
    size_t nout = mcr_2d.outer_cycle_.size();
    for(size_t vi =0; vi<nout; ++vi){
      unsigned u = static_cast<unsigned>(mcr_2d.outer_cycle_[vi].x()), v =static_cast<unsigned>( mcr_2d.outer_cycle_[vi].y());
      vgl_point_3d<double> p3d = dem_to_world(u, v);
      vgl_point_3d<double> proj_z;
      if(!proj_vertical(p3d, pl, proj_z)){
        std::cout << " index = " << indx << std::endl;
        continue;
      }
      mcr_3d.outer_cycle_.push_back(proj_z);
    }
    for(std::map<size_t, std::vector<vgl_point_2d<double> > >::const_iterator hit = mcr_2d.holes_.begin();
        hit != mcr_2d.holes_.end(); ++hit){
      size_t hindx = hit->first;
      const std::vector<vgl_point_2d<double> >& hole = hit->second;
      size_t nh = hole.size();
      for(size_t vh = 0; vh<nh; ++vh){
        unsigned u = static_cast<unsigned>(hole[vh].x()), v =static_cast<unsigned>(hole[vh].y());
        vgl_point_3d<double> p3d = dem_to_world(u, v);
        vgl_point_3d<double> proj_z;
        if(!proj_vertical(p3d, pl, proj_z)){
          std::cout << " index = " << indx << std::endl;
          continue;
        }
        mcr_3d.holes_[hindx].push_back(proj_z);
      }
    }
    regions_3d_[indx]=mcr_3d;
  }
}


void vkm_ground_truth::convert_to_meshes()
{
  for(std::map<size_t, mc_region_3d >::iterator rit = regions_3d_.begin();
      rit != regions_3d_.end(); ++rit){
    size_t index = rit->first;
    const mc_region_3d& mcr_3d = rit->second;
    // the outer boundary cycle
	size_t n = mcr_3d.outer_cycle_.size();
    if(n == 0){
      std::cout << "Region " << index << " has no boundaries - skip" << std::endl;
      continue;
    }
    PolyMesh mesh;
    std::vector<PolyMesh::VertexHandle> vhandles;
    for(size_t i = 0; i<n; ++i){
      const vgl_point_3d<double>& p = mcr_3d.outer_cycle_[i];
      OpenMesh::Vec3f vf(static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()));
      vhandles.push_back(mesh.add_vertex(vf));
    }
    mesh.add_face(vhandles);
    region_meshes_[index][index] = mesh;
    // add in the holes
    for(std::map<size_t, std::vector<vgl_point_3d<double> > >::const_iterator hit = mcr_3d.holes_.begin();
        hit != mcr_3d.holes_.end(); ++hit){
      PolyMesh hmesh;
      vhandles.clear();
      size_t hindx = hit->first;
      const std::vector<vgl_point_3d<double> >& hole = hit->second;;
      size_t nh = hole.size();
      for(unsigned i =0; i<nh; ++i){
        const vgl_point_3d<double>& p = hole[i];
        OpenMesh::Vec3f vf(static_cast<float>(p.x()), static_cast<float>(p.y()), static_cast<float>(p.z()));
        vhandles.push_back(hmesh.add_vertex(vf));
      }
      hmesh.add_face(vhandles);
      region_meshes_[index][hindx] = hmesh;
    }
  }
}


bool vkm_ground_truth::region_contained(
    mc_region_3d const& reg, vgl_polygon<double> const& perim)
{
  mc_region_2d mcr_2d = reg.xy_region();
  const std::vector<vgl_point_2d<double> >& verts = mcr_2d.outer_cycle_;
  bool inside = true;
  for(std::vector<vgl_point_2d<double> >::const_iterator vit = verts.begin();
      vit != verts.end()&&inside; ++vit)
    inside = perim.contains(*vit);
  return inside;
}


// map 3d regions onto x-y plane
bool vkm_ground_truth::write_xy_polys(
    std::string const& site_name, std::string const& path)
{
  const vgl_polygon<double>& perim = perimeters_[site_name];
  if(perim.num_sheets() == 0){
    std::cout << "no perimeter defined for " << site_name << std::endl;
    return false;
  }
  std::ofstream ostr(path.c_str());
    if(!ostr){
      std::cout << "Failed to open " << path << " for writing xy regions" << std::endl;
      return false;
    }
	std::vector < vgl_polygon<double> > out_polys;
    for(std::map<size_t, mc_region_3d>::iterator rit =  regions_3d_.begin();
        rit != regions_3d_.end(); ++rit){
      mc_region_2d mcr_2d = rit->second.xy_region();
      const std::vector<vgl_point_2d<double> >& verts = mcr_2d.outer_cycle_;
      bool inside = true;
      for(std::vector<vgl_point_2d<double> >::const_iterator vit = verts.begin();
          vit != verts.end()&&inside; ++vit)
        inside = perim.contains(*vit);
      if(!inside)
        continue;
      vgl_polygon<double> poly = mcr_2d.poly();
	  out_polys.push_back(poly);
    }
    ostr << out_polys.size() << std::endl;
    for (std::vector < vgl_polygon<double> >::iterator oit = out_polys.begin();
         oit != out_polys.end(); ++oit)
      ostr << *oit;

    ostr.close();
    return true;
}


bool vkm_ground_truth::write_processed_ground_truth(std::string const& path)
{
  // TODO - write holes, deal with intersection over union and holes
  // for now, only write 1st minor mesh (parent object)

  // flatten structure (map of MeshT pointers)
  std::stringstream ss;
  std::map<std::string, const PolyMesh*> groups;
  for (const auto& major : region_meshes_ ) {
    for (const auto& minor : major.second ) {
      ss.str("");
      ss << major.first << "_" << minor.first << " "
         << type_to_string(surface_types_[major.first]);
      groups[ss.str()] = &(minor.second);
      break;
    }
  }

  // write via private function
  return vkm_obj_io::write_group_pointer_obj_file(path,groups);
}


bool vkm_ground_truth::load_processed_ground_truth(
    std::string const& path,
    std::map<size_t, std::map<size_t, PolyMesh> >& region_meshes,
    std::map<size_t, surface_t>& surface_types)
{
  // get region meshes
  std::map<size_t, std::map<size_t, std::string> > properties;
  if (!vkm_obj_io::read_composite_obj_file(path,region_meshes,properties))
    return false;

  // custom parse surface-type property
  for (const auto& p1 : properties)
    for (const auto& p2 : p1.second)
      surface_types[p1.first] = string_to_type(p2.second);

  return true;


  // std::ifstream istr(path.c_str());
  // if(!istr){
  //   std::cout << "Failed to open " << path << " for reading obj file" << std::endl;
  //   return false;
  // }
  // std::string gs, gids, str;
  // size_t start_id = 1;
  // size_t cum_id = 1;
  // std::string surface_type_str;
  // std::map<size_t, size_t> vmap;
  // std::map<size_t, PolyMesh>  cur_mesh_group;
  // while (istr >> gs >> gids >> surface_type_str) {
  //   if (gs != "g") {
  //     std::cout << "Parse failed " << gs << ' ' << gids << std::endl;
  //     return false;
  //   }
  //   std::string majs, mins;
  //   size_t outer_id = 0, hole_id =0;
  //   bool reading_major = true;
  //   for (std::string::iterator sit = gids.begin();
  //        sit != gids.end(); ++sit) {
  //     if (reading_major&& *sit != '_') {
  //       majs.push_back(*sit);
  //       continue;
  //     }
  //     else if (*sit == '_') {
  //       reading_major = false;
  //       continue;
  //     }
  //     else if (reading_major == false) {
  //       mins.push_back(*sit);
  //     }
  //   }
  //   std::stringstream maj_ss(majs), min_ss(mins);
  //   maj_ss >> outer_id;
  //   min_ss >> hole_id;
  //   surface_types[outer_id] = string_to_type(surface_type_str);
  //   PolyMesh mesh;
  //   // read the mesh vertices
  //   std::vector< PolyMesh::VertexHandle> in_vhandles;
  //   std::string type;
  //   while(true){
  //     double x, y, z;
  //     istr >> type;
  //     if(type != "v"){
  //       break;
  //     }
  //     istr >> x >> y >> z;
  //     PolyMesh::Point p(x, y, z);
  //     in_vhandles.push_back(mesh.add_vertex(p));
  //     cum_id++;
  //   }//end while(true)

  //   // read the face
  //   if(type != "f"){
  //     std::cout << "Parse failed " << type << std::endl;
  //     return false;
  //   }
  //   std::vector< PolyMesh::VertexHandle> vhandles;
  //   size_t k = 0;
  //   for(size_t i = start_id; i<cum_id ; ++i, ++k){
  //     size_t vh;
  //     istr >> vh;
  //     vhandles.push_back(in_vhandles[k]);
  //   }
  //   start_id = cum_id;
  //   mesh.add_face(vhandles);
  //   // end of face read
  //   region_meshes[outer_id][hole_id] = mesh;
  // } // end of file

  // return true;
}


void vkm_ground_truth::convert_img_regions_to_meshes()
{
  for(std::map<size_t, std::pair<vgl_polygon<double>, vgl_box_2d<double> >  >::iterator rit = img_regions_.begin();
      rit != img_regions_.end(); ++rit){
    size_t index = rit->first;
    const std::vector<vgl_point_2d<double> >& verts = rit->second.first[0];
    size_t n = verts.size();
    if(n == 0)
      continue;
    PolyMesh mesh;
    std::vector<PolyMesh::VertexHandle> vhandles;
    for(size_t i = 0; i<n; ++i){
      unsigned u = static_cast<unsigned>(verts[i].x()), v =static_cast<unsigned>(verts[i].y());
      vgl_point_3d<double> p3d = dem_to_world(u, v);
      OpenMesh::Vec3f vf(static_cast<float>(p3d.x()), static_cast<float>(p3d.y()), static_cast<float>(p3d.z()));
      vhandles.push_back(mesh.add_vertex(vf));
    }
    mesh.add_face(vhandles);
    img_meshes_[index][index] = mesh;
  }
}


bool vkm_ground_truth::load_site_perimeter(std::string const& site_name, std::string const& bwm_ptset_path)
{
  // currently defined by a bwm pointset ascii file in dem image coordinates - replace later with Kitware format
  std::ifstream istr(bwm_ptset_path.c_str());
  if(!istr){
    std::cout << "Failed to open " << bwm_ptset_path << " for reading footprint" << std::endl;
    return false;
  }
  size_t npts;
  istr >> npts;
  if(npts == 0){
    std::cout << "no perimeter points - fatal" << std::endl;
    return false;
  }
  std::vector<vgl_point_2d<double> > verts;
  double avg_z = 0.0;
  for(size_t i = 0; i<npts; ++i){
    double ud , vd;
    istr >> ud >> vd >> std::ws;
    unsigned u = static_cast<unsigned>(ud), v = static_cast<unsigned>(vd);
    vgl_point_3d<double> p3d = dem_to_world(u, v);
    verts.push_back(vgl_point_2d<double>(p3d.x(), p3d.y()));
    avg_z += p3d.z();
  }
  avg_z /= npts;
  local_gnd_planes_[site_name] = vgl_plane_3d<double>(0.0, 0.0, 1.0, -avg_z);
  vgl_polygon<double> poly(verts);
  perimeters_[site_name] = poly;
  return true;
}


vkm_ground_truth::PolyMesh::FaceHandle vkm_ground_truth::extrude_base_pts(
    const std::vector<vgl_point_3d<double> >& base_points,
    vgl_plane_3d<double> ground_plane, PolyMesh& mesh)
{
  size_t npts = base_points.size();
  std::vector<PolyMesh::Point> base_ompts, end_ompts;
  for(size_t i = 0; i<npts; ++i){
    vgl_point_3d<double> cp = vgl_closest_point(base_points[i], ground_plane);
    if(cp.z()>=base_points[i].z()){
      std::cout << "warning, in ExtrudedShape to gnd plane, closest point on ground plane " << cp << std::endl
                << "is above the base point " << base_points[i] << " clipping base point elevation " << std::endl;

      base_ompts.push_back(PolyMesh::Point(static_cast<float>(base_points[i].x()),
                                           static_cast<float>(base_points[i].y()),
                                           static_cast<float>(cp.z())));
    }else
      base_ompts.push_back(PolyMesh::Point(static_cast<float>(base_points[i].x()),
                                           static_cast<float>(base_points[i].y()),
                                           static_cast<float>(base_points[i].z())));

    end_ompts.push_back(PolyMesh::Point(static_cast<float>(cp.x()),
                                        static_cast<float>(cp.y()),
                                        static_cast<float>(cp.z())));
  }
  std::vector<PolyMesh::VertexHandle> base_vhandles, end_vhandles, temp;
  for(size_t i = 0; i<npts; ++i){
    PolyMesh::VertexHandle bvh = mesh.add_vertex(base_ompts[i]);
    base_vhandles.push_back(bvh);
    PolyMesh::VertexHandle evh = mesh.add_vertex(end_ompts[i]);
    end_vhandles.push_back(evh);
  }

  PolyMesh::FaceHandle fh = mesh.add_face(base_vhandles);

  for(int i = static_cast<int>(npts-1); i>=0; --i)
      temp.push_back(end_vhandles[i]);
  mesh.add_face(temp);

  for(size_t i = 0; i<npts; ++i){
    temp.clear();
    temp.push_back(base_vhandles[i]);
    if(base_vhandles[i]!=end_vhandles[i])
      temp.push_back(end_vhandles[i]);

    temp.push_back(end_vhandles[(i+1)%npts]);
    if(base_vhandles[(i+1)%npts]!=end_vhandles[(i+1)%npts])
      temp.push_back(base_vhandles[(i+1)%npts]);

    mesh.add_face(temp);
  }
  return fh;
}


bool vkm_ground_truth::construct_extruded_gt_model(std::string const& name)
{
  const vgl_polygon<double>& perim = perimeters_[name];
  if(perim.num_sheets() == 0){
    std::cout << "no perimeter defined for site " << name << std::endl;
    return false;
  }

  std::map<std::string, vgl_plane_3d<double> >::iterator git;
  git =  local_gnd_planes_.find(name);
  if(git == local_gnd_planes_.end()){
    std::cout << "no ground plane for site " << name << std::endl;
    return false;
  }
  vgl_plane_3d<double> local_ground_plane = git->second;
  std::map<size_t, PolyMesh> mesh_map;

  // only outer cycle for now - handle holes with mc_mesh later
  for(std::map<size_t, mc_region_3d>::iterator rit =  regions_3d_.begin();
      rit != regions_3d_.end(); ++rit){
    size_t index = rit->first;
    if(!region_contained(rit->second, perim))
      continue;
    PolyMesh mesh;

    const std::vector<vgl_point_3d<double> >& base_points = rit->second.outer_cycle_;
    size_t npts = base_points.size();
    if(npts<3){
      std::cout << " warning - insufficient number of points " << npts << " for index " << index <<  " can't create mesh " << std::endl;
      continue;
    }
    vgl_plane_3d<double> ground_plane = local_ground_plane;
    size_t rid = rit->second.enclosing_region_;
    if(rid != -1){
      std::map<size_t, vgl_plane_3d<double> >::iterator pit = region_planes_.find(rid);
      if(pit == region_planes_.end()){
        std::cout << " warning - enclosing plane doesn't exist  for region index " << index <<  " can't create mesh " << std::endl;
        continue;
      }
      ground_plane = pit->second;
    }
    extrude_base_pts(base_points, ground_plane, mesh);
    mesh_map[index] = mesh;
  }
  extruded_regions_[name] = mesh_map;
  return true;
}


bool vkm_ground_truth::write_extruded_gt_model(std::string const& name, std::string const& path)
{
  std::map<std::string, std::map<size_t, PolyMesh> >::iterator mit = extruded_regions_.find(name);
  if(mit == extruded_regions_.end()){
    std::cout << "site name " << name << " doesn't exist - can't write " << std::endl;
    return false;
  }
  std::map<size_t, PolyMesh>& mesh_map = mit->second;
  if(!vkm_obj_io::write_composite_obj_file(path, mesh_map))
    return false;
  return true;
}


bool vkm_ground_truth::write_ground_planes(std::string const& gnd_plane_path) const
{
  std::ofstream ostr(gnd_plane_path.c_str());
  if(!ostr){
    std::cout << "can't open " << gnd_plane_path << " to save site ground planes" << std::endl;
    return false;
  }
  for(std::map<std::string, vgl_plane_3d<double> >::const_iterator git = local_gnd_planes_.begin();
      git != local_gnd_planes_.end(); ++git)
    ostr << git->first << ' ' << git->second << std::endl;
  ostr.close();
  return true;
}
