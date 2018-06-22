#include <sstream>
#include <assert.h>
#include <fstream>
#include <iomanip>

// VXL ------------------------------------
#include <vgl/algo/vgl_fit_plane_3d.h>
#include <vgl/vgl_intersection.h>
#include <vgl/vgl_clip.h>
#include <vgl_area.h>
#include <vnl/vnl_math.h>
//-----------------------------------------

#include "vkm_metrics.h"


xy_region::xy_region(std::map<size_t, PolyMesh> const& map)
{
  assert(map.size()>0);
  std::map<size_t, PolyMesh>::const_iterator mit = map.begin();
  PolyMesh mesh = mit->second;
  std::vector<vgl_point_3d<double> > verts_3d;
  for(PolyMesh::VertexIter vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit){
    OpenMesh::Vec3f p = mesh.point(*vit);
    verts_3d.push_back(vgl_point_3d<double>(p[0], p[1], p[2]));
  }
  // compute the planar surface
  size_t n = verts_3d.size();
  std::vector<vgl_homg_point_3d<double> > plane_pts;
  for (size_t i = 0; i < n; ++i){
    vgl_homg_point_3d<double> hp3d(verts_3d[i]);
    plane_pts.push_back(hp3d);
  }
  vgl_fit_plane_3d<double> fitter(plane_pts);
  double tol = 0.25;
  assert(fitter.fit(tol, &std::cout));
  vgl_homg_plane_3d<double> plf = fitter.get_plane();
  plf.normalize();
  plane_.set(plf.a(), plf.b(), plf.c(), plf.d());
  // form the (potentially multiply connected) boundary
  vgl_polygon<double> poly;
  for(std::map<size_t, PolyMesh>::const_iterator mit = map.begin();
      mit != map.end(); ++mit){
    std::vector<vgl_point_2d<double> > verts_2d;
    for(PolyMesh::VertexIter vit=mesh.vertices_begin(); vit!=mesh.vertices_end(); ++vit){
      OpenMesh::Vec3f p = mesh.point(*vit);
      vgl_point_2d<double> p2d(p[0], p[1]);
      verts_2d.push_back(p2d);
      bounding_box_.add(p2d);
    }
    poly.push_back(verts_2d);
  }
  region_ = poly;
}


typename xy_region::PolyMesh xy_region::mesh()
{
  PolyMesh ret;
  for(size_t js = 0; js<region_.num_sheets(); ++js){
    std::vector<OpenMesh::VertexHandle> vhandles;
    for(size_t i = 0; i< region_[js].size(); ++i){
      const vgl_point_2d<double>& vert = region_[js][i];
      double a = plane_.a(), b = plane_.b(), c = plane_.c(), d = plane_.d();
      double z = 0.0;
      if(fabs(c)>0.001)//vertical plane so z is undefined
        z = -(a*vert.x() + b*vert.y() + d)/c;
      OpenMesh::Vec3f v(vert.x(), vert.y(), z);
      vhandles.push_back(ret.add_vertex(v));
    }
    ret.add_face(vhandles);
  }
  return ret;
}


std::vector<vgl_point_3d<double> > xy_region::verts_3d() const
{
  std::vector<vgl_point_3d<double> > ret;
  size_t n = region_[0].size();
  for (size_t i = 0; i < n; ++i) {
    const vgl_point_2d<double>& vert = region_[0][i];
    double a = plane_.a(), b = plane_.b(), c = plane_.c(), d = plane_.d();
    double z = 0.0;
    if(fabs(c)>0.001)//vertical plane so z is undefined
      z = -(a*vert.x() + b*vert.y() + d)/c;
    ret.push_back(vgl_point_3d<double>(vert.x(), vert.y(), z));
  }
  return ret;
}


void xy_region::print_poly(vgl_polygon<double> const& poly)
{
  for(size_t js = 0; js<poly.num_sheets(); ++js){
    std::cout << "sheet " << js << std::endl;
    for(size_t i = 0; i< poly[js].size(); ++i)
      std::cout << poly[js][i].x() << ' ' << poly[js][i].y() << std::endl;
  }
}


xy_region xy_region::translate(double tx, double ty, double tz) const
{
  xy_region ret;
  //translate verts, compute new bounding box
  for(size_t js = 0; js<region_.num_sheets(); ++js){
    std::vector<vgl_point_2d<double> > tsheet;
    for(size_t i = 0; i< region_[js].size(); ++i){
      const vgl_point_2d<double>& vert = region_[js][i];
      vgl_point_2d<double> tvert(vert.x()+tx, vert.y()+ty);
      ret.bounding_box_.add(tvert);
      tsheet.push_back(tvert);
    }
    ret.region_.push_back(tsheet);
  }
  double a = plane_.a(), b = plane_.b(), c = plane_.c(), d = plane_.d();
  ret.plane_.set(a, b, c, (d -(a*tx + b*ty + c*tz)));
  return ret;
}


vgl_point_2d<double> xy_region::xy_centroid() const
{
  vgl_point_2d<double> ret(0.0, 0.0);
  double xc = 0, yc = 0, nd = 0;
  for(size_t js = 0; js<region_.num_sheets(); ++js){
    std::vector<vgl_point_2d<double> > tsheet;
    for(size_t i = 0; i< region_[js].size(); ++i){
      const vgl_point_2d<double>& vert = region_[js][i];
      xc += vert.x(); yc += vert.y(); nd +=1.0;
    }
  }
  if(nd == 0.0)
    return ret;
  xc /= nd;  yc /= nd;
  ret.set(xc, yc);
  return ret;
}


double xy_region::z_off(
    vgl_plane_3d<double> const& pl0,
    vgl_plane_3d<double> const& pl1,
    vgl_point_2d<double> const& p)
{
    double n0x = pl0.a(), n0y = pl0.b(), n0z = pl0.c(), d0 = pl0.d();
    double n1x = pl1.a(), n1y = pl1.b(), n1z = pl1.c(), d1 = pl1.d();
    double tempx = p.x()*(n1z*n0x-n0z*n1x);
    double tempy = p.y()*(n1z*n0y-n0z*n1y);
    double tempd = (n1z*d0 - n0z*d1);
    double z = (tempx + tempy +tempd)/(n1z*n0z);
    return z;
}


// encode json string for score instance
void score_json_helper(std::stringstream &ss, std::string name, double val) {
  ss << "\"" << name << "\": ";
  if (std::isfinite(val))
    ss << val;
  else
    ss << "null";
}

std::string score::json() const
{
  std::stringstream ss;
  ss << "{ ";
  ss << "\"gt_region_id\": " << gt_region_id_ << ", ";
  ss << "\"test_region_id\": " << test_region_id_ << ", ";
  score_json_helper(ss,"area",area_); ss << ", ";
  score_json_helper(ss,"comp",comp_); ss << ", ";
  score_json_helper(ss,"corr",corr_); ss << ", ";
  score_json_helper(ss,"iou",iou_); ss << ", ";
  score_json_helper(ss,"normal_ang_diff",normal_ang_diff_); ss << ", ";
  score_json_helper(ss,"z_error",z_error_);
  ss << " }";
  return ss.str();
}


score operator +(score const& sa, score const& sb)
{
  score ret;
  ret.comp_ = sa.comp_ + sb.comp_;
  ret.corr_ = sa.corr_ + sb.corr_;
  ret.iou_ = sa.iou_ + sb.iou_;
  return ret;
}


score operator /=(score& a, double val)
{
  a.comp_ /=val;
  a.corr_/= val;
  a.iou_/= val;
  return a;
}


std::ostream& operator<<(std::ostream& str, score const& s)
{
  str <<"(comp corr iou) " << s.comp_ << ' ' << s.corr_ << ' ' << s.iou_;
  return str;
}


bool vkm_metrics::load_ground_planes(std::string const& path)
{
  std::ifstream istr(path.c_str());
  if(!istr){
    std::cout << "Can't load site ground planes from " << path << std::endl;
    return false;
  }
  std::string site_name;
  vgl_plane_3d<double> gnd_plane;
  while(istr >> site_name >> gnd_plane)
    local_gnd_planes_[site_name] = gnd_plane;
  return true;
}


void vkm_metrics::construct_xy_regions()
{
  // for now test model is only simply connected
  for(std::map<size_t, std::map<size_t,PolyMesh> >::iterator mit = test_model_.begin();
      mit != test_model_.end(); ++mit)
    test_regions_[mit->first] = xy_region(mit->second);

  // the gt could be multiply connected
  for (std::map<size_t, std::map<size_t, PolyMesh> >::iterator git = gt_region_meshes_.begin();
	  git != gt_region_meshes_.end(); ++git) {
	  gt_regions_[git->first] = xy_region(git->second);
  }
}


void vkm_metrics::translate_test_model_xy()
{
  //  registered_test_footprint_ = test_footprint_.translate(tx_, ty_);
  for(std::map<size_t, xy_region>::const_iterator rit = test_regions_.begin();
      rit != test_regions_.end(); ++rit)
    registered_test_regions_[rit->first] = (rit->second).translate(tx_, ty_);
}


void vkm_metrics::match_xy_regions()
{
  for(std::map<size_t, xy_region>::const_iterator git = gt_regions_.begin();
      git != gt_regions_.end(); ++git){
    size_t gindx = git->first;
    const vgl_box_2d<double>& gbb = git->second.bounding_box_;
    const vgl_polygon<double>& gpoly = git->second.region_;
    for(std::map<size_t, xy_region>::const_iterator rit = registered_test_regions_.begin();
      rit != registered_test_regions_.end(); ++rit){
      size_t tindx = rit->first;
      const vgl_box_2d<double>& tbb = rit->second.bounding_box_;
      const vgl_polygon<double>& tpoly = rit->second.region_;
      vgl_box_2d<double> bint = vgl_intersection(gbb, tbb);
      if(bint.is_empty())
        continue;
      // potential intersection
      vgl_polygon<double> int_poly = vgl_clip(gpoly, tpoly);
      double int_area = vgl_area(int_poly);
      if (int_area == 0.0)
        continue;
      vgl_polygon<double> union_poly = vgl_clip(gpoly, tpoly, vgl_clip_type_union);
       double union_area = vgl_area(union_poly);
      if(union_area < vgl_area(gpoly) || union_area < vgl_area(tpoly))
        continue;
      double iou = int_area/union_area;
      if(iou<params_.min_iou_)
        continue;
      double tp_area = int_area;
      double gpoly_area = vgl_area(gpoly);
      double comp = tp_area/gpoly_area;
      double tpoly_area = vgl_area(tpoly);
      double corr = tp_area/tpoly_area;
      score s;
      s.comp_ = comp;
      s.corr_ = corr;
      s.iou_ = iou;
      gt_to_test_[gindx][tindx]= s;
      test_to_gt_[tindx][gindx] = s;
    }
  }
}


void vkm_metrics::compute_best_match_2d_score_stats()
{
  double dn = 0.0;
  score avg_score;
  score score_h;
  for(std::map<size_t, std::map<size_t, score> >::iterator git = gt_to_test_.begin();
      git != gt_to_test_.end(); ++git, dn+=1.0){
    size_t j = git->first; // ground truth region index
    double max_iou = 0.0;
    score max_score;
    size_t max_region_indx;
    for(std::map<size_t, score>::iterator mit = git->second.begin();
        mit != git->second.end(); ++mit){
      size_t indx = mit->first; // test region index
      score s = mit->second;
      if(s.iou_ > max_iou){
        max_iou = s.iou_;
        max_score = s;
        max_region_indx = indx; // test region index
      }
    }
    max_score.gt_region_id_ = j;
    max_score.test_region_id_ = max_region_indx;
    //                           test region
    std::pair<size_t, score> pr(max_region_indx, max_score);
    max_scores_[j] = pr;
    avg_score = avg_score + max_score;
  }
  avg_score /= dn;
  avg_best_match_score_ = avg_score;
  std::cout << " avg_best_match " << avg_score << std::endl;
}


void vkm_metrics::print_score_array()
{
  size_t width = 12;
  std::ostringstream oss;

  oss << std::setw(width) << "gt_indx"
      << std::setw(width) << "test_indx"
      << std::setw(width) << "area"
      << std::setw(width) << "comp"
      << std::setw(width) << "corr"
      << std::setw(width) << "iou"
      << std::setw(width) << "ang_er"
      << std::setw(width) << "z_error"
      ;
  std::cout << oss.str() << std::endl;

  for(std::map<size_t, std::pair<size_t, score > >::iterator sit = max_scores_.begin();
      sit != max_scores_.end(); ++sit){
    const score& s = sit->second.second;

    oss.str(""); oss.clear();

    oss << std::fixed;
    oss << std::setw(width) << s.gt_region_id_;
    oss << std::setw(width) << s.test_region_id_;
    oss << std::setprecision(3);
    oss << std::setw(width) << std::fixed << s.area_;
    oss << std::setprecision(6);
    oss << std::setw(width) << s.comp_;
    oss << std::setw(width) << s.corr_;
    oss << std::setw(width) << s.iou_;
    oss << std::setw(width) << s.normal_ang_diff_;
    oss << std::setw(width) << s.z_error_;
    std::cout << oss.str() << std::endl;
  }
}


void vkm_metrics::find_z_offset()
{
  z_off_ = 0.0;
  double iou_sum = 0.0;
  std::vector<std::pair<double, double> > area_zoff_vals;
  for(std::map<size_t, std::pair<size_t, score> >::iterator mit = max_scores_.begin();
      mit != max_scores_.end(); ++mit){
    size_t gt_indx = mit->first;
    size_t test_indx = mit->second.first;
    double iou = mit->second.second.iou_;
    const xy_region& gt_reg = gt_regions_[gt_indx];
    const xy_region& test_reg = registered_test_regions_[test_indx];
    vgl_plane_3d<double> pl_gt = gt_reg.plane_, pl_test = test_reg.plane_;
    // the intersection region
    const vgl_polygon<double>& gt_poly = gt_reg.region_;
    const vgl_polygon<double>& test_poly = test_reg.region_;
    xy_region temp;
    vgl_polygon<double> int_poly = vgl_clip(gt_poly, test_poly);
    temp.region_ = int_poly;
    // centroid of intersection region
    vgl_point_2d<double> c = temp.xy_centroid();
    // compute z offset
    double z = xy_region::z_off(pl_test, pl_gt, c);
    std::pair<double, double> pr(vgl_area(int_poly), z);
    area_zoff_vals.push_back(pr);
    z_off_ += z*iou;
    iou_sum += iou;
  }
  if(iou_sum == 0.0){
    z_off_ = 0.0;
    std::cout << "z offset failed " <<  std::endl;
    return;
  }
  std::cout << "debugging: z-off/iou_sum " << z_off_ << "/" << iou_sum << std::endl;
  z_off_ = z_off_/iou_sum;
  std::cout << "z offset " << z_off_ << " meters" << std::endl;
  // transform regions. assume  x-y trans already known.
  for(std::map<size_t, xy_region>::iterator rit = test_regions_.begin();
      rit != test_regions_.end(); ++rit){
    size_t indx = rit->first;
    xy_region trans = rit->second.translate(tx_, ty_, z_off_);
    registered_test_regions_[indx] = trans;
  }
  // compute avg z error
  double area_sum = 0.0, esum = 0.0;
  for(std::vector<std::pair<double, double> >::iterator ait = area_zoff_vals.begin();
      ait != area_zoff_vals.end(); ++ait){
    area_sum += ait->first;
    esum += fabs(ait->second - z_off_)*ait->first;
  }
  double avg_z_er = esum/area_sum;
  std::cout << "Average z error after registration " << avg_z_er << " meters (based on TP regions only)"<< std::endl;
}


bool vkm_metrics::save_transformed_regions_as_meshes(std::string const& path, std::string site_name) const
{
  std::map<std::string, vgl_plane_3d<double> >::const_iterator git;
  git = local_gnd_planes_.find(site_name);
  if(git == local_gnd_planes_.end()){
    std::cout << "site name " << site_name << " not found - can't save transformed regions" << std::endl;
    return false;
  }
  vgl_plane_3d<double> gnd_plane = git->second;
  std::map<size_t, PolyMesh> meshes;
  for(std::map<size_t, xy_region>::const_iterator rit =  registered_test_regions_.begin();
      rit != registered_test_regions_.end(); ++rit){
    size_t indx = rit->first;

    std::vector<vgl_point_3d<double> > pts3d = rit->second.verts_3d();
    PolyMesh mesh;
    vkm_ground_truth::extrude_base_pts(pts3d, gnd_plane, mesh);
    meshes[indx] = mesh;
  }
  std::string mat = "cyanMtl";
  std::string mat_file = "junk.mtl";
  return vkm_obj_io::write_composite_obj_file(path, meshes, mat_file, mat);
}


void vkm_metrics::compute_best_match_3d_score_stats()
{
  for(std::map<size_t, std::pair<size_t, score> >::iterator mit =  max_scores_.begin();
      mit !=  max_scores_.end();  ++mit){

    size_t gt_indx = mit->first;
    size_t test_indx = mit->second.first;
    score s = mit->second.second;
    double iou = mit->second.second.iou_;
    const xy_region& gt_reg = gt_regions_[gt_indx];
    const xy_region& test_reg = registered_test_regions_[test_indx];
    // the intersection region
    const vgl_polygon<double>& gt_poly = gt_reg.region_;
    double gt_area = vgl_area(gt_poly);
    if(gt_indx == 22)
      std::cout << ' ';
    s.area_ = gt_area;
    const vgl_polygon<double>& test_poly = test_reg.region_;
    xy_region temp;
    vgl_polygon<double> int_poly = vgl_clip(gt_poly, test_poly);
    temp.region_ = int_poly;
    // centroid of intersection region
    vgl_point_2d<double> c = temp.xy_centroid();
    // region planes
    vgl_plane_3d<double> pl_gt = gt_reg.plane_, pl_test = test_reg.plane_;
    s.z_error_ = xy_region::z_off(pl_test, pl_gt, c);
    // plane normal error - angle between vectors
    vgl_vector_3d<double> n_test = pl_test.normal(), n_gt = pl_gt.normal();
    if(n_test.z()<0) n_test *= -1.0;
    if(n_gt.z()<0) n_gt *= -1.0;
    double ang_rad = angle(n_test, n_gt);
    s.normal_ang_diff_ = ang_rad*vnl_math::deg_per_rad;
    max_scores_[gt_indx]=std::pair<size_t, score>(test_indx, s);
  }
}


std::string vkm_metrics::json() const
{
  std::stringstream ss;
  ss << "[" << std::endl;
  size_t n = max_scores_.size();
  size_t i = 0;
  for(std::map<size_t, std::pair<size_t, score > >::const_iterator sit = max_scores_.begin();
      sit != max_scores_.end(); ++sit, ++i)
  {
    ss << sit->second.second.json();
    if(i<(n-1)) {ss << ", ";}
    ss << std::endl;
  }
  ss << std::endl << "]" << std::endl;
  return ss.str();
}
