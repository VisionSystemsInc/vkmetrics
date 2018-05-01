#include <iostream>
#include <testlib/testlib_test.h>
#include <string>
#include <cmath>
#include <fstream>

// VXL ------------------------------------
#include <vgl/vgl_polygon.h>
#include <vil/vil_save.h>
// ----------------------------------------

#include "../vkm_binary_regions.h"


static void test_binary_regions()
{
  std::string site = "square_complex";
  //std::string site = "small_square_complex";
  //std::string site = "yasu_house";
  //std::string site = "firehouse";
  //std::string site = "directorate";

  std::string d2_dir = "D:/core3d_experiments/phase1a_D2/";
  std::string gt_dir = "D:/core3d_experiments/phase1a_D2/GTruth_data_v2/";

  std::map<std::string, std::string> xy_poly_path, gt_xy_poly_path, bin_image_path, gt_bin_image_path,  overlay_path;

  xy_poly_path["square_complex"] = d2_dir + "square_complex/xy_polys.txt";
  xy_poly_path["small_square_complex"] = d2_dir + "small_square_complex/xy_polys.txt";
  xy_poly_path["yasu_house"] = d2_dir + "yasu_house/xy_polys.txt";
  xy_poly_path["firehouse"] = d2_dir + "firehouse/xy_polys.txt";
  xy_poly_path["directorate"] = d2_dir + "directorate/xy_polys.txt";
  
  gt_xy_poly_path["square_complex"]= gt_dir + "square_complex_xy_polys.txt";
  gt_xy_poly_path["small_square_complex"]= gt_dir + "small_square_complex_xy_polys.txt";
  gt_xy_poly_path["yasu_house"]= gt_dir + "yasu_house_xy_polys.txt";
  gt_xy_poly_path["firehouse"]= gt_dir + "firehouse_xy_polys.txt";
  gt_xy_poly_path["directorate"]= gt_dir + "directorate_xy_polys.txt";
  
  bin_image_path["square_complex"] = d2_dir + "square_complex/footprint_bin_image.tif";
  bin_image_path["small_square_complex"] = d2_dir + "small_square_complex/footprint_bin_image.tif";
  bin_image_path["yasu_house"] = d2_dir + "yasu_house/footprint_bin_image.tif";
  bin_image_path["firehouse"] = d2_dir + "firehouse/footprint_bin_image.tif";
  bin_image_path["directorate"] = d2_dir + "directorate/footprint_bin_image.tif";

  gt_bin_image_path["square_complex"] = gt_dir + "square_complex_footprint_bin_image.tif";
  gt_bin_image_path["small_square_complex"] = gt_dir + "small_square_complex_footprint_bin_image.tif";
  gt_bin_image_path["yasu_house"] = gt_dir + "yasu_house_footprint_bin_image.tif";
  gt_bin_image_path["firehouse"] = gt_dir + "firehouse_footprint_bin_image.tif";
  gt_bin_image_path["directorate"] = gt_dir + "directorate_footprint_bin_image.tif";

  overlay_path["square_complex"] = gt_dir + "square_complex_overlay_image.tif";
  overlay_path["small_square_complex"] = gt_dir + "small_square_complex_overlay_image.tif";
  overlay_path["yasu_house"] = gt_dir + "yasu_house_overlay_image.tif";
  overlay_path["firehouse"] = gt_dir + "firehouse_overlay_image.tif";
  overlay_path["directorate"] = gt_dir + "directorate_overlay_image.tif";
  
  std::vector<vgl_polygon<double> > polys, gt_polys;
  std::ifstream istr(xy_poly_path[site].c_str());
  if(!istr)
    return;
  size_t n;
  istr >> n;
  for(size_t i = 0;i<n; ++i){
    vgl_polygon<double> poly;
    istr >> poly;
    size_t ns = poly.num_sheets();
    polys.push_back(poly);
  }
  istr.close();
  std::ifstream gt_istr(gt_xy_poly_path[site].c_str());
  if(!gt_istr)
    return;
  gt_istr >> n;
  for(size_t i = 0;i<n; ++i){
    vgl_polygon<double> poly;
    gt_istr >> poly;
    size_t ns = poly.num_sheets();
    gt_polys.push_back(poly);
  }
  gt_istr.close();
  vkm_binary_regions br;
 br.create_binary_image(site,polys);
 vil_save(br.bin_image(site), bin_image_path[site].c_str());
  //br.extract_region_boundaries("square_complex");
  //vil_save(br.region_image("square_complex"), region_img_path.c_str());

  br.create_binary_image("gt_"+site, gt_polys);
  vil_save(br.bin_image("gt_"+site), gt_bin_image_path[site].c_str());
  //  br.extract_region_boundaries("gt_square_complex");
  //vil_save(br.region_image("gt_square_complex"), gt_region_img_path.c_str());
  double r = 8.0;
  double max_iou = 0.0;
  double max_comp = 0.0, max_corr = 0.0;
  double tx_max = 0.0, ty_max = 0.0;
  for (double tx = -r; tx <= r; tx += 1)
    for (double ty = -r; ty <= r; ty += 1) {
      double iou, comp, corr;
      br.region_metrics(site, "gt_"+site, iou, corr, comp, tx, ty);
      std::cout << tx << ' ' << ty << ' ' << iou << std::endl;
      if(iou > max_iou){
        max_iou = iou;
        max_comp = comp;
        max_corr = corr;
        tx_max = tx;
        ty_max = ty;
      }
    }
  std::cout << "Footprint metrics: iou " << max_iou << " completeness " << max_comp << " correctness " << max_corr << std::endl;
  vil_image_view<unsigned char> overlay = br.generate_overlay_image(site, "gt_"+site, tx_max, ty_max);
  vil_save(overlay, overlay_path[site].c_str());
}
TESTMAIN(test_binary_regions);
