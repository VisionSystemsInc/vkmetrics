import os
import json
import argparse
import copy
import glob

import openmesh
import vkm


# parse ground truth
def run_metrics(inputpath,truthpath,outputpath):

  # check path existance
  if not os.path.isdir(inputpath):
    raise IOError('Input directory not found <{}>'.format(inputpath))
  elif not os.path.isdir(truthpath):
    raise IOError('Ground truth directory not found <{}>'.format(truthpath))
  elif not os.path.isdir(outputpath):
    raise IOError('Output directory not found <{}>'.format(outputpath))


  # ground truth surface files
  filegt = [os.path.join(truthpath,f) for f in os.listdir(truthpath)
            if f.endswith('_surfaces.obj')]
  if len(filegt) == 0:
    raise IOError('Missing ground truth "*_surfaces.obj"')

  # ground truth data
  keys = ('name','surface','ground','valid')
  datagt = [dict.fromkeys(keys,None) for f in filegt]

  for file,item in zip(filegt,datagt):
    name = os.path.splitext(os.path.basename(file))[0]
    name = name.rpartition('_surfaces')[0]

    path = os.path.dirname(file)

    item['name'] = name
    item['surface'] = file
    item['valid'] = False

    # ground planes
    f = os.path.join(path,name+'_ground_planes.txt')
    if not os.path.isfile(f): continue
    item['ground'] = f

    # set to valid
    item['valid'] = True

  # handle only one ground truth file at this time
  datagt = [item for item in datagt if item['valid']]
  if len(datagt) > 1:
    raise IOError('Too many ground truth OBJ files')
  datagt = datagt[0]


  # DEBUGGING: READ OBJ FILE
  fileout = os.path.join(outputpath,'groups.obj')
  groups = vkm.obj_io.read_group(datagt['surface'])
  vkm.obj_io.write(fileout, groups)

  fileout = os.path.join(outputpath,'meshes.obj')
  meshes = vkm.obj_io.read_composite(datagt['surface'])
  vkm.obj_io.write(fileout, meshes)

  raise ValueError('DEBUGGING')


  # locate input files
  fileinput = [os.path.join(inputpath,f) for f in os.listdir(inputpath)
              if f.endswith(('.obj','.OBJ'))]
  if len(fileinput) == 0:
    raise IOError('Missing input "*.obj"')



  # verbose report
  print('\nGROUND TRUTH:')
  print(json.dumps(datagt,indent=2))

  print('\nINPUT FILES:')
  for f in fileinput: print(f)


  # metrics for each input file
  z_off = 232.111831665;
  z_gnd_elev = 17.9843;

  for file in fileinput:

    met = vkm.metrics(z_off,z_gnd_elev)
    # met.set_translation(tr[0], tr[1])

    met.load_gnd_truth_model(datagt['surface'])
    met.load_ground_planes(datagt['ground']);

    met.load_simply_connected_test_model(file);
    met.construct_xy_regions();
    met.translate_test_model_xy();
    met.match_xy_regions();

    met.compute_best_match_2d_score_stats();
    met.find_z_offset();
    # met.save_transformed_regions_as_meshes(site, registered_model_path[site]);
    met.compute_best_match_3d_score_stats();
    met.print_score_array();






 #  std::string site = "square_complex";
 #  //std::string site = "small_square_complex";
 #  //std::string site = "yasu_house";
 #  //std::string site = "firehouse";
 #  //std::string site = "directorate";

 #  std::string d2_dir = "D:/core3d_experiments/phase1a_D2/";
 #  std::string gt_dir = "D:/core3d_experiments/phase1a_D2/GTruth_data_v2/";

 #  std::map<std::string, std::string> xy_poly_path, gt_xy_poly_path, bin_image_path, gt_bin_image_path,  overlay_path;

 #  xy_poly_path["square_complex"] = d2_dir + "square_complex/xy_polys.txt";
 #  xy_poly_path["small_square_complex"] = d2_dir + "small_square_complex/xy_polys.txt";
 #  xy_poly_path["yasu_house"] = d2_dir + "yasu_house/xy_polys.txt";
 #  xy_poly_path["firehouse"] = d2_dir + "firehouse/xy_polys.txt";
 #  xy_poly_path["directorate"] = d2_dir + "directorate/xy_polys.txt";
  
 #  gt_xy_poly_path["square_complex"]= gt_dir + "square_complex_xy_polys.txt";
 #  gt_xy_poly_path["small_square_complex"]= gt_dir + "small_square_complex_xy_polys.txt";
 #  gt_xy_poly_path["yasu_house"]= gt_dir + "yasu_house_xy_polys.txt";
 #  gt_xy_poly_path["firehouse"]= gt_dir + "firehouse_xy_polys.txt";
 #  gt_xy_poly_path["directorate"]= gt_dir + "directorate_xy_polys.txt";
  
 #  bin_image_path["square_complex"] = d2_dir + "square_complex/footprint_bin_image.tif";
 #  bin_image_path["small_square_complex"] = d2_dir + "small_square_complex/footprint_bin_image.tif";
 #  bin_image_path["yasu_house"] = d2_dir + "yasu_house/footprint_bin_image.tif";
 #  bin_image_path["firehouse"] = d2_dir + "firehouse/footprint_bin_image.tif";
 #  bin_image_path["directorate"] = d2_dir + "directorate/footprint_bin_image.tif";

 #  gt_bin_image_path["square_complex"] = gt_dir + "square_complex_footprint_bin_image.tif";
 #  gt_bin_image_path["small_square_complex"] = gt_dir + "small_square_complex_footprint_bin_image.tif";
 #  gt_bin_image_path["yasu_house"] = gt_dir + "yasu_house_footprint_bin_image.tif";
 #  gt_bin_image_path["firehouse"] = gt_dir + "firehouse_footprint_bin_image.tif";
 #  gt_bin_image_path["directorate"] = gt_dir + "directorate_footprint_bin_image.tif";

 #  overlay_path["square_complex"] = gt_dir + "square_complex_overlay_image.tif";
 #  overlay_path["small_square_complex"] = gt_dir + "small_square_complex_overlay_image.tif";
 #  overlay_path["yasu_house"] = gt_dir + "yasu_house_overlay_image.tif";
 #  overlay_path["firehouse"] = gt_dir + "firehouse_overlay_image.tif";
 #  overlay_path["directorate"] = gt_dir + "directorate_overlay_image.tif";
  
 #  std::vector<vgl_polygon<double> > polys, gt_polys;
 #  std::ifstream istr(xy_poly_path[site].c_str());
 #  if(!istr)
 #    return;
 #  size_t n;
 #  istr >> n;
 #  for(size_t i = 0;i<n; ++i){
 #    vgl_polygon<double> poly;
 #    istr >> poly;
 #    size_t ns = poly.num_sheets();
 #    polys.push_back(poly);
 #  }
 #  istr.close();
 #  std::ifstream gt_istr(gt_xy_poly_path[site].c_str());
 #  if(!gt_istr)
 #    return;
 #  gt_istr >> n;
 #  for(size_t i = 0;i<n; ++i){
 #    vgl_polygon<double> poly;
 #    gt_istr >> poly;
 #    size_t ns = poly.num_sheets();
 #    gt_polys.push_back(poly);
 #  }
 #  gt_istr.close();
 #  binary_regions br;
 # br.create_binary_image(site,polys);
 # vil_save(br.bin_image(site), bin_image_path[site].c_str());
 #  //br.extract_region_boundaries("square_complex");
 #  //vil_save(br.region_image("square_complex"), region_img_path.c_str());

 #  br.create_binary_image("gt_"+site, gt_polys);
 #  vil_save(br.bin_image("gt_"+site), gt_bin_image_path[site].c_str());
 #  //  br.extract_region_boundaries("gt_square_complex");
 #  //vil_save(br.region_image("gt_square_complex"), gt_region_img_path.c_str());
 #  double r = 8.0;
 #  double max_iou = 0.0;
 #  double max_comp = 0.0, max_corr = 0.0;
 #  double tx_max = 0.0, ty_max = 0.0;
 #  for (double tx = -r; tx <= r; tx += 1)
 #    for (double ty = -r; ty <= r; ty += 1) {
 #      double iou, comp, corr;
 #      br.region_metrics(site, "gt_"+site, iou, corr, comp, tx, ty);
 #      std::cout << tx << ' ' << ty << ' ' << iou << std::endl;
 #      if(iou > max_iou){
 #        max_iou = iou;
 #        max_comp = comp;
 #        max_corr = corr;
 #        tx_max = tx;
 #        ty_max = ty;
 #      }
 #    }
 #  std::cout << "Footprint metrics: iou " << max_iou << " completeness " << max_comp << " correctness " << max_corr << std::endl;
 #  vil_image_view<unsigned char> overlay = br.generate_overlay_image(site, "gt_"+site, tx_max, ty_max);
 #  vil_save(overlay, overlay_path[site].c_str());



#   //std::string site = "yasu_house";
#   //std::string site = "firehouse";
#  //std::string site = "directorate";
#   //std::string site = "square_complex";
#   //std::string site = "small_square_complex";
#   //std::string site = "multi_storey";
#   std::string site = "flat_roof_complex";

#   std::string gt_dir ="D:/core3d_experiments/phase1a_D2/GTruth_data_v2/";
#   std::string model_dir ="D:/core3d_experiments/phase1a_D2/";
#   std::map<std::string, std::string> test_model_path, registered_model_path;
#   test_model_path["yasu_house"] = model_dir + "yasu_house/alpha_meshes/alpha_hull_planar_meshes.obj";
#   test_model_path["firehouse"] = model_dir + "firehouse/alpha_meshes/alpha_hull_planar_meshes.obj";
#   test_model_path["directorate"] = model_dir + "directorate/alpha_meshes/alpha_hull_planar_meshes.obj";
#   test_model_path["square_complex"] = model_dir + "square_complex/alpha_meshes/alpha_hull_planar_meshes.obj";
#   test_model_path["small_square_complex"] = model_dir + "small_square_complex/alpha_meshes/alpha_hull_planar_meshes.obj";
#   test_model_path["multi_storey"] = model_dir + "multi_storey/alpha_meshes/alpha_hull_planar_meshes.obj";
#   test_model_path["flat_roof_complex"] = model_dir + "flat_roof_complex/alpha_meshes/alpha_hull_planar_meshes.obj";

#   registered_model_path["yasu_house"] = model_dir + "yasu_house/registered_model.obj";
#   registered_model_path["firehouse"] = model_dir + "firehouse/registered_model.obj";
#   registered_model_path["directorate"] = model_dir + "directorate/registered_model.obj";
#   registered_model_path["square_complex"] = model_dir + "square_complex/registered_model.obj";
#   registered_model_path["small_square_complex"] = model_dir + "small_square_complex/registered_model.obj";
#   registered_model_path["multi_storey"] = model_dir + "multi_storey/registered_model.obj";
#   registered_model_path["flat_roof_complex"] = model_dir + "flat_roof_complex/registered_model.obj";

#   std::string gt_surf_path = gt_dir + "CORE3D-Phase1a.kw18_surfaces.obj";
#   std::string ground_plane_path = gt_dir + "site_ground_planes.txt";
#   // translation found by basics/binary_regions, align binary images 
#   // see basics/tests/test_binary_regions
#   std::map<std::string, std::pair<double, double> > tr;
#   // more robust that intersection/union of footprint polygons
#   tr["square_complex"] = std::pair<double, double>(2.0, -3.0);
#   tr["yasu_house"] = std::pair<double, double>(2.0, -4.0);
#   tr["firehouse"] = std::pair<double, double>(3.0, 1.0);
#   tr["small_square_complex"] = std::pair<double, double>(2.0, 0.0);
#   tr["multi_storey"] = std::pair<double, double>(3.0, 1.0);

#   // offset due to lvcs
#   double z_off =   232.111831665;

#   // should be replaced by modl_ground_plane
#   double z_gnd_elev = 17.9843;

#   mesh_utils_metrics met(z_off, z_gnd_elev);
#   met.set_translation(tr[site].first, tr[site].second); 
#   met.load_gnd_truth_model(gt_surf_path);
#   met.load_ground_planes(ground_plane_path);
#   met.load_simply_connected_test_model(test_model_path[site]);
#   met.construct_xy_regions();
#   met.translate_test_model_xy();
#   met.match_xy_regions();

#   met.compute_best_match_2d_score_stats();
#   met.find_z_offset();
#   met.save_transformed_regions_as_meshes(site, registered_model_path[site]);
#   met.compute_best_match_3d_score_stats();
#   met.print_score_array();
#   //std::cout << "======JSON ENCODING======\n"<< met.json() << std::endl;







# command line function
def main(args=None):

  # setup input parser
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', dest='input',
      help='Input model directory', required=True)
  parser.add_argument('-g', '--groundtruth', dest='truth',
      help='Ground truth model directory', required=True)
  parser.add_argument('-o', '--output', dest='output',
      help='Output directory', required=True)

  # parse arguments
  args = parser.parse_args(args)
  print(args)

  # gather arguments
  kwargs = {
    'inputpath': args.input,
    'truthpath': args.truth,
    'outputpath': args.output,
  }

  # run function
  run_metrics(**kwargs)


if __name__ == "__main__":
  main()
