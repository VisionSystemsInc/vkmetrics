import os
import json
import argparse
import copy
import glob
from collections import OrderedDict

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
    name = os.path.splitext(os.path.basename(file))[0]

    met = vkm.metrics(z_off,z_gnd_elev)
    # met.set_translation(tr[0], tr[1]) # TODO - REGISTRATION

    met.load_ground_truth_model(datagt['surface'])
    met.load_simply_connected_test_model(file,name);
    met.construct_xy_regions();
    met.translate_test_model_xy();
    met.match_xy_regions();

    met.compute_best_match_2d_score_stats();
    met.find_z_offset();

    met.compute_best_match_3d_score_stats();

    scores = json.loads(met.json(), object_pairs_hook=OrderedDict)
    met.print_score_array();

    met.load_ground_planes(datagt['ground']);
    fileout = os.path.join(outputpath,name+"_registered.obj")
    met.save_transformed_regions_as_meshes(fileout, name);





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
