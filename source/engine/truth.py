import os
import json
import argparse
import copy
import glob

import vkm


# parse ground truth
def run_groundtruth(inputpath,outputpath):

  # check path existance
  if not os.path.isdir(inputpath):
    raise IOError('Input path not found <{}>'.format(inputpath))
  elif not os.path.isdir(outputpath):
    raise IOError('Output path not found <{}>'.format(outputpath))

  # locate file(s) of interest
  keys = ('region','type','region_pc','dsm')
  tmp = dict.fromkeys(keys)
  tmp['valid'] = False

  data = {}
  for root, dirs, files in os.walk(inputpath):
    for file in files:

      if file.endswith(('.regions','.REGIONS')):
        key = 'region'
      elif file.endswith(('.types','.TYPES')):
        key = 'type'
      else:
        continue

      id = os.path.splitext(file)[0]
      if id not in data: data[id] = copy.deepcopy(tmp)
      data[id][key] = os.path.join(root,file)

  # process each item
  for id in data:
    item = data[id]
    path = os.path.dirname(item['region'])

    # region file
    f = item['region']
    if not f or not os.path.isfile(f): continue

    # region types
    f = item['type']
    if not f or not os.path.isfile(f): continue

    # region_pc file
    f = item['region'] + '.xyz'
    if not os.path.isfile(f): continue
    item['region_pc'] = f

    # lidar dsm
    f = os.path.join(path,'LiDAR.tif')
    if not os.path.isfile(f): continue
    item['dsm'] = f

    # glob for perimeter files
    perim = []
    for file in glob.glob(os.path.join(path,'*_perimeter.txt')):
      name = os.path.basename(file).replace('_perimeter.txt','')
      perim.append((name,file))
    item['perim'] = perim

    # set to valid
    item['valid'] = True

  # verbose discovered data
  print(json.dumps(data,indent=2))

  # process each item
  for id in data:
    item = data[id]
    if not item['valid']: continue

    z_off = 232.111831665

    gt = vkm.ground_truth(z_off)
    gt.load_ground_truth_img_regions(item['region'])
    gt.load_surface_types(item['type'])
    gt.load_dem_image(item['dsm'])

    gt.load_ground_truth_pc_regions(item['region_pc'])
    gt.compute_img_to_xy_trans()

    gt.snap_image_region_vertices()
    gt.convert_img_regions_to_meshes()
    gt.process_region_containment()
    gt.fit_region_planes()
    gt.construct_polygon_soup()
    gt.convert_to_meshes()

    fileout = os.path.join(outputpath,id+'_surfaces.obj')
    gt.write_processed_ground_truth(fileout)

    fileout = os.path.join(outputpath,id+'_ground_planes.txt')
    gt.write_ground_planes(fileout)

    # process subregions
    for p in item['perim']:
      gt.load_site_perimeter(p[0],p[1])

      fileout = os.path.join(outputpath,p[0]+'_xy_polys.txt')
      gt.write_xy_polys(p[0],fileout)

      gt.construct_extruded_gt_model(p[0])

      fileout = os.path.join(outputpath,p[0]+'_extruded.obj')
      gt.write_extruded_gt_model(p[0],fileout)



# command line function
def main(args=None):

  # setup input parser
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', dest='input',
      help='Input path', required=True)
  parser.add_argument('-o', '--output', dest='output',
      help='Output path', required=True)

  # parse arguments
  args = parser.parse_args(args)
  print(args)

  # gather arguments
  kwargs = {
    'inputpath': args.input,
    'outputpath': args.output,
  }

  # run function
  run_groundtruth(**kwargs)


if __name__ == "__main__":
  main()
