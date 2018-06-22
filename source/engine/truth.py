import os
import json
import argparse
import copy
import glob
import numpy as np

import openmesh
import vkm
from vxl import vnl

from osgeo import gdal, osr


# parse ground truth
def run_groundtruth(truthpath,outputpath):

  # check path existance
  if not os.path.isdir(truthpath):
    raise IOError('Input path not found <{}>'.format(truthpath))
  elif not os.path.isdir(outputpath):
    raise IOError('Output path not found <{}>'.format(outputpath))

  # locate region files
  filereg = [os.path.join(truthpath,f) for f in os.listdir(truthpath)
             if f.endswith(('.regions','.REGIONS'))]

  # initialize data
  keys = ('name','region','type','region_pc','dsm','valid')
  data = [dict.fromkeys(keys,None) for f in filereg]

  for file,item in zip(filereg,data):
    name = os.path.splitext(os.path.basename(file))[0]
    path = os.path.dirname(file)

    item['name'] = name
    item['region'] = file
    item['valid'] = False

    # type file
    f = os.path.join(path,name+'.types')
    if not os.path.isfile(f): continue
    item['type'] = f

    # region_pc file (optional)
    f = os.path.join(path,name+'.regions.xyz')
    if os.path.isfile(f): item['region_pc'] = f

    # lidar dsm
    f = os.path.join(path,'LiDAR.tif')
    if not os.path.isfile(f): continue
    item['dsm'] = f

    # glob for perimeter files
    perim = []
    for f in glob.glob(os.path.join(path,'*_perimeter.txt')):
      n = os.path.basename(f).replace('_perimeter.txt','')
      perim.append((n,f))
    item['perim'] = perim

    # set to valid
    item['valid'] = True

  # verbose discovered data
  print(json.dumps(data,indent=2))

  if not any((item['valid'] for item in data)):
    raise ValueError('No grouth truth regions discovered!')


  # process each item
  z_off = 232.111831665

  for item in data:
    if not item['valid']: continue

    # load DSM metadata
    dataset = gdal.Open(item['dsm'],gdal.GA_ReadOnly)
    band = dataset.GetRasterBand(1)
    meta = {
        'RasterXSize':  dataset.RasterXSize,
        'RasterYSize':  dataset.RasterYSize,
        'RasterCount':  dataset.RasterCount,
        'DataType':     band.DataType,
        'Projection':   dataset.GetProjection(),
        'GeoTransform': list(dataset.GetGeoTransform()),
        'NoDataValue':  band.GetNoDataValue(),
    }
    dataset = None

    # confirm expected DSM
    if meta['RasterCount'] != 1:
      raise ValueError('RasterCount is not 1')

    # srs = osr.SpatialReference(wkt=meta['Projection'])
    # if not srs.IsProjected or 'UTM' not in srs.GetAttrValue('projcs'):
    #   func_raise_error('Inputfile is not UTM')

    # ground truth pyvkm object
    gt = vkm.ground_truth(z_off)

    # load base information
    gt.load_ground_truth_img_regions(item['region'])
    gt.load_dsm_image(item['dsm'])
    gt.load_surface_types(item['type'])

    # compute image to x/y transformation
    if item['region_pc']:
      tmp = item['region_pc']
    else:
      G = meta['GeoTransform']
      img_to_xy = np.array([[G[1],G[2],0],[G[4],G[5],0],[0,0,1]],dtype=np.float)
      tmp = vnl.matrix_fixed_3x3(img_to_xy)

    gt.img_to_xy_trans(tmp)

    # processing
    gt.snap_image_region_vertices()
    # gt.convert_img_regions_to_meshes()
    gt.process_region_containment()
    gt.fit_region_planes()
    gt.construct_polygon_soup()
    gt.convert_to_meshes()

    # base output
    fileout = os.path.join(outputpath,item['name']+'_surfaces.obj')
    gt.write_processed_ground_truth(fileout)

    # load all site perimeters
    for p in item['perim']:
      gt.load_site_perimeter(p[0],p[1])

    # site perimiter ground planes
    fileout = os.path.join(outputpath,item['name']+'_ground_planes.txt')
    gt.write_ground_planes(fileout)

    # subregion output
    for p in item['perim']:
      gt.construct_extruded_gt_model(p[0])

      fileout = os.path.join(outputpath,p[0]+'_xy_polys.txt')
      gt.write_xy_polys(p[0],fileout)

      fileout = os.path.join(outputpath,p[0]+'_extruded.obj')
      gt.write_extruded_gt_model(p[0],fileout)



# command line function
def main(args=None):

  # setup input parser
  parser = argparse.ArgumentParser()
  parser.add_argument('-g', '--groundtruth', dest='truth',
      help='Ground truth directory', required=True)
  parser.add_argument('-o', '--output', dest='output',
      help='Output directory', required=True)

  # parse arguments
  args = parser.parse_args(args)
  print(args)

  # gather arguments
  kwargs = {
    'truthpath': args.truth,
    'outputpath': args.output,
  }

  # run function
  run_groundtruth(**kwargs)


if __name__ == "__main__":
  main()
