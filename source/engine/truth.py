import os
import json
import argparse
import copy
import glob
import numpy as np
import traceback

import vkm
from vxl import vnl, vpgl

from osgeo import gdal, osr


# parse ground truth
def run_groundtruth(truth_json,output_path):

  # check path existance
  if not os.path.isfile(truth_json):
    raise IOError('Ground truth JSON file not found <{}>'.format(truth_json))
  elif not os.path.isdir(output_path):
    raise IOError('Output path not found <{}>'.format(output_path))

  # helper function for supporting truth files
  truth_dir = os.path.dirname(truth_json)
  def truth_file(file):
    return os.path.join(truth_dir,file)


  # read truth data
  with open(truth_json,'r') as fid:
    data = json.load(fid)

  # confirm required fields
  keys = ['region','type','dsm']
  missing_keys = [k for k in keys if k not in data]
  if missing_keys:
      raise IOError('JSON input missing {}'.format(missing_keys))

  # check for truth files
  for key in keys:
    if not os.path.isfile(truth_file(data[key])):
      raise IOError('Cannot find "{}" file <{}>'.format(key,data[key]))

  # check perimeter files
  if "perim" in data:
    for item in data['perim']:
      if not os.path.isfile(truth_file(item['file'])):
        raise IOError('Cannot find "perim" file <{}>'.format(item['file']))


  # load DSM metadata
  dataset = gdal.Open(truth_file(data['dsm']),gdal.GA_ReadOnly)
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

  # calculate elevation range
  dsm = band.ReadAsArray()
  NDV = meta['NoDataValue']
  if NDV is not None:
    mn = np.amin(dsm[dsm!=NDV])
    mx = np.amax(dsm[dsm!=NDV])
  else:
    mn = np.amin(dsm)
    mx = np.amax(dsm)

  meta.update({
      'MinimumElevation': float(mn),
      'MaximumElevation': float(mx),
  })

  # clear DSM and GDAL dataaset
  del dsm
  del dataset

  # confirm expected DSM
  if meta['RasterCount'] != 1:
    raise IOError('Ground truth RasterCount must be 1')

  # confirm expected UTM
  srs = osr.SpatialReference(wkt=meta['Projection'])
  if not srs.IsProjected or 'UTM' not in srs.GetAttrValue('projcs'):
    raise IOError('Ground truth DSM must be UTM')

  # utm zone
  meta['UTMZone'] = srs.GetUTMZone()

  # report metadata
  print(json.dumps(meta,indent=2))


  # ground truth pyvkm object
  gt = vkm.ground_truth()

  # load base data (dsm, region, type)
  gt.load_dsm_image(truth_file(data['dsm']))
  gt.load_ground_truth_img_regions(truth_file(data['region']))
  gt.load_surface_types(truth_file(data['type']))


  # image to local x/y transformation
  # --geotransform locates the center of the upper left pixel at [.5,.5],
  #   while region files locate the center of the upper left pixel at [0,0].
  #   The matrix "B" accomodates this half pixel shift.
  # --this transform does not include the geotransform translational
  #   shift [G[0],G[3]], which is instead reported in "*_info.json"
  G = meta['GeoTransform']
  A = np.array( [ [G[1],G[2],0], [G[4],G[5],0], [0,0,1] ], dtype=np.float)
  B = np.array( [ [1,0,0.5], [0,1,0.5], [0,0,1]], dtype=np.float)
  img_to_xy = np.matmul(A,B)
  img_to_xy_vnl = vnl.matrix_fixed_3x3(img_to_xy)

  gt.img_to_xy_trans(img_to_xy_vnl)

  # z-offset = minimum valid DSM value
  gt.z_off(meta['MinimumElevation'])

  # UTM origin
  origin = [G[0],G[3],meta['MinimumElevation']]

  # VXL local vertical coordinate system (LVCS)
  # useful for VXL-based point projection
  utm_zone = abs(meta['UTMZone'])
  is_south = meta['UTMZone']<0

  utmobj = vpgl.utm()
  (lon,lat) = utmobj.utm2lonlat(origin[0],origin[1],utm_zone,is_south)
  elev = origin[2]

  lvcs = vpgl.lvcs(
    lat, lon, elev, vpgl.lvcs.cs_names.utm,
    vpgl.lvcs.AngUnits.DEG, vpgl.lvcs.LenUnits.METERS
  )


  # ground truth processing
  gt.snap_image_region_vertices()
  # gt.convert_img_regions_to_meshes()
  gt.process_region_containment()
  gt.fit_region_planes()
  gt.construct_polygon_soup()
  gt.convert_to_meshes()


  # output header
  output_header = os.path.join(output_path,
      os.path.splitext(os.path.basename(truth_json))[0])

  # write surfaces
  filesurf = output_header + '_surfaces.obj'
  gt.write_processed_ground_truth(filesurf)

  # ground truth info
  gtinfo = {
    "name": data["name"],
    "projection": meta['Projection'],
    "origin": origin,
    "lvcs": lvcs.writes().replace('\n',' ').strip(),
    "surfaces": os.path.basename(filesurf),
  }

  # write ground truth info to JSON
  fileinfo = output_header + "_truth.json"
  with open(fileinfo,'w') as fid:
    json.dump(gtinfo,fid,indent=2)


  # load all site perimeters
  for item in data['perim']:
    gt.load_site_perimeter(item['name'],truth_file(item['file']))

  # site perimeter ground planes
  fileground = output_header + '_ground_planes.txt'
  gt.write_ground_planes(fileground)

  # subregion output
  for item in data['perim']:
    name = item['name']

    gt.construct_extruded_gt_model(name)

    fileout = os.path.join(output_path, name+'_xy_polys.txt')
    gt.write_xy_polys(name,fileout)

    fileout = os.path.join(output_path, name+'_extruded.obj')
    gt.write_extruded_gt_model(name,fileout)




# command line function
def main(args=None):

  # setup input parser
  parser = argparse.ArgumentParser()
  parser.add_argument('-g', '--groundtruth', dest='truth_json',
      help='Ground truth JSON file', required=True)
  parser.add_argument('-o', '--output', dest='output_path',
      help='Output directory', required=True)

  # parse arguments as dict
  kwargs = vars(parser.parse_args(args))
  print(kwargs)

  # run function
  run_groundtruth(**kwargs)


if __name__ == "__main__":
  main()
