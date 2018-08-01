import os
import json
import argparse
import copy
import glob
from collections import OrderedDict

import vkm
from vxl import vpgl


# POINT CONVERSION ROUTINES

# LVCS point conversion:
# Convert point through local vertical coordinate system (LVCS)
# this function will take a local point in a source lvcs,
# convert the point to global WGS84 space (lat,lon,elev),
# then convert to the destination lvcs

# convert coordinates through WGS84 space
transform_point_convert = (vpgl.lvcs.cs_names.wgs84,vpgl.lvcs.AngUnits.DEG,
                           vpgl.lvcs.LenUnits.METERS)

# helper function
def print_tuple(tup,fmt='{:6.2f}'):
  return ','.join([fmt.format(e) for e in tup])

# conversion function
def convert_point_lvcs(pt_src,lvcs_src,lvcs_dst,verbose=False):

  # local_src to global
  pt_global = lvcs_src.local_to_global(*pt_src,*transform_point_convert)

  # global to local_dst
  pt_dst = lvcs_dst.global_to_local(*pt_global,*transform_point_convert)

  # verbose output
  if verbose:
    print('PT_SRC ({}) -> PT_WGS84 ({}) -> PT_DST ({})'.format(
        print_tuple(pt_src),
        print_tuple(pt_global),
        print_tuple(pt_dst))
    )

  # cleanup
  return pt_dst



# parse ground truth
def run_metrics(truth_json,input_json,output_path):

  # check path existance
  if not os.path.isfile(truth_json):
    raise IOError('Ground truth JSON file not found <{}>'.format(truth_json))
  elif not os.path.isfile(input_json):
    raise IOError('Input JSON file not found <{}>'.format(input_json))
  elif not os.path.isdir(output_path):
    raise IOError('Output directory not found <{}>'.format(output_path))

  # helper function for supporting truth files
  truth_dir = os.path.dirname(truth_json)
  def truth_file(file):
    return os.path.join(truth_dir,file)

  # helper function for supporting input files
  input_dir = os.path.dirname(input_json)
  def input_file(file):
    return os.path.join(input_dir,file)


  # ground truth
  with open(truth_json,'r') as fid:
    truth_data = json.load(fid)

  # confirm required fields
  keys = ['projection','origin','surfaces','lvcs']
  missing_keys = [k for k in keys if k not in truth_data]
  if missing_keys:
      raise IOError('Ground truth JSON missing {}'.format(missing_keys))

  # check for truth support files
  keys = ['surfaces']
  for key in keys:
    if not os.path.isfile(truth_file(truth_data[key])):
      raise IOError('GCannot find "{}" file <{}>'.format(key,truth_data[key]))


  # input data
  with open(input_json,'r') as fid:
    input_data = json.load(fid)

  # confirm required fields
  keys = ['coordinate_system','models']
  missing_keys = [k for k in keys if k not in input_data]
  if missing_keys:
      raise IOError('Input JSON missing {}'.format(missing_keys))

  # check coordiante system
  if 'type' not in input_data['coordinate_system']:
    raise IOError('Input JSON does not specify coordinate system "type"')

  # check OBJ files
  for item in input_data['models']:
    if not os.path.isfile(input_file(item['file'])):
      raise IOError('Cannot find "models" file <{}>'.format(item['file']))


  # verbose report
  print('\nGROUND TRUTH:')
  print(json.dumps(truth_data,indent=2))

  print('\nINPUT:')
  print(json.dumps(input_data,indent=2))


  # data structure for input data in ground truth coordinate system
  converted_data = copy.deepcopy(input_data)
  def converted_file(file):
    return os.path.join(output_path,file)

  converted_data['name'] = converted_data.get('name','') + '_converted'

  for item in converted_data['models']:
    n,e = os.path.splitext(os.path.basename(item['file']))
    item['file'] = "{}_converted{}".format(n,e)


  # input->truth point conversion routine
  cs_type = input_data['coordinate_system']['type'].lower()

  if cs_type == 'lvcs':
    input_lvcs = vpgl.lvcs()
    input_lvcs.reads(input_data['coordinate_system']['parameters'])

    truth_lvcs = vpgl.lvcs()
    truth_lvcs.reads(truth_data['lvcs'])

    def convert_point(input_pt):
      return convert_point_lvcs(input_pt,input_lvcs,truth_lvcs,False)

    converted_data['coordinate_system'] = {
        "type": "LVCS",
        "parameters": truth_lvcs.writes().replace('\n',' ').strip(),
    }

  # define other conversion types...
  # elif cs_type == 'other':
  #   ....

  # unhandled coordinate system
  else:
    raise IOError('Unhandled input coordinate system "{}"'.format(cs_type))


  # write converted JSON
  fileout = os.path.splitext(os.path.basename(input_json))[0] + "_converted.json"
  fileout = converted_file(fileout)
  with open(fileout,'w') as fid:
    json.dump(converted_data,fid,indent=2)


  # convert model files to ground truth coordinate system
  for itemin,itemout in zip(input_data['models'],converted_data['models']):

    filein = input_file(itemin['file'])
    fileout = converted_file(itemout['file'])
    print("Converting: {} -> {}".format(filein,fileout))

    # convert OBJ vertices
    # manually read the file, as other readers (e.g., OpenMesh)
    # may change the file structure or remove comments/group information
    if filein.endswith(('.obj','.OBJ')):
      with open(filein,'r') as fin, open(fileout,'w') as fout:

        for cnt,line in enumerate(fin):

          # convert vertices
          if line[0] == 'v':
            vals = line.split()
            pt_src = [float(v) for v in vals[1:]]
            pt_dst = convert_point(pt_src)

            pt_dst_str = ["{:.6g}".format(v) for v in pt_dst]
            buf = "v {}\n".format(' '.join(pt_dst_str))
            fout.write(buf)

          # write all other lines directly
          else:
            fout.write(line)

          # report conversion progress
          if cnt>0 and cnt%1e4 == 0:
            print('  line {}...'.format(cnt))

      print('  Conversion complete ({} lines)'.format(cnt))


    # define other conversion types...
    # elif filein.endswith(('.ext','.EXT')):
    #   ....

    # unhandled file type
    else:
      raise IOError('Unhandled input file extension "{}"'.format(e))



  # metrics for each input file
  for item in converted_data['models']:
    print(item)

    # initialize object
    met = vkm.metrics()

    # load models
    met.load_ground_truth_model(truth_file(truth_data['surfaces']))
    met.load_simply_connected_test_model(
        converted_file(item['file']),item['name'])

    # prepare models
    # met.delete_isolated_vertices() # not typically necessary
    met.construct_xy_regions()
    # met.set_translation(tr[0], tr[1]) # TODO - REGISTRATION
    met.translate_test_model_xy()
    met.match_xy_regions()

    # compute scores
    met.compute_best_match_2d_score_stats()
    met.find_z_offset()
    met.compute_best_match_3d_score_stats()

    # report scores as printed array
    met.print_score_array()

    # save scores as json
    scores = json.loads(met.json(), object_pairs_hook=OrderedDict)
    fileout = converted_file('metrics.json')
    with open(fileout,'w') as fid:
      json.dump(scores,fid,indent=2)

    # met.load_ground_planes(datagt['ground'])
    # fileout = os.path.join(output_path,name+"_registered.obj")
    # met.save_transformed_regions_as_meshes(fileout, name)

    # cleanup
    del met


# command line function
def main(args=None):

  # setup input parser
  parser = argparse.ArgumentParser()
  parser.add_argument('-g', '--groundtruth', dest='truth_json',
      help='Ground truth model JSON', required=True)
  parser.add_argument('-i', '--input', dest='input_json',
      help='Input model JSON', required=True)
  parser.add_argument('-o', '--output', dest='output_path',
      help='Output directory', required=True)

  # parse arguments as dict
  kwargs = vars(parser.parse_args(args))
  print(kwargs)

  # run function
  run_metrics(**kwargs)


if __name__ == "__main__":
  main()
