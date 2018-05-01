#ifndef vkm_obj_io_h
#define vkm_obj_io_h

#include <vector>
#include <map>
#include <string>
#include <limits>
#include <algorithm>
#include <iostream>

// VXL ------------------------------------
#include <vgl/vgl_polygon.h>
#include <vgl/vgl_point_2d.h>
//-----------------------------------------

// OpenMesh--------------------------------
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
//-----------------------------------------

#include "vkm_vrml_color.h"


class vkm_obj_io{

public:
  // write out a set of mesh objects into a single file. The input is a set of mesh groups.
  // The mesh groups are given a unique group id to enable selection and manipulation
  // by mesh visualization tools (.e.g, CloudCompare)
  // The entire scene can be assigned a material for visualization purposes
  // template over the type of mesh
  template< typename MeshT >
  static bool write_composite_obj_file(std::string const& obj_path,
                                       std::map<size_t, std::map<size_t,MeshT> >& scene,
                                       std::string mat_file = "", std::string mat = ""){
    std::ofstream ostr(obj_path.c_str());
    if(!ostr){
      std::cout << "Failed to open " << obj_path << " for writing obj file" << std::endl;
      return false;
    }
    if(mat_file != ""&& mat != ""){
      ostr << "mtllib " +  mat_file << std::endl;
      ostr << "usemtl "<< mat << std::endl;
    }
    size_t voff = 1;
    for(typename std::map<size_t, std::map<size_t, MeshT> >::iterator sit = scene.begin();
        sit != scene.end(); ++sit){
      size_t gid = sit->first;
      typename std::map<size_t, MeshT>& meshes = sit->second;
      for(typename std::map<size_t, MeshT>::iterator mit = meshes.begin();
          mit!= meshes.end(); ++mit){
        size_t n_verts = 0;
        std::stringstream ss;
        ss << "g " << gid << '_' << mit->first;
        ostr <<  ss.str() << std::endl;
        MeshT& pmesh = mit->second;
        for(typename MeshT::VertexIter vit = pmesh.vertices_begin();
            vit != pmesh.vertices_end(); ++vit){
          typename MeshT::Point p = pmesh.point(*vit);
          ostr << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
          n_verts++;
        }
        for (typename MeshT::FaceIter fit = pmesh.faces_begin();
             fit!=pmesh.faces_end(); ++fit){
          ostr << "f ";

          for (typename MeshT::FaceVertexIter fv_it = pmesh.fv_iter(*fit); fv_it.is_valid(); ++fv_it){
          typename OpenMesh::VertexHandle vh = *fv_it;
          unsigned int id = vh.idx();
          id += voff;
          ostr << id << ' ';
          }
          ostr << std::endl;
        }
        voff += n_verts;
      }
    }
    return true;
  }
  // using vector instead of map
  template< typename MeshT >
  static bool write_composite_obj_file(std::string const& obj_path,
                                       std::map<size_t, std::vector<MeshT> >& scene,
                                       std::string mat_file = "", std::string mat = ""){
    std::ofstream ostr(obj_path.c_str());
    if(!ostr){
      std::cout << "Failed to open " << obj_path << " for writing obj file" << std::endl;
      return false;
    }
    if(mat_file != ""&& mat != ""){
      ostr << "mtllib " +  mat_file << std::endl;
      ostr << "usemtl "<< mat << std::endl;
    }
    size_t voff = 1;
    for(typename std::map<size_t, std::vector< MeshT> >::iterator sit = scene.begin();
        sit != scene.end(); ++sit){
      size_t gid = sit->first;
      typename std::vector<MeshT>& meshes = sit->second;
      size_t sub_g = 0;
      for(typename std::vector<MeshT>::iterator mit = meshes.begin();
          mit!= meshes.end(); ++mit, sub_g++){
        size_t n_verts = 0;
        std::stringstream ss;
        ss << "g " << gid << '_' << sub_g;
        ostr <<  ss.str() << std::endl;
        MeshT& pmesh = (*mit);

        for(typename MeshT::VertexIter vit = pmesh.vertices_begin();
            vit != pmesh.vertices_end(); ++vit){
          typename MeshT::Point p = pmesh.point(*vit);
          ostr << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
          n_verts++;
        }
        pmesh.release_vertex_normals();

        for (typename MeshT::FaceIter fit = pmesh.faces_begin();
             fit!=pmesh.faces_end(); ++fit){
          ostr << "f ";

          for (typename MeshT::FaceVertexIter fv_it = pmesh.fv_iter(*fit); fv_it.is_valid(); ++fv_it){
          typename OpenMesh::VertexHandle vh = *fv_it;
          unsigned int id = vh.idx();
          id += voff;
          ostr << id << ' ';
          }
          ostr << std::endl;
        }
        voff += n_verts;
      }
    }
    return true;
  }

  template< typename MeshT >
  static bool read_composite_obj_file(std::string const& obj_path,
                                      std::map<size_t, std::map<size_t,MeshT> >& scene){

    std::ifstream istr(obj_path.c_str());
    if(!istr){
      std::cout << "Failed to open " << obj_path << " for writing obj file" << std::endl;
      return false;
    }
    std::string gs, gids;
    size_t gid = 0;
    size_t start_id = 1;
    size_t cum_id = 1;
    size_t minor_id = 0;
    std::map<size_t, size_t> vmap;
    std::map<size_t, MeshT>  cur_mesh_group;
    while(istr >> gs >> gids){
      if(gs != "g"){
        std::cout << "Parse failed " << gs << ' ' << gids << std::endl;
        return false;
      }
      std::string majs, mins;
      size_t temp;
      bool reading_major = true;
      for(std::string::iterator sit = gids.begin();
          sit != gids.end(); ++sit){
        if(reading_major&& *sit != '_'){
          majs.push_back(*sit);
		  continue;
		}else if(*sit == '_'){
          reading_major = false;
          continue;
        }
        mins.push_back(*sit);
      }

      std::stringstream maj_ss(majs), min_ss(mins);
      maj_ss >> temp;
      min_ss >> minor_id;
      MeshT mesh;
      // read the mesh vertices
      std::vector<typename MeshT::VertexHandle> in_vhandles;
      std::string type;
      while(true){
        double x, y, z;
        istr >> type;
        if(type != "v" && type != "vn" ){
          break;
        }
        istr >> x >> y >> z;
        if(type == "v"){
          typename MeshT::Point p(x, y, z);
          in_vhandles.push_back(mesh.add_vertex(p));
          cum_id++;
        }
      }//end while(true)

      // read the face
      if(type != "f"){
       std::cout << "Parse failed " << type << std::endl;
        return false;
      }
      std::vector<typename MeshT::VertexHandle> vhandles;
      size_t k = 0;
      for(size_t i = start_id; i<cum_id ; ++i, ++k){
        size_t vh;
        istr >> vh;
        vhandles.push_back(in_vhandles[k]);
      }
      start_id = cum_id;
      mesh.add_face(vhandles);
      // end of face read

      // add to scene[gid]
      if(temp == gid){
        cur_mesh_group[minor_id]=mesh;
      }else{//end of gid meshes
        if(cur_mesh_group.size()>0)
		 scene[gid] = cur_mesh_group;
        cur_mesh_group.clear();
        cur_mesh_group[minor_id] = mesh;
        gid = temp;
      }
    } // end of file

    // add last set of meshes to the scene
    if(cur_mesh_group.size())
      scene[gid] = cur_mesh_group;
	return true;

  }

  // write out a set of mesh objects into a single file. Each mesh object is assigned a group id
  // for manipulation by mesh visualization tools.
  //template over the type of mesh
  template< typename MeshT >
    static bool write_composite_obj_file(std::string const& obj_path,
                                       std::map<size_t, MeshT >& scene,
                                       std::string mat_file = "", std::string mat = ""){
    std::ofstream ostr(obj_path.c_str());
    if(!ostr){
      std::cout << "Failed to open " << obj_path << " for writing obj file" << std::endl;
      return false;
    }
    if(mat_file != ""&& mat != ""){
      ostr << "mtllib " +  mat_file << std::endl;
      ostr << "usemtl "<< mat << std::endl;
    }
    size_t voff = 1;
    for(typename std::map<size_t, MeshT >::iterator sit = scene.begin();
        sit != scene.end(); ++sit){
      size_t gid = sit->first;
      size_t n_verts = 0;
      std::stringstream ss;
      ss << "g " << gid;
      ostr <<  ss.str() << std::endl;
      MeshT& pmesh = sit->second;

	  for (typename MeshT::VertexIter vit = pmesh.vertices_begin();
		  vit != pmesh.vertices_end(); ++vit) {
		  typename MeshT::Point p = pmesh.point(*vit);
		  ostr << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
		  n_verts++;
	  }

      for (typename MeshT::FaceIter fit = pmesh.faces_begin();
           fit!=pmesh.faces_end(); ++fit){
        ostr << "f ";
        for (typename MeshT::FaceVertexIter fv_it = pmesh.fv_iter(*fit); fv_it.is_valid(); ++fv_it){
          OpenMesh::VertexHandle vh = *fv_it;
          unsigned int id = vh.idx();
          id += voff;
          ostr << id << ' ';
        }
        ostr << std::endl;
      }
      voff += n_verts;
    }
    return true;
  }

  // output composite objects with two parts, e.g. a) sliced cylinder and b) base extrusion
  // the two parts are assigned unique group ids, so for example the extrusion base can be removed
  // for visualization purposes.
  // template over the type of mesh
  template< typename MeshT >
  static bool write_composite_obj_file(std::string const& obj_path,
                                       const std::map<size_t, MeshT> & parts_a,
                                       const std::map<size_t, MeshT> & parts_b,
                                       std::string mat_file = "", std::string mat = ""){
    std::ofstream ostr(obj_path.c_str());
    if(!ostr){
      std::cout << "Failed to open " << obj_path << " for writing obj file" << std::endl;
      return false;
    }
    if(mat_file != ""&& mat != ""){
      ostr << "mtllib " +  mat_file << std::endl;
      ostr << "usemtl "<< mat << std::endl;
    }
    size_t voff = 1;
    typename std::map<size_t, MeshT >::const_iterator ait = parts_a.begin();
    typename std::map<size_t, MeshT >::const_iterator bit = parts_b.begin();
    for(; ait != parts_a.end()&& bit!= parts_b.end(); ++ait, ++bit){
      size_t n_verts = 0;
      std::stringstream ssa;
      ssa << "g " << ait->first << "_part_a";
      ostr <<  ssa.str() << std::endl;
      MeshT pmesh_a = ait->second;
      pmesh_a.request_vertex_normals();
      for(typename MeshT::VertexIter vit = pmesh_a.vertices_begin();
        vit != pmesh_a.vertices_end(); ++vit){
        typename MeshT::Point p = pmesh_a.point(*vit);
        ostr << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
        n_verts++;
        typename MeshT::Point n = pmesh_a.normal(*vit);
        ostr << "vn " << n[0] << ' ' << n[1] << ' ' << n[2] << std::endl;
        //std::cout<< "vn " << n[0] << ' ' << n[1] << ' ' << n[2] << std::endl;
      }
      pmesh_a.release_vertex_normals();
      for (typename MeshT::FaceIter fit = pmesh_a.faces_begin();
           fit!=pmesh_a.faces_end(); ++fit){
        ostr << "f ";
        for (typename MeshT::FaceVertexIter fv_it = pmesh_a.fv_iter(*fit); fv_it.is_valid(); ++fv_it){
          OpenMesh::VertexHandle vh = *fv_it;
          unsigned int id = vh.idx();
          id += voff;
          ostr << id << ' ';
        }
        ostr << std::endl;
      }
      voff += n_verts;
      n_verts = 0;
      std::stringstream ssb;
      ssb << "g " << bit->first << "_part_b";
      ostr <<  ssb.str() << std::endl;
      MeshT pmesh_b = bit->second;
      pmesh_b.request_vertex_normals();
      for(typename MeshT::VertexIter vit = pmesh_b.vertices_begin();
        vit != pmesh_b.vertices_end(); ++vit){
        typename MeshT::Point p = pmesh_b.point(*vit);
        ostr << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
        n_verts++;
        typename MeshT::Point n = pmesh_b.normal(*vit);
        ostr << "vn " << n[0] << ' ' << n[1] << ' ' << n[2] << std::endl;
        //std::cout<< "vn " << n[0] << ' ' << n[1] << ' ' << n[2] << std::endl;
      }
      pmesh_b.release_vertex_normals();
      for (typename MeshT::FaceIter fit = pmesh_b.faces_begin();
           fit!=pmesh_b.faces_end(); ++fit){
        ostr << "f ";
        for (typename MeshT::FaceVertexIter fv_it = pmesh_b.fv_iter(*fit); fv_it.is_valid(); ++fv_it){
          OpenMesh::VertexHandle vh = *fv_it;
          unsigned int id = vh.idx();
          id += voff;
          ostr << id << ' ';
        }
        ostr << std::endl;
      }
      voff += n_verts;
    }
    return true;
  }

  //writing composite obj from multiple vgl_polygon faces. This assumes 2D polygons and places them at a
  //set z value
  static bool write_composite_obj_file(std::string const& obj_path,
                                       std::map<size_t, vgl_polygon<double> > polyMap,
                                       std::string mat_file = "", std::string mat = "")
  {
    double heightVal = 1.0;
    std::ofstream ostr(obj_path.c_str());
    if(!ostr){
      std::cout << "Failed to open " << obj_path << " for writing obj file" << std::endl;
      return false;
    }
    size_t groupInd = 0;
    for(std::map<size_t,vgl_polygon<double> >::iterator it = polyMap.begin(); it != polyMap.end(); ++it)
    {
      ostr << "g group"<< groupInd <<std::endl;
      if(mat_file != ""&& mat != ""){
        ostr << "mtllib " +  mat_file << std::endl;
        ostr << "usemtl "<< mat << std::endl;
      }
      vgl_polygon<double> currPoly = it->second;
      size_t n_verts =  currPoly.num_vertices();
      if(n_verts==0){
        std::cout << "null alpha boundary" << std::endl;
        return false;
      }
      std::vector<vgl_point_2d<double> > verts = currPoly[0];
      for(int i = 0;i<static_cast<int>(verts.size());i++)
      {
        vgl_point_2d<double> pt = verts.at(i);
        ostr << "v " << pt.x()<< " " << pt.y() << " " << heightVal << std::endl;
      }
      ostr << "f ";
      int converted= static_cast<int>(verts.size());
      for(int i = 0;i<static_cast<int>(verts.size());i++)
      {
        ostr << (i - converted) << " ";
      }
      ostr<<std::endl;
      groupInd ++;
    }
    return true;
  }

  // map a property with real value on the range [0.0 1.0] to colors
  static void map_prop_to_color(double prop, double& r, double& g, double& b){
    if(prop<0.0) prop=0.0;
    else if(prop>1.0) prop=1.0;
    double ncolors = static_cast<double>(vkm_vrml_color::heatmap_custom_size);
    size_t color_index = static_cast<size_t>(prop*(ncolors-1));
    if(color_index>=vkm_vrml_color::heatmap_custom_size)
      color_index = vkm_vrml_color::heatmap_custom_size - 1;
    r = static_cast<double>(vkm_vrml_color::heatmap_custom[color_index][0])/255.0;
    g = static_cast<double>(vkm_vrml_color::heatmap_custom[color_index][1])/255.0;
    b = static_cast<double>(vkm_vrml_color::heatmap_custom[color_index][2])/255.0;
  }

  // write out a mesh with each vertex assigned a color according to the value of the designated property
  // if the requested property is not present on the mesh, the method returns false.
  template< typename MeshT >
  static bool write_obj_file_vertex_prop_as_color(std::string const& obj_path , MeshT& mesh, std::string const& prop_name){
    std::ofstream ostr(obj_path.c_str());
    if(!ostr){
      std::cout << "Failed to open " << obj_path << " for writing obj file" << std::endl;
      return false;
    }
    OpenMesh::VPropHandleT< double > prop;
    bool prop_valid = mesh.get_property_handle(prop, prop_name);
    if(!prop_valid){
      std::cout << "Property " << prop_name << " not valid on mesh" << std::endl;
      return false;
    }
    // normalize property values
    std::vector<double> properties;
    double min_v = std::numeric_limits<double>::max(), max_v = -min_v;
    for(typename MeshT::VertexIter vit = mesh.vertices_begin();
        vit != mesh.vertices_end(); ++vit){
      double pro = mesh.property(prop, *vit);
      if(pro<min_v) min_v = pro;
      if(pro>max_v) max_v = pro;
      properties.push_back(pro);
    }
    size_t np = properties.size();
    double low_val = min_v, high_val = max_v;
    if(properties.size()>1000){
    // compute tails to eliminate outliers
      std::vector<double> temp = properties;
      std::sort(temp.begin(), temp.end());
      double tail = 0.025;//the tail cutoff fraction
      size_t itail = static_cast<size_t>(np*tail);
      low_val = temp[itail]; high_val = temp[np-1-itail];
    }
    double norm = 1.0;
    double ndiff = 2.0*(high_val-low_val)/(high_val + fabs(low_val));
    if(ndiff>0.00)
      norm = 1.0/(high_val-low_val);
    for(std::vector<double>::iterator pit = properties.begin();
        pit != properties.end(); ++pit){
      if((*pit)>=low_val && (*pit)<=high_val)
        *pit = (*pit-low_val)*norm;
      if((*pit)<low_val) *pit = 0.0;
      if((*pit)>high_val) *pit = 1.0;
    }
    // write out the vertices
    size_t k = 0;
    mesh.request_vertex_normals();
    for(typename MeshT::VertexIter vit = mesh.vertices_begin();
        vit != mesh.vertices_end(); ++vit, k++){
      typename MeshT::Point p = mesh.point(*vit);
      double r, g, b;
      map_prop_to_color(properties[k], r, g, b);
      // note that the obj format supports vertex colors as three float values appended to the vertex coordinates
      ostr << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << r << ' ' << g << ' ' << b << std::endl;
      typename MeshT::Point n = mesh.normal(*vit);
      ostr << "vn " << n[0] << ' ' << n[1] << ' ' << n[2] << std::endl;
    }
    mesh.release_vertex_normals();
    for (typename MeshT::FaceIter fit = mesh.faces_begin();
         fit!=mesh.faces_end(); ++fit){
        ostr << "f ";
        for (typename MeshT::FaceVertexIter fv_it = mesh.fv_iter(*fit); fv_it.is_valid(); ++fv_it){
          OpenMesh::VertexHandle vh = *fv_it;
          unsigned int id = vh.idx()+1; // in obj format, vertex numbering for faces starts at 1
          ostr << id << ' ';
        }
        ostr << std::endl;
      }
    ostr.close();
    return true;
  }

  static void map_label_to_color(size_t label, double& r, double& g, double& b){
    if(label>=vkm_vrml_color::heatmap_custom_size)
      label = vkm_vrml_color::heatmap_custom_size - 1;
    r = static_cast<double>(vkm_vrml_color::heatmap_custom[label][0])/255.0;
    g = static_cast<double>(vkm_vrml_color::heatmap_custom[label][1])/255.0;
    b = static_cast<double>(vkm_vrml_color::heatmap_custom[label][2])/255.0;
  }

  template< typename MeshT >
  static bool write_obj_file_vertex_label_as_color(std::string const& obj_path , MeshT& mesh, std::string const& prop_name){
    std::ofstream ostr(obj_path.c_str());
    if(!ostr){
      std::cout << "Failed to open " << obj_path << " for writing obj file" << std::endl;
      return false;
    }
    OpenMesh::VPropHandleT< size_t > prop;
    bool prop_valid = mesh.get_property_handle(prop, prop_name);
    if(!prop_valid){
      std::cout << "Property " << prop_name << " not valid on mesh" << std::endl;
      return false;
    }
    std::set<size_t> unique_labels;
    for(typename MeshT::VertexIter vit = mesh.vertices_begin();
    vit != mesh.vertices_end(); ++vit){
      size_t label = mesh.property(prop, *vit);
      unique_labels.insert(label);
    }
    std::map<size_t, size_t> label_map;
    for(std::set<size_t>::iterator sit = unique_labels.begin();
        sit != unique_labels.end(); ++sit)
      label_map[*sit]= (2860486313*(*sit))%vkm_vrml_color::heatmap_custom_size;
    // write out the vertices
    for(typename MeshT::VertexIter vit = mesh.vertices_begin();
        vit != mesh.vertices_end(); ++vit){
      typename MeshT::Point p = mesh.point(*vit);
      size_t label = mesh.property(prop, *vit);
      size_t mapped_label = label_map[label];
      double r, g, b;
      map_label_to_color(mapped_label, r, g, b);
      // note that the obj format supports vertex colors as three float values appended to the vertex coordinates
      ostr << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << r << ' ' << g << ' ' << b << std::endl;
    }
    for (typename MeshT::FaceIter fit = mesh.faces_begin();
         fit!=mesh.faces_end(); ++fit){
      ostr << "f ";
      for (typename MeshT::FaceVertexIter fv_it = mesh.fv_iter(*fit); fv_it.is_valid(); ++fv_it){
        OpenMesh::VertexHandle vh = *fv_it;
        unsigned int id = vh.idx()+1; // in obj format, vertex numbering for faces starts at 1
        ostr << id << ' ';
      }
      ostr << std::endl;
    }
    ostr.close();
    return true;
  }

private:

  vkm_obj_io();

};

#endif //vkm_obj_io_h
