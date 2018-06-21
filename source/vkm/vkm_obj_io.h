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

class vkm_obj_io{

public:

  // Read a composite OBJ file, saving to a map of maps of meshes
  template< typename MeshT >
  static bool read_composite_obj_file(std::string const& obj_path,
                                      std::map<size_t, std::map<size_t,MeshT> >& scene){

    // open file
    std::ifstream istr(obj_path.c_str());
    if(!istr){
      std::cout << "Failed to open " << obj_path << " for writing obj file" << std::endl;
      return false;
    }

    std::string line;
    std::string type;

    std::string temp;
    bool reading_major;
    std::string smajor, sminor;

    size_t major=0, minor=0;

    MeshT mesh;

    double x,y,z;
    size_t i, offset=0;
    std::vector<typename MeshT::VertexHandle> vhandles;

    // loop through all available lines
    while (std::getline(istr,line)) {

      std::istringstream iss(line);
      iss >> type;

      // parse group identifier
      if (type == "g") {

        // save previous mesh
        if (mesh.n_faces() > 0)
          scene[major][minor] = mesh;

        // clear mesh for new values
        offset += mesh.n_vertices();
        mesh.clear();

        // update group identifiers
        smajor.clear(); sminor.clear();
        iss >> temp;

        reading_major = true;
        for(std::string::iterator sit = temp.begin(); sit != temp.end(); ++sit)
        {
          if (*sit == '_') {
            reading_major = false;
          } else if (reading_major) {
            smajor.push_back(*sit);
          } else {
            sminor.push_back(*sit);
          }
        }

        major = 0;
        if (smajor.length() > 0)
          major = std::stoi(smajor);

        minor = 0;
        if (sminor.length() > 0)
          minor = std::stoi(sminor);

      } //end type=='g'

      // parse vertex
      else if (type == "v") {
        iss >> x >> y >> z;
        mesh.add_vertex(typename MeshT::Point(x,y,z));

      } //end type=='v'

      // parse face
      else if (type == "f") {
        vhandles.clear();
        while (iss >> i) {
          vhandles.push_back(mesh.vertex_handle(i-1-offset));
        }
        mesh.add_face(vhandles);

      } // end type=='f'

    } // end line loop

    // save last mesh
    if (mesh.n_faces() > 0)
      scene[major][minor] = mesh;

    return true;

  }


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

private:

  vkm_obj_io();

};

#endif //vkm_obj_io_h
