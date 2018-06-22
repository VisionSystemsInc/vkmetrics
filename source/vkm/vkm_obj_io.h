#ifndef vkm_obj_io_h
#define vkm_obj_io_h

#include <vector>
#include <map>
#include <string>
#include <limits>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>

// VXL ------------------------------------
#include <vgl/vgl_polygon.h>
#include <vgl/vgl_point_2d.h>
//-----------------------------------------

// OpenMesh--------------------------------
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
//-----------------------------------------

class vkm_obj_io{

public:


  /* read composite OBJ file (map of meshes)
   * map key = group identifier as string
   * map element = Openmesh object of type MeshT
   *
   * Note this function supports at most one group without an identifier.
   * Multiple groups without identifiers will result in an unspecified
   * operation.
  */
  template< typename MeshT >
  static bool read_composite_obj_file(
      std::string const& obj_path,
      std::map<std::string, MeshT>& groups)
  {
    // open file
    std::ifstream istr(obj_path.c_str());
    if(!istr){
      std::cerr << "Failed to open " << obj_path << " for writing obj file" << std::endl;
      return false;
    }

    std::string line, type, id = "unknown";
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
          groups[id] = mesh;

        // clear mesh for new values
        offset += mesh.n_vertices();
        mesh.clear();

        // save new identifier
        std::getline(iss,id);
        while (!id.empty() && id[0]==' ') {id.erase(0,1);}
        if (id.empty()) {id = "unknown";}

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
      groups[id] = mesh;

    return true;
  }


  /* read map of map of meshes
   * map 1st level key = major identifier
   * map 2nd level key = minor identifier
   * map element = Openmesh object of type MeshT
   *
   * group identifiers not of the form "<MAJOR>_<MINOR>", with <MAJOR>
   * and <MINOR> of type unsigned integer may result in unspecifed operation.
   *
   * optinally return any un-parsed property information for each group
   * for later custom processing.
  */
  template< typename MeshT >
  static bool read_composite_obj_file(
      std::string const& obj_path,
      std::map<size_t, std::map<size_t, MeshT> >& scene)
  {
    std::map<size_t, std::map<size_t, std::string> > properties;
    return read_composite_obj_file(obj_path,scene,properties);
  }


  template< typename MeshT >
  static bool read_composite_obj_file(
      std::string const& obj_path,
      std::map<size_t, std::map<size_t, MeshT> >& scene,
      std::map<size_t, std::map<size_t, std::string> >& properties)
  {
    std::map<std::string, MeshT> groups;
    if (!read_composite_obj_file(obj_path,groups))
      return false;

    std::istringstream iss_group, iss_id;
    std::string smajor, sminor, id, prop;
    size_t major, minor;

    // parse scene
    for (auto& g : groups) {

      // clear variables
      iss_group.clear(); iss_id.clear();
      smajor.clear(); sminor.clear(); id.clear(); prop.clear();
      major = 0; minor = 0;

      // split into ID and remaining property
      iss_group.str(g.first);
      iss_group >> id;

      std::getline(iss_group,prop);
      while (!prop.empty() && prop[0]==' ') {prop.erase(0,1);}

      std::cout << "Group ID: " << id << std::endl;

      // parse major/minor indices
      iss_id.str(id);
      std::getline(iss_id,smajor,'_');
      std::getline(iss_id,sminor,'_');

      if (smajor.length() > 0) {major = std::stoi(smajor);}
      if (sminor.length() > 0) {minor = std::stoi(sminor);}

      // add to scene
      scene[major][minor] = g.second;
      properties[major][minor] = prop;
    }

    return true;
  }


  /* write composite OBJ file (map of meshes)
   * map key = group identifier as string
   * map val = Openmesh object of type MeshT
   *
   * --The entire scene can be assigned a material for visualization purposes
   * --templated over the type of mesh
  */
  template< typename MeshT >
  static bool write_composite_obj_file(
      std::string const& obj_path,
      std::map<std::string, MeshT> const& scene,
      std::string const& mat_file = "", std::string const& mat = "")
  {
    // groups (map of MeshT pointers)
    std::map<std::string, const MeshT*> groups;
    for (const auto& item : scene)
      groups[item.first] = &(item.second);

    // write via pointer function
    return write_group_pointer_obj_file(obj_path,groups,mat_file,mat);
  }

  /* write map of meshes
   * map key = group identifier as size_t
  */
  template< typename MeshT >
  static bool write_composite_obj_file(
      std::string const& obj_path,
      std::map<size_t, MeshT >& scene,
      std::string const& mat_file = "", std::string const& mat = "")
  {
    // groups (map of MeshT pointers)
    std::stringstream ss;
    std::map<std::string, const MeshT*> groups;
    for (const auto& item : scene) {
      ss.str("");
      ss << item.first;
      groups[ss.str()] = &(item.second);
    }

    // write via pointer function
    return write_group_pointer_obj_file(obj_path,groups,mat_file,mat);
  }


  /* write map of map of meshes
   * map 1st level key = major identifier (size_t)
   * map 2nd level key = minor identifier (size_t)
   * map val = Openmesh object of type MeshT
   *
   * group identifiers will be written in the "<MAJOR>_<MINOR>"
  */
  template< typename MeshT >
  static bool write_composite_obj_file(
      std::string const& obj_path,
      std::map<size_t, std::map<size_t, MeshT> > const& scene,
      std::string const& mat_file = "", std::string const& mat = "")
  {
    // flatten structure (map of MeshT pointers)
    std::stringstream ss;
    std::map<std::string, const MeshT*> groups;
    for (const auto& major : scene) {
      for (const auto& minor : major.second ) {
        ss.str("");
        ss << major.first << "_" << minor.first;
        groups[ss.str()] = &(minor.second);
      }
    }

    // write via pointer function
    return write_group_pointer_obj_file(obj_path,groups,mat_file,mat);
  }

  /* write map of vector of meshes
   * map 1st level key = major identifier (size_t)
   * map 2nd level key = vector index [0:vec.length)
   * map val = Openmesh object of type MeshT
   *
   * group identifiers will be written in the "<MAJOR>_<INDEX>"
  */
  template< typename MeshT >
  static bool write_composite_obj_file(
      std::string const& obj_path,
      std::map<size_t, std::vector<MeshT> >& scene,
      std::string const& mat_file = "", std::string const& mat = "")
  {
    // flatten structure (map of MeshT pointers)
    std::stringstream ss;
    std::map<std::string, const MeshT*> groups;
    for (const auto& major : scene) { // map of vectors
      for (int idx=0; idx < major.second.size(); idx++) { // vector item
        ss.str("");
        ss << major.first << "_" << idx;
        groups[ss.str()] = &(major.second[idx]);
      }
    }

    // write via pointer function
    return write_group_pointer_obj_file(obj_path,groups,mat_file,mat);
  }


  /* Main function to write OBJ file
   * accepts map of pointers to Open mesh objects
   *
   * map key = group string identifier
   * map element = pointer to Openmesh object of type MeshT
   *
   * The entire scene can be assigned a material for visualization purposes
   *
   * templated over the type of mesh
  */
  template< typename MeshT >
  static bool write_group_pointer_obj_file(
      std::string const& obj_path,
      std::map<std::string, MeshT*> const& groups,
      std::string const& mat_file = "", std::string const& mat = "")
  {
    // open file
    std::ofstream ostr(obj_path.c_str());
    if(!ostr){
      std::cout << "Failed to open " << obj_path << " for writing obj file" << std::endl;
      return false;
    }

    // add material stuff
    if(mat_file != "" && mat != ""){
      ostr << "mtllib " + mat_file << std::endl;
      ostr << "usemtl "<< mat << std::endl;
    }

    // write groups
    size_t voff = 1;

    for (const auto& g : groups) {
      ostr << "g " << g.first << std::endl;
      // std::cout << g.first << " has " << g.second->n_vertices() << " vertices" << std::endl;

      for (const auto& vh : g.second->vertices())
        ostr << "v " << g.second->point(vh) << std::endl;

      for (const auto& fh : g.second->faces()) {
        ostr << "f";
        for (const auto& vh : g.second->fv_range(fh))
          ostr << " " << vh.idx()+voff;
        ostr << std::endl;
      }

      voff += g.second->n_vertices();
    } // end mesh loop

    return true;
  }


private:

  vkm_obj_io();

};

#endif //vkm_obj_io_h
