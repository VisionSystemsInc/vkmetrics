#include "pyvkm_obj_io.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <vkm_obj_io.h>

namespace py = pybind11;

namespace pyvkm {

  template<class MeshT>
  std::map<std::string, MeshT> wrap_read_composite_obj_str_map(std::string const& obj_path)
  {
    std::map<std::string, MeshT> scene;
    if(!vkm_obj_io::read_composite_obj_file<MeshT>(obj_path, scene))
      throw py::value_error("Read group OBJ failed");
    return scene;
  }

  template<class MeshT>
  std::map<size_t, std::map<size_t,MeshT> > wrap_read_composite_obj_int_map_map(std::string const& obj_path)
  {
    std::map<size_t, std::map<size_t,MeshT> > scene;
    if(!vkm_obj_io::read_composite_obj_file<MeshT>(obj_path, scene))
      throw py::value_error("Read composite OBJ failed");
    return scene;
  }

  template< typename MeshT >
  void wrap_write_composite_obj_str_map(
      std::string const& obj_path,
      std::map<std::string, MeshT> const& scene,
      std::string const& mat_file, std::string const& mat)
  {
    if(!vkm_obj_io::write_composite_obj_file<MeshT>(obj_path, scene, mat_file, mat))
      throw py::value_error("Write composite OBJ failed");
  }

  template< typename MeshT >
  void wrap_write_composite_obj_int_map_map(
      std::string const& obj_path,
      std::map<size_t, std::map<size_t, MeshT> > const& scene,
      std::string const& mat_file, std::string const& mat)
  {
    if(!vkm_obj_io::write_composite_obj_file<MeshT>(obj_path, scene, mat_file, mat))
      throw py::value_error("Write composite OBJ failed");
  }


  // main wrapper (called in header)
  void wrap_obj_io(py::module &m){

    // Convenience typedef
    typedef OpenMesh::TriMesh_ArrayKernelT<> TriMesh;
    typedef OpenMesh::PolyMesh_ArrayKernelT<> PolyMesh;

    // main class
    py::class_<vkm_obj_io> obj_io(m, "obj_io");

    // functions
    obj_io
      .def_static("read_group", &wrap_read_composite_obj_str_map<PolyMesh>)
      .def_static("read_group_tri", &wrap_read_composite_obj_str_map<TriMesh>)

      .def_static("read_composite", &wrap_read_composite_obj_int_map_map<PolyMesh>)
      .def_static("read_composite_tri", &wrap_read_composite_obj_int_map_map<TriMesh>)

      .def_static("write", &wrap_write_composite_obj_str_map<PolyMesh>,
          py::arg("filename"), py::arg("meshes"),
          py::arg("mat_file")="", py::arg("mat")="")
      .def_static("write", &wrap_write_composite_obj_str_map<TriMesh>,
          py::arg("filename"), py::arg("meshes"),
          py::arg("mat_file")="", py::arg("mat")="")
      .def_static("write", &wrap_write_composite_obj_int_map_map<PolyMesh>,
          py::arg("filename"), py::arg("meshes"),
          py::arg("mat_file")="", py::arg("mat")="")
      .def_static("write", &wrap_write_composite_obj_int_map_map<TriMesh>,
          py::arg("filename"), py::arg("meshes"),
          py::arg("mat_file")="", py::arg("mat")="")
      ;

  } // end wrap_obj_io

} // end namespace pyvkm

