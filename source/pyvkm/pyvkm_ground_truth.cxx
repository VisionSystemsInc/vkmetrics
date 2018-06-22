#include "pyvkm_ground_truth.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <vkm_obj_io.h>
#include <vkm_ground_truth.h>

namespace py = pybind11;

namespace pyvkm {

  // wrap bool functions to instead throw errors
  void load_ground_truth_img_regions(vkm_ground_truth &_self, std::string const& path) {
    if (!_self.load_ground_truth_img_regions(path))
      throw pybind11::value_error("load_ground_truth_img_regions failure!");
  }

  void load_surface_types(vkm_ground_truth &_self, std::string const& path) {
    if (!_self.load_surface_types(path))
      throw pybind11::value_error("load_surface_types failure!");
  }

  void load_dem_image(vkm_ground_truth &_self, std::string const& path) {
    if (!_self.load_dem_image(path))
      throw pybind11::value_error("load_dem_image failure!");
  }

  void load_site_perimeter(vkm_ground_truth &_self, std::string const& name, std::string const& path) {
    if (!_self.load_site_perimeter(name,path))
      throw pybind11::value_error("load_site_perimeter failure!");
  }

  void write_xy_polys(vkm_ground_truth &_self, std::string const& name, std::string const& path) {
    if (!_self.write_xy_polys(name,path))
      throw pybind11::value_error("write_xy_polys failure!");
  }

  void write_processed_ground_truth(vkm_ground_truth &_self, std::string const& path) {
    if (!_self.write_processed_ground_truth(path))
      throw pybind11::value_error("write_processed_ground_truth failure!");
  }

  void compute_img_to_xy_trans(vkm_ground_truth &_self, std::string const& path) {
    if (!_self.compute_img_to_xy_trans(path))
      throw pybind11::value_error("compute_img_to_xy_trans failure!");
  }

  void ensure_consistent_topology(vkm_ground_truth &_self, size_t outer_index, mc_region_2d& mcr) {
    if (!_self.ensure_consistent_topology(outer_index,mcr))
      throw pybind11::value_error("ensure_consistent_topology failure!");
  }

  void construct_extruded_gt_model(vkm_ground_truth &_self, std::string const& name) {
    if (!_self.construct_extruded_gt_model(name))
      throw pybind11::value_error("construct_extruded_gt_model failure!");
  }

  void write_extruded_gt_model(vkm_ground_truth &_self, std::string const& name, std::string const& path) {
    if (!_self.write_extruded_gt_model(name,path))
      throw pybind11::value_error("write_extruded_gt_model failure!");
  }

  void write_ground_planes(vkm_ground_truth &_self, std::string const& path) {
    if (!_self.write_ground_planes(path))
      throw pybind11::value_error("write_ground_planes failure!");
  }

  // main wrapper (called in header)
  void wrap_ground_truth(py::module &m){

    // main ground truth class
    py::class_<vkm_ground_truth> ground_truth(m, "ground_truth");

    // enumeration(s) attached to class
    py::enum_<vkm_ground_truth::surface_t>(ground_truth, "surface_type")
      .value("FLAT",    vkm_ground_truth::FLAT)
      .value("SLOPED",  vkm_ground_truth::SLOPED)
      .value("ARCHED",  vkm_ground_truth::ARCHED)
      .value("DOMED",   vkm_ground_truth::DOMED)
      .value("MISC",    vkm_ground_truth::MISC)
      ;

    // functions
    ground_truth
      .def(py::init<double>(), py::arg("z_off")=0.0)
      .def("set_verbose", &vkm_ground_truth::set_verbose)
      .def("load_ground_truth_img_regions", &load_ground_truth_img_regions)
      .def("load_surface_types", &load_surface_types)
      .def("load_dem_image", &load_dem_image)
      .def("load_site_perimeter", &load_site_perimeter)
      .def("write_xy_polys", &write_xy_polys)
      .def("write_processed_ground_truth", &write_processed_ground_truth)
      // .def_static("load_processed_ground_truth", &vkm_ground_truth::load_processed_ground_truth)
      .def("img_to_xy_trans", &compute_img_to_xy_trans)
      .def("img_to_xy_trans", &vkm_ground_truth::set_img_to_xy_trans)
      .def("snap_image_region_vertices", &vkm_ground_truth::snap_image_region_vertices,
          py::arg("tol")=2.5)
      .def("process_region_containment", &vkm_ground_truth::process_region_containment)
      .def("fit_region_planes", &vkm_ground_truth::fit_region_planes)
      .def("ensure_consistent_topology", &ensure_consistent_topology)
      .def("construct_polygon_soup", &vkm_ground_truth::construct_polygon_soup)
      .def("convert_to_meshes", &vkm_ground_truth::convert_to_meshes)
      // .def("region_meshes", &vkm_ground_truth::region_meshes)
      .def("construct_extruded_gt_model", &construct_extruded_gt_model)
      .def("write_extruded_gt_model", &write_extruded_gt_model)
      // .def_static("extrude_base_pts", &vkm_ground_truth::extrude_base_pts)
      .def("convert_img_regions_to_meshes", &vkm_ground_truth::convert_img_regions_to_meshes)
      //.def("img_meshes", &vkm_ground_truth::img_meshes)
      .def("write_ground_planes", &write_ground_planes)

      ;

  } // end wrap_ground_truth

} // end namespace pyvkm
