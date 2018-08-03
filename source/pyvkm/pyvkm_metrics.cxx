#include "pyvkm_metrics.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <vkm_obj_io.h>
#include <vkm_metrics.h>

namespace py = pybind11;

namespace pyvkm {

  // wrap bool functions to instead throw errors
  void load_ground_truth_model(vkm_metrics &_self, std::string const& path) {
    if (!_self.load_ground_truth_model(path))
      throw pybind11::value_error("load_ground_truth_model failure!");
  }

  void load_ground_planes(vkm_metrics &_self, std::string const& path) {
    if (!_self.load_ground_planes(path))
      throw pybind11::value_error("load_ground_planes failure!");
  }

  void load_simply_connected_test_model(vkm_metrics &_self, std::string const& path, std::string site_name = "") {
    if (!_self.load_simply_connected_test_model(path,site_name))
      throw pybind11::value_error("load_simply_connected_test_model failure!");
  }

  void save_transformed_regions_as_meshes(vkm_metrics &_self, std::string const& path, std::string site_name) {
    if (!_self.save_transformed_regions_as_meshes(path,site_name))
      throw pybind11::value_error("save_transformed_regions_as_meshes failure!");
  }


  // main wrapper (called in header)
  void wrap_metrics(py::module &m){

    // main metrics class
    py::class_<vkm_metrics> metrics(m, "metrics");

    // functions
    metrics
      .def(py::init<>())
      .def("set_verbose", &vkm_metrics::set_verbose)
      .def("set_translation", &vkm_metrics::set_translation)
      .def("load_ground_truth_model", &load_ground_truth_model)
      .def("load_ground_planes", &load_ground_planes)
      .def("load_simply_connected_test_model", &load_simply_connected_test_model)
      .def("delete_isolated_vertices", &vkm_metrics::delete_isolated_vertices)
      .def("construct_xy_regions", &vkm_metrics::construct_xy_regions)
      .def("match_xy_regions", &vkm_metrics::match_xy_regions)
      .def("translate_test_model_xy", &vkm_metrics::translate_test_model_xy)
      .def("find_z_offset", &vkm_metrics::find_z_offset)
      .def("compute_best_match_2d_score_stats", &vkm_metrics::compute_best_match_2d_score_stats)
      .def("save_transformed_regions_as_meshes", &save_transformed_regions_as_meshes)
      .def("compute_best_match_3d_score_stats", &vkm_metrics::compute_best_match_3d_score_stats)
      .def("json", &vkm_metrics::json)
      .def("print_score_array", &vkm_metrics::print_score_array)
      ;

  }

} // end namespace pyvkm

