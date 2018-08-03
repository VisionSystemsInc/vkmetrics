#include <pybind11/pybind11.h>
#include "pyvkm_ground_truth.h"
#include "pyvkm_metrics.h"
#include "pyvkm_obj_io.h"

namespace py = pybind11;

// helper function to check if py::module import exists
bool import_exists(std::string const& library_name)
{
  py::module importlib = py::module::import("importlib");
  return (!importlib.attr("find_loader")(library_name.c_str()).is_none());
}

PYBIND11_MODULE(vkm, m)
{
  m.doc() = "Python bindings for VK Metrics (VKM)";
  pyvkm::wrap_ground_truth(m);
  pyvkm::wrap_metrics(m);
  pyvkm::wrap_obj_io(m);

  // include openmesh submodule (wrap_obj_io returns openmesh objects)
  if (import_exists("openmesh"))
    py::module::import("openmesh");
}

