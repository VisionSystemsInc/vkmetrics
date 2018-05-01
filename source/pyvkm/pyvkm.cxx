#include <pybind11/pybind11.h>
#include "pyvkm_ground_truth.h"
#include "pyvkm_metrics.h"

namespace py = pybind11;

PYBIND11_MODULE(vkm, m)
{
  m.doc() =  "Python bindings for VK Metrics (VKM)";
  pyvkm::wrap_ground_truth(m);
  pyvkm::wrap_metrics(m);
}
