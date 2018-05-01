#ifndef pyvkm_ground_truth_h_
#define pyvkm_ground_truth_h_

#include <pybind11/pybind11.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace pyvkm {
  void wrap_ground_truth(pybind11::module &m);
}

#endif
