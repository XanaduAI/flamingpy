#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <lemon/connectivity.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>

typedef lemon::ListGraph Graph;
typedef Graph::EdgeMap<float> Weights;
typedef lemon::MaxWeightedPerfectMatching<Graph,Weights> MWPM;

std::vector<int> mwpm(const std::vector<float> & adjacency, int num_nodes)
{
  Graph g;
  Weights w(g);

  for (int n=0; n < num_nodes; n++){
    g.addNode();
  }
 
  for (int i = 0; i < num_nodes; i++){
    for (int j = i+1; j < num_nodes; j++){
      // if (adjacency[i*num_nodes + j] != 0.0){
        Graph::Edge e = g.addEdge(g.nodeFromId(i), g.nodeFromId(j));
        w[e] = adjacency[i*num_nodes + j];
      // }
    }
  }
  
  MWPM pm(g, w);  
  bool success = pm.run();
  std::vector<int> matching_vec;

  for (Graph::EdgeIt e(g); e != lemon::INVALID; ++e) {    
    if (pm.matching(e)){
      matching_vec.push_back(g.id(g.u(e)));
      matching_vec.push_back(g.id(g.v(e)));
    }
  }
  
  return matching_vec;
}

// ----------------
// Python interface
// ----------------

namespace py = pybind11;

// wrap C++ function with NumPy array IO
py::array py_mwpm(py::array_t<float, py::array::c_style | py::array::forcecast> array)
{
  // check input dimensions
  if ( array.ndim()     != 2 )
    throw std::runtime_error("Input should be 2-D NumPy array");

  // allocate std::vector (to pass to the C++ function)
  std::vector<float> pos(array.size());

  // copy py::array -> std::vector
  std::memcpy(pos.data(),array.data(),array.size()*sizeof(float));

  // call pure C++ function
  std::vector<int> result = mwpm(pos, array.shape()[0]);

  Py_ssize_t              ndim    = 2;
  std::vector<Py_ssize_t> shape   = { static_cast<Py_ssize_t>(result.size() / 2) , 2 };
  std::vector<Py_ssize_t> strides = { sizeof(int)*2 , sizeof(int) };

  // return 2-D NumPy array
  return py::array(py::buffer_info(
    result.data(),                           /* data as contiguous array  */
    sizeof(int),                          /* size of one scalar        */
    py::format_descriptor<int>::format(), /* data type                 */
    ndim,                                    /* number of dimensions      */
    shape,                                   /* shape of the matrix       */
    strides                                  /* strides for each axis     */
  ));
}

PYBIND11_MODULE(lemonpy,m)
{
  m.doc() = "pybind11 lemon plugin";
  m.def("mwpm", &py_mwpm, "GetMaxWeightedPerfectMatching");
}
