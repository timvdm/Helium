#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/enumeratesubgraphs.h"
#include "../common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

DATA_MEMBER_TO_FUNCTION(Helium::Subgraph, std::vector<bool>, atoms);
DATA_MEMBER_TO_FUNCTION(Helium::Subgraph, std::vector<bool>, bonds);

struct SubgraphCallback
{
  virtual void call(const Helium::Subgraph &subgraph) = 0;

  void operator()(const Helium::Subgraph &subgraph)
  {
    call(subgraph);
  }
};

struct SubgraphCallbackWrapper : SubgraphCallback, wrapper<SubgraphCallback>
{
  void call(const Helium::Subgraph &subgraph)
  {
    this->get_override("__call__")(boost::ref(subgraph));
  }
};

void enumerate_subgraphs(const Molecule &mol, SubgraphCallback &callback,
    int maxSize, bool trees = false)
{
  Helium::enumerate_subgraphs(mol, callback, maxSize, trees);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(enumerate_subgraphs_overload, enumerate_subgraphs, 3, 4);

void export_enumerate_subgraphs()
{

  class_<SubgraphCallbackWrapper, boost::noncopyable>("SubgraphCallback")
    .def("__call__", &SubgraphCallback::call);

  class_<Helium::Subgraph>("Subgraph", no_init)
    .def("hashable", &Helium::Subgraph::hashable)
    .def("atoms", make_function(&atoms, return_value_policy<copy_const_reference>()))
    .def("bonds", make_function(&bonds, return_value_policy<copy_const_reference>()))
    ;

  def("enumerate_subgraphs", enumerate_subgraphs, enumerate_subgraphs_overload());

}
