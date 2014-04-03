#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/components.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

void export_components()
{

  def("connected_bond_components", &Helium::connected_bond_components<Molecule>);
  def("connected_atom_components", &Helium::connected_atom_components<Molecule>);
  def("num_connected_components", &Helium::num_connected_components<Molecule>);

}
