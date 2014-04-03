#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/floydwarshall.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

Helium::DistanceMatrix* floyd_warshall(const Molecule &mol)
{
  return new Helium::DistanceMatrix(Helium::floyd_warshall(mol));
}

void export_floyd_warshall()
{

  def("floyd_warshall", &floyd_warshall, return_value_policy<manage_new_object>());

}
