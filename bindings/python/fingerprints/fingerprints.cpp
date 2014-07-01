#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/fingerprints/fingerprints.h"
#include "../common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

Fingerprint* path_fingerprint(const Molecule &mol, int size = 7, int numWords = 16, int hashPrime = 1021)
{
  Helium::Word *fp = new Helium::Word[numWords];
  Helium::path_fingerprint(mol, fp, size, numWords, hashPrime);
  return new Fingerprint(fp, numWords, true);
}

Fingerprint* tree_fingerprint(const Molecule &mol, int size = 7, int numWords = 16, int hashPrime = 1021)
{
  Helium::Word *fp = new Helium::Word[numWords];
  Helium::tree_fingerprint(mol, fp, size, numWords, hashPrime);
  return new Fingerprint(fp, numWords, true);
}

Fingerprint* subgraph_fingerprint(const Molecule &mol, int size = 7, int numWords = 16, int hashPrime = 1021)
{
  Helium::Word *fp = new Helium::Word[numWords];
  Helium::subgraph_fingerprint(mol, fp, size, numWords, hashPrime);
  return new Fingerprint(fp, numWords, true);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(path_fingerprint_overloads, path_fingerprint, 1, 4)
BOOST_PYTHON_FUNCTION_OVERLOADS(tree_fingerprint_overloads, tree_fingerprint, 1, 4)
BOOST_PYTHON_FUNCTION_OVERLOADS(subgraph_fingerprint_overloads, subgraph_fingerprint, 1, 4)

void export_fingerprints()
{

  def("path_fingerprint", path_fingerprint, path_fingerprint_overloads()[return_value_policy<manage_new_object>()]);
  def("tree_fingerprint", tree_fingerprint, tree_fingerprint_overloads()[return_value_policy<manage_new_object>()]);
  def("subgraph_fingerprint", subgraph_fingerprint, subgraph_fingerprint_overloads()[return_value_policy<manage_new_object>()]);

}
