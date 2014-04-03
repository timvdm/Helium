#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/extendedconnectivities.h"
#include "../common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

template<typename AtomInvariantType>
std::vector<unsigned long> extended_connectivities(const Molecule &mol,
    const AtomInvariantType &atomInvariant)
{
  return Helium::extended_connectivities(mol, atomInvariant);
}

void export_extended_connectivities()
{

  def("extended_connectivities", &extended_connectivities<AtomInvariant>);
  def("extended_connectivities", &extended_connectivities<Helium::DefaultAtomInvariant>);

}
