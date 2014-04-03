#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/canonical.h"
#include "../common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;
  
template<typename AtomInvariantType, typename BondInvariantType>
PyObject* canonicalize_component(const Molecule &mol, const list &symmetry,
      const AtomInvariantType &atomInvariant, const BondInvariantType &bondInvariant)
{
  std::pair<std::vector<Helium::Index>, std::vector<unsigned long> > result;
  result = Helium::canonicalize_component(mol, vector_from_list<unsigned long>(symmetry),
      atomInvariant, bondInvariant);
  tuple *t = new tuple(make_tuple(result.first, result.second));
  return t->ptr();
}

template<typename AtomInvariantType, typename BondInvariantType>
PyObject* canonicalize(const Molecule &mol, const list &symmetry,
    const AtomInvariantType &atomInvariant, const BondInvariantType &bondInvariant,
    const list &atomComponents, const list &bondComponents)
{
  std::pair<std::vector<Helium::Index>, std::vector<unsigned long> > result;
  result = Helium::canonicalize(mol, vector_from_list<unsigned long>(symmetry),
      atomInvariant, bondInvariant, vector_from_list<unsigned int>(atomComponents),
      vector_from_list<unsigned int>(bondComponents));
  tuple *t = new tuple(make_tuple(result.first, result.second));
  return t->ptr();
}

void export_canonical()
{

  def("canonicalize_component", &canonicalize_component<Helium::DefaultAtomInvariant,
      Helium::DefaultBondInvariant>);
  def("canonicalize_component", &canonicalize_component<AtomInvariant,
      Helium::DefaultBondInvariant>);
  def("canonicalize_component", &canonicalize_component<Helium::DefaultAtomInvariant,
      BondInvariant>);
  def("canonicalize_component", &canonicalize_component<AtomInvariant, BondInvariant>);

  def("canonicalize", &canonicalize<Helium::DefaultAtomInvariant, Helium::DefaultBondInvariant>);
  def("canonicalize", &canonicalize<AtomInvariant, Helium::DefaultBondInvariant>);
  def("canonicalize", &canonicalize<Helium::DefaultAtomInvariant, BondInvariant>);
  def("canonicalize", &canonicalize<AtomInvariant, BondInvariant>);

}
