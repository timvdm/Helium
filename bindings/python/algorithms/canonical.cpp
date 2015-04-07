#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/canonical.h"
#include "../common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

template<typename AtomInvariantType, typename BondInvariantType>
object canonicalize_component(const Molecule &mol, const list &symmetry,
      const AtomInvariantType &atomInvariant, const BondInvariantType &bondInvariant)
{
  if (len(symmetry) != mol.numAtoms())
    throw std::runtime_error("Invalid symmetry parameter, size should be equal to the number of atoms");

  std::pair<std::vector<Helium::Index>, std::vector<unsigned long> > result;
  result = Helium::canonicalize_component(mol, vector_from_list<unsigned long>(symmetry),
      atomInvariant, bondInvariant);

  return boost::python::make_tuple(result.first, result.second);
}

template<typename AtomInvariantType, typename BondInvariantType>
object canonicalize(const Molecule &mol, const list &symmetry,
    const AtomInvariantType &atomInvariant, const BondInvariantType &bondInvariant,
    const list &atomComponents, const list &bondComponents)
{
  if (len(symmetry) != mol.numAtoms())
    throw std::runtime_error("Invalid symmetry parameter, size should be equal to the number of atoms");
  if (len(atomComponents) != mol.numAtoms())
    throw std::runtime_error("Invalid atomComponents parameter, size should be equal to the number of atoms");
  if (len(bondComponents) != mol.numBonds())
    throw std::runtime_error("Invalid bondComponents parameter, size should be equal to the number of bonds");

  std::pair<std::vector<Helium::Index>, std::vector<unsigned long> > result;
  result = Helium::canonicalize(mol, vector_from_list<unsigned long>(symmetry),
      atomInvariant, bondInvariant, vector_from_list<unsigned int>(atomComponents),
      vector_from_list<unsigned int>(bondComponents));

  return boost::python::make_tuple(result.first, result.second);
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
