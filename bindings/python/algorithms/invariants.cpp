#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/invariants.h"
#include "../common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

struct AtomInvariantWrapper : AtomInvariant, wrapper<AtomInvariant>
{
  unsigned long call_1(const Molecule &mol, const Molecule::atom_type &atom) const
  {
    return this->get_override("__call__")(boost::ref(mol), boost::ref(atom));
  }
};

struct BondInvariantWrapper : BondInvariant, wrapper<BondInvariant>
{
  unsigned long call_1(const Molecule &mol, const Molecule::bond_type &bond) const
  {
    return this->get_override("__call__")(boost::ref(mol), boost::ref(bond));
  }
};

unsigned long DefaultAtomInvariant_call(const Helium::DefaultAtomInvariant &invariant, const Molecule &mol,
    const Molecule::atom_type &atom)
{
  return invariant(mol, atom);
}

unsigned long DefaultBondInvariant_call(const Helium::DefaultBondInvariant &invariant, const Molecule &mol,
    const Molecule::bond_type &bond)
{
  return invariant(mol, bond);
}

void export_default_atom_invariant()
{
  scope in_DefaultAtomInvariant = class_<Helium::DefaultAtomInvariant, boost::noncopyable>("DefaultAtomInvariant", init<>())
    .def(init<int>())
    .def("__call__", &DefaultAtomInvariant_call)
    ;

  enum_<Helium::DefaultAtomInvariant::Invariants>("Invariants")
    .value("None", Helium::DefaultAtomInvariant::None)
    .value("Element", Helium::DefaultAtomInvariant::Element)
    .value("Mass", Helium::DefaultAtomInvariant::Mass)
    .value("Charge", Helium::DefaultAtomInvariant::Charge)
    .value("Degree", Helium::DefaultAtomInvariant::Degree)
    .value("Aromatic", Helium::DefaultAtomInvariant::Aromatic)
    .value("All", Helium::DefaultAtomInvariant::All)
    ;
}

void export_default_bond_invariant()
{
  scope in_DefaultBondInvariant = class_<Helium::DefaultBondInvariant, boost::noncopyable>("DefaultBondInvariant", init<>())
    .def(init<int>())
    .def("__call__", &DefaultBondInvariant_call)
    ;

  enum_<Helium::DefaultBondInvariant::Invariants>("Invariants")
    .value("None", Helium::DefaultBondInvariant::None)
    .value("Order", Helium::DefaultBondInvariant::Order)
    .value("Aromatic", Helium::DefaultBondInvariant::Aromatic)
    .value("All", Helium::DefaultBondInvariant::All)
    ;
}

void export_invariants()
{

  class_<AtomInvariantWrapper, boost::noncopyable>("AtomInvariant")
    .def("__call__", &AtomInvariant::call_1);
  class_<BondInvariantWrapper, boost::noncopyable>("BondInvariant")
    .def("__call__", &BondInvariant::call_1);

  export_default_atom_invariant();
  export_default_bond_invariant();

}
