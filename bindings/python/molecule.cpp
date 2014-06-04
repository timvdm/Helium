#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/smiles.h"
#include "common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

bool atom_equal(const Molecule::atom_type &self, const Molecule::atom_type &other)
{
  return self.index() == other.index();
}

bool bond_equal(const Molecule::bond_type &self, const Molecule::bond_type &other)
{
  return self.index() == other.index();
}

Molecule* Molecule_ctor_1(const Molecule &source, const list &atoms, const list &bonds, bool adjustHydrogens)
{
  return new Molecule(source, vector_from_list<bool>(atoms), vector_from_list<bool>(bonds), adjustHydrogens);
}

Molecule* Molecule_ctor_2(const Molecule &source, const list &atoms, const list &bonds)
{
  return Molecule_ctor_1(source, atoms, bonds, true);
}

Molecule::bond_type bond_1(const Molecule &mol, Helium::Index index)
{
  return mol.bond(index);
}

Molecule::bond_type bond_2(const Molecule &mol, const Molecule::atom_type &source, const Molecule::atom_type &target)
{
  return mol.bond(source, target);
}

void reset_implicit_hydrogens_1(Molecule &mol)
{
  reset_implicit_hydrogens(mol);
}

void reset_implicit_hydrogens_2(Molecule &mol, const Molecule::atom_type &atom)
{
  reset_implicit_hydrogens(mol, atom);
}

void export_molecule()
{
  IteratorWrapper<Molecule::atom_type, Molecule::atom_iter>().wrap("AtomIterator");
  IteratorWrapper<Molecule::bond_type, Molecule::bond_iter>().wrap("BondIterator");
  IteratorWrapper<Molecule::atom_type, Molecule::nbr_iter>().wrap("NbrIterator");
  IteratorWrapper<Molecule::bond_type, Molecule::incident_iter>().wrap("IncidentIterator");

  class_<Molecule::atom_type>("Atom", no_init)
    .def("index", &Molecule::atom_type::index)
    .def("degree", &Molecule::atom_type::degree)
    .def("bonds", &Molecule::atom_type::bonds)
    .def("nbrs", &Molecule::atom_type::nbrs)
    .def("isAromatic", &Molecule::atom_type::isAromatic)
    .def("setAromatic", &Molecule::atom_type::setAromatic)
    .def("element", &Molecule::atom_type::element)
    .def("setElement", &Molecule::atom_type::setElement)
    .def("mass", &Molecule::atom_type::mass)
    .def("setMass", &Molecule::atom_type::setMass)
    .def("hydrogens", &Molecule::atom_type::hydrogens)
    .def("setHydrogens", &Molecule::atom_type::setHydrogens)
    .def("charge", &Molecule::atom_type::charge)
    .def("setCharge", &Molecule::atom_type::setCharge)
    .def("isHydrogen", &Molecule::atom_type::isHydrogen)
    .def("isCarbon", &Molecule::atom_type::isCarbon)
    .def("isNitrogen", &Molecule::atom_type::isNitrogen)
    .def("isOxygen", &Molecule::atom_type::isOxygen)
    .def("isPhosphorus", &Molecule::atom_type::isPhosphorus)
    .def("isSulfur", &Molecule::atom_type::isSulfur)
    .def("heavyDegree", &Molecule::atom_type::heavyDegree)
    .def("boSum", &Molecule::atom_type::boSum)
    .def("valence", &Molecule::atom_type::valence)
    .def("connectivity", &Molecule::atom_type::connectivity)
    .def("__eq__", &atom_equal)
    ;

  class_<Molecule::bond_type>("Bond", no_init)
    .def("index", &Molecule::bond_type::index)
    .def("source", &Molecule::bond_type::source)
    .def("target", &Molecule::bond_type::target)
    .def("other", &Molecule::bond_type::other)
    .def("isAromatic", &Molecule::bond_type::isAromatic)
    .def("setAromatic", &Molecule::bond_type::setAromatic)
    .def("order", &Molecule::bond_type::order)
    .def("setOrder", &Molecule::bond_type::setOrder)
    .def("__eq__", &bond_equal)
    ;

  class_<Molecule, boost::shared_ptr<Molecule>, boost::noncopyable>("Molecule")
    .def(init<const Molecule&>())
    .def("__init__", make_constructor(&Molecule_ctor_1))
    .def("__init__", make_constructor(&Molecule_ctor_2))
    .def("numAtoms", &Molecule::numAtoms)
    .def("numBonds", &Molecule::numBonds)
    .def("atoms", &Molecule::atoms)
    .def("bonds", &Molecule::bonds)
    .def("atom", &Molecule::atom)
    .def("bond", &bond_1)
    .def("bond", &bond_2)
    .def("addAtom", &Molecule::addAtom)
    .def("removeAtom", &Molecule::removeAtom)
    .def("addBond", &Molecule::addBond)
    .def("removeBond", &Molecule::removeBond)
    .def("clear", &Molecule::clear)
    ;

  def("make_hydrogens_explicit", &Helium::make_hydrogens_explicit<Molecule>);
  def("make_hydrogens_implicit", &Helium::make_hydrogens_implicit<Molecule>);
  def("reset_implicit_hydrogens", &reset_implicit_hydrogens_1);
  def("reset_implicit_hydrogens", &reset_implicit_hydrogens_2);

}
