#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/smiles.h"
#include "common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

Molecule* create_molecule()
{
  return new Molecule();
}

bool atom_equal(const Molecule::atom_type &self, const Molecule::atom_type &other)
{
  return self.index() == other.index();
}

bool bond_equal(const Molecule::bond_type &self, const Molecule::bond_type &other)
{
  return self.index() == other.index();
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

  class_<Molecule, boost::noncopyable>("Molecule")
    .def("__init__", make_constructor(&create_molecule))
    .def("numAtoms", &Molecule::numAtoms)
    .def("numBonds", &Molecule::numBonds)
    .def("atoms", &Molecule::atoms)
    .def("bonds", &Molecule::bonds)
    .def("atom", &Molecule::atom)
    .def("bond", &Molecule::bond)
    .def("addAtom", &Molecule::addAtom)
    .def("removeAtom", &Molecule::removeAtom)
    .def("addBond", &Molecule::addBond)
    .def("removeBond", &Molecule::removeBond)
    .def("clear", &Molecule::clear)
    ;

}
