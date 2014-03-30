#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/cycles.h"
#include "../../src/ring.h"
#include "common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

Helium::Size cyclomatic_number_1(const Molecule &mol, Helium::Size numComponents)
{
  return Helium::cyclomatic_number(mol, numComponents);
}

Helium::Size cyclomatic_number_2(const Molecule &mol)
{
  return Helium::cyclomatic_number(mol);
}

object cycle_membership(const Molecule &mol)
{
  std::vector<bool> cycleAtoms;
  std::vector<bool> cycleBonds;

  Helium::cycle_membership(mol, cycleAtoms, cycleBonds);

  return make_tuple(cycleAtoms, cycleBonds);
}

Helium::RingSet<Molecule> relevant_cycles_1(const Molecule &mol, Helium::Size cyclomaticNumber,
    const list &cyclicAtoms, const list &cyclicBonds)
{
  return Helium::relevant_cycles(mol, cyclomaticNumber,
      vector_from_list<bool>(cyclicAtoms), vector_from_list<bool>(cyclicBonds));
}

Helium::RingSet<Molecule> relevant_cycles_2(const Molecule &mol)
{
  return Helium::relevant_cycles(mol);
}

void export_rings()
{

  def("cyclomatic_number", &cyclomatic_number_1);
  def("cyclomatic_number", &cyclomatic_number_2);
  def("cycle_membership", &cycle_membership);
  def("relevant_cycles", &relevant_cycles_1);
  def("relevant_cycles", &relevant_cycles_2);

  class_<Helium::Ring<Molecule> >("Ring", no_init)
    .def("size", &Helium::Ring<Molecule>::size)
    .def("__len__", &Helium::Ring<Molecule>::size)
    .def("atoms", &Helium::Ring<Molecule>::atoms, return_value_policy<copy_const_reference>())
    .def("atom", &Helium::Ring<Molecule>::atom)
    .def("bonds", &Helium::Ring<Molecule>::bonds, return_value_policy<copy_const_reference>())
    .def("bond", &Helium::Ring<Molecule>::bond)
    .def("containsAtom", &Helium::Ring<Molecule>::containsAtom)
    .def("containsBond", &Helium::Ring<Molecule>::containsBond)
    ;

  class_<Helium::RingSet<Molecule> >("RingSet", init<const Molecule&>())
    .def("size", &Helium::RingSet<Molecule>::size)
    .def("__len__", &Helium::RingSet<Molecule>::size)
    .def("rings", &Helium::RingSet<Molecule>::rings, return_value_policy<copy_const_reference>())
    .def("ring", &Helium::RingSet<Molecule>::ring, return_internal_reference<>())
    .def("isAtomInRing", &Helium::RingSet<Molecule>::isAtomInRing)
    .def("isBondInRing", &Helium::RingSet<Molecule>::isBondInRing)
    .def("isAtomInRingSize", &Helium::RingSet<Molecule>::isAtomInRingSize)
    .def("isBondInRingSize", &Helium::RingSet<Molecule>::isBondInRingSize)
    .def("numRingBonds", &Helium::RingSet<Molecule>::numRingBonds)
    .def("numRingNbrs", &Helium::RingSet<Molecule>::numRingNbrs)
    .def("numRings", &Helium::RingSet<Molecule>::numRings)
    ;
}
