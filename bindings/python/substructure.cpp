#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/substructure.h"
#include "common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

Helium::Substructure<Molecule>* create_substructure(const Molecule &mol, const list &atoms, const list &bonds)
{
  return new Helium::Substructure<Molecule>(mol, vector_from_list<bool>(atoms),
      vector_from_list<bool>(bonds));
}

bool atom_equal(const Helium::Substructure<Molecule>::atom_type &self, const Helium::Substructure<Molecule>::atom_type &other);
bool bond_equal(const Helium::Substructure<Molecule>::bond_type &self, const Helium::Substructure<Molecule>::bond_type &other);

void export_substructure()
{
  IteratorWrapper<Helium::Substructure<Molecule>::atom_type, Helium::Substructure<Molecule>::atom_iter>().wrap("AtomIterator");
  IteratorWrapper<Helium::Substructure<Molecule>::bond_type, Helium::Substructure<Molecule>::bond_iter>().wrap("BondIterator");
  IteratorWrapper<Helium::Substructure<Molecule>::atom_type, Helium::Substructure<Molecule>::nbr_iter>().wrap("NbrIterator");
  IteratorWrapper<Helium::Substructure<Molecule>::bond_type, Helium::Substructure<Molecule>::incident_iter>().wrap("IncidentIterator");

  /*
  class_<Helium::Substructure<Molecule>::atom_type>("Atom", no_init)
    .def("index", &Helium::Substructure<Molecule>::atom_type::index)
    .def("degree", &Helium::Substructure<Molecule>::atom_type::degree)
    .def("bonds", &Helium::Substructure<Molecule>::atom_type::bonds)
    .def("nbrs", &Helium::Substructure<Molecule>::atom_type::nbrs)
    .def("isAromatic", &Helium::Substructure<Molecule>::atom_type::isAromatic)
    .def("element", &Helium::Substructure<Molecule>::atom_type::element)
    .def("mass", &Helium::Substructure<Molecule>::atom_type::mass)
    .def("hydrogens", &Helium::Substructure<Molecule>::atom_type::hydrogens)
    .def("charge", &Helium::Substructure<Molecule>::atom_type::charge)
    .def("__eq__", &atom_equal)
    ;

  class_<Helium::Substructure<Molecule>::bond_type>("Bond", no_init)
    .def("index", &Helium::Substructure<Molecule>::bond_type::index)
    .def("source", &Helium::Substructure<Molecule>::bond_type::source)
    .def("target", &Helium::Substructure<Molecule>::bond_type::target)
    .def("other", &Helium::Substructure<Molecule>::bond_type::other)
    .def("isAromatic", &Helium::Substructure<Molecule>::bond_type::isAromatic)
    .def("order", &Helium::Substructure<Molecule>::bond_type::order)
    .def("__eq__", &bond_equal)
    ;
    */

  class_<Helium::Substructure<Molecule>, boost::noncopyable>("Substructure", no_init)
    .def("__init__", make_constructor(&create_substructure))
    .def("numAtoms", &Helium::Substructure<Molecule>::numAtoms)
    .def("numBonds", &Helium::Substructure<Molecule>::numBonds)
    .def("atoms", &Helium::Substructure<Molecule>::atoms)
    .def("bonds", &Helium::Substructure<Molecule>::bonds)
    .def("atom", &Helium::Substructure<Molecule>::atom)
    .def("bond", &Helium::Substructure<Molecule>::bond)
    ;

}
