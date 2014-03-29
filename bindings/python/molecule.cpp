#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"

using namespace boost::python;

class PyMolecule : public Helium::Chemist::Molecule
{
  public:
    PyMolecule() : Molecule()
    {
    }

    atom_iter beginAtoms() const
    {
      return atoms().begin();
    }

    atom_iter endAtoms() const
    {
      return atoms().end();
    }
};

void export_molecule()
{
  class_<Helium::Chemist::Molecule::atom_type>("Atom", "Helium atom datastructure")
    .def("element", &Helium::Chemist::Molecule::atom_type::element, "Get the atom's element")

    ;

  class_<PyMolecule, boost::noncopyable>("Molecule", "Helium molecule datastructure")

    .def("numAtoms", &PyMolecule::numAtoms, "Get the number of atoms in the molecule")
    .def("numBonds", &PyMolecule::numBonds, "Get the number of bonds in the molecule")

    .def("atoms", range(&PyMolecule::beginAtoms, &PyMolecule::endAtoms))


    .def("addAtom", &PyMolecule::addAtom, "Add an atom to the molecule") 



    ;
}
