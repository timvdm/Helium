#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/smiles.h"
#include "common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

std::string write_1(Helium::Smiles &self, const Molecule &mol)
{
  return self.write(mol);
}

std::string write_2(Helium::Smiles &self, const Molecule &mol, int flags)
{
  return self.write(mol, flags);
}

std::string write_3(Helium::Smiles &self, const Molecule &mol)
{
  return self.writeCanonical(mol);
}

std::string write_4(Helium::Smiles &self, const Molecule &mol, const list &order)
{
  return self.write(mol, vector_from_list<Helium::Index>(order));
}

std::string write_5(Helium::Smiles &self, const Molecule &mol, const list &order, int flags)
{
  return self.write(mol, vector_from_list<Helium::Index>(order), flags);
}

void export_smiles()
{

  scope in_Smiles = class_<Helium::Smiles>("Smiles")
    .def("read", &Helium::Smiles::read<Molecule>)
    .def("write", &write_1)
    .def("write", &write_2)
    .def("writeCanonical", &write_3)
    .def("write", &write_4)
    .def("write", &write_5)
    .def("error", &Helium::Smiles::error, return_internal_reference<>())
    ;

  enum_<Helium::Smiles::Flags>("Flags")
    .value("None", Helium::Smiles::None)
    .value("Mass", Helium::Smiles::Mass)
    .value("Charge", Helium::Smiles::Charge)
    .value("Hydrogens", Helium::Smiles::Hydrogens)
    .value("Order", Helium::Smiles::Order)
    .value("All", Helium::Smiles::All)
    ;

}
