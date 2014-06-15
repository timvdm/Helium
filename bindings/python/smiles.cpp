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

std::string write_3(Helium::Smiles &self, const Molecule &mol, const dict &atomClasses)
{
  return self.write(mol, map_from_dict<Helium::Index, int>(atomClasses));
}

std::string write_4(Helium::Smiles &self, const Molecule &mol, const dict &atomClasses, int flags)
{
  return self.write(mol, map_from_dict<Helium::Index, int>(atomClasses), flags);
}

std::string write_5(Helium::Smiles &self, const Molecule &mol, const list &order)
{
  if (len(order) != mol.numAtoms())
    throw std::runtime_error("Invalid order parameter, size not equal to the number of atoms");
  return self.write(mol, vector_from_list<Helium::Index>(order));
}

std::string write_6(Helium::Smiles &self, const Molecule &mol, const list &order, int flags)
{
  if (len(order) != mol.numAtoms())
    throw std::runtime_error("Invalid order parameter, size not equal to the number of atoms");
  return self.write(mol, vector_from_list<Helium::Index>(order), flags);
}

std::string write_7(Helium::Smiles &self, const Molecule &mol, const list &order,
    const dict &atomClasses)
{
  if (len(order) != mol.numAtoms())
    throw std::runtime_error("Invalid order parameter, size not equal to the number of atoms");
  return self.write(mol, vector_from_list<Helium::Index>(order), map_from_dict<Helium::Index, int>(atomClasses));
}

std::string write_8(Helium::Smiles &self, const Molecule &mol, const list &order,
    const dict &atomClasses, int flags)
{
  if (len(order) != mol.numAtoms())
    throw std::runtime_error("Invalid order parameter, size not equal to the number of atoms");
  return self.write(mol, vector_from_list<Helium::Index>(order), map_from_dict<Helium::Index, int>(atomClasses), flags);
}

std::string write_canonical_1(Helium::Smiles &self, const Molecule &mol)
{
  return self.writeCanonical(mol);
}

std::string write_canonical_2(Helium::Smiles &self, const Molecule &mol, int flags)
{
  return self.writeCanonical(mol, flags);
}

std::string write_canonical_3(Helium::Smiles &self, const Molecule &mol, const dict &atomClasses)
{
  return self.writeCanonical(mol, map_from_dict<Helium::Index, int>(atomClasses));
}

std::string write_canonical_4(Helium::Smiles &self, const Molecule &mol, const dict &atomClasses, int flags)
{
  return self.writeCanonical(mol, map_from_dict<Helium::Index, int>(atomClasses), flags);
}

void export_smiles()
{

  scope in_Smiles = class_<Helium::Smiles>("Smiles")
    .def("read", &Helium::Smiles::read<Molecule>)
    .def("write", &write_1)
    .def("write", &write_2)
    .def("write", &write_3)
    .def("write", &write_4)
    .def("write", &write_5)
    .def("write", &write_6)
    .def("write", &write_7)
    .def("write", &write_8)
    .def("writeCanonical", &write_canonical_1)
    .def("writeCanonical", &write_canonical_2)
    .def("writeCanonical", &write_canonical_3)
    .def("writeCanonical", &write_canonical_4)
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
