#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/fileio/moleculefile.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

bool readMolecule_1(Helium::MoleculeFile &file, Molecule &mol)
{
  return file.readMolecule(mol);
}

bool readMolecule_2(Helium::MoleculeFile &file, unsigned int index, Molecule &mol)
{
  if (index >= file.numMolecules())
    throw std::runtime_error("Invalid molecule index");
  return file.readMolecule(index, mol);
}

bool readMolecule_3(Helium::MemoryMappedMoleculeFile &file, unsigned int index, Molecule &mol)
{
  if (index >= file.numMolecules())
    throw std::runtime_error("Invalid molecule index");
  return file.readMolecule(index, mol);
}

void export_moleculefile()
{

  class_<Helium::MoleculeOutputFile, boost::noncopyable>("MoleculeOutputFile", init<const std::string&>())
    .def("writeMolecule", &Helium::MoleculeOutputFile::writeMolecule<Molecule>)
    .def("error", &Helium::MoleculeOutputFile::error, return_value_policy<copy_const_reference>())
    ;

  class_<Helium::MoleculeFile, boost::noncopyable>("MoleculeFile")
    .def(init<const std::string&>())
    .def("load", &Helium::MoleculeFile::load)
    .def("numMolecules", &Helium::MoleculeFile::numMolecules)
    .def("readMolecule", &readMolecule_1)
    .def("readMolecule", &readMolecule_2)
    .def("close", &Helium::MoleculeFile::close)
    .def("error", &Helium::MoleculeFile::error, return_internal_reference<>())
    ;

  class_<Helium::MemoryMappedMoleculeFile, boost::noncopyable>("MemoryMappedMoleculeFile")
    .def(init<const std::string&>())
    .def("load", &Helium::MemoryMappedMoleculeFile::load)
    .def("numMolecules", &Helium::MemoryMappedMoleculeFile::numMolecules)
    .def("readMolecule", &readMolecule_3)
    .def("error", &Helium::MemoryMappedMoleculeFile::error, return_internal_reference<>())
    ;
}
