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
  return file.readMolecule(index, mol);
}

void export_moleculefile()
{

  class_<Helium::MoleculeFile, boost::noncopyable>("MoleculeFile")
    .def(init<const std::string&>())
    .def("load", &Helium::MoleculeFile::load)
    .def("numMolecules", &Helium::MoleculeFile::numMolecules)
    .def("readMolecule", &readMolecule_1)
    .def("readMolecule", &readMolecule_2)
    .def("close", &Helium::MoleculeFile::close)
    ;

  class_<Helium::MemoryMappedMoleculeFile, boost::noncopyable>("MemoryMappedMoleculeFile")
    .def(init<const std::string&>())
    .def("load", &Helium::MemoryMappedMoleculeFile::load)
    .def("numMolecules", &Helium::MemoryMappedMoleculeFile::numMolecules)
    .def("readMolecule", &Helium::MemoryMappedMoleculeFile::readMolecule<Molecule>)
    ;
}
