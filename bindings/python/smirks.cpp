#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/smirks.h"

using Helium::Chemist::Molecule;
using namespace boost::python;
      
bool init_1(Helium::Smirks &SMIRKS, const std::string &smirks)
{
  return SMIRKS.init(smirks);
}

bool init_2(Helium::Smirks &SMIRKS, const std::string &reactant, const std::string &product)
{
  return SMIRKS.init(reactant, product);
}

Helium::Error error(const Helium::Smirks &smirks)
{
  if (smirks.error().type() == Helium::SmirksError::None)
    return Helium::Error();
  return Helium::Error(smirks.error().what());
}

bool apply_1(Helium::Smirks &smirks, Molecule &mol, const Helium::RingSet<Molecule> &rings)
{
  return smirks.apply(mol, rings);
}

bool apply_2(Helium::Smirks &smirks, Molecule &mol)
{
  return smirks.apply(mol, Helium::RingSet<Molecule>(mol));
}

std::vector<boost::shared_ptr<Molecule> > react_1(Helium::Smirks &smirks, const Molecule &mol, const Helium::RingSet<Molecule> &rings, int min, int max)
{
  return smirks.react(mol, rings, min, max);
}

std::vector<boost::shared_ptr<Molecule> > react_2(Helium::Smirks &smirks, const Molecule &mol, const Helium::RingSet<Molecule> &rings, int min)
{
  return smirks.react(mol, rings, min);
}

std::vector<boost::shared_ptr<Molecule> > react_3(Helium::Smirks &smirks, const Molecule &mol, const Helium::RingSet<Molecule> &rings)
{
  return smirks.react(mol, rings);
}

std::vector<boost::shared_ptr<Molecule> > react_4(Helium::Smirks &smirks, const Molecule &mol, int min, int max)
{
  return smirks.react(mol, min, max);
}

std::vector<boost::shared_ptr<Molecule> > react_5(Helium::Smirks &smirks, const Molecule &mol, int min)
{
  return smirks.react(mol, min);
}

std::vector<boost::shared_ptr<Molecule> > react_6(Helium::Smirks &smirks, const Molecule &mol)
{
  return smirks.react(mol);
}

void export_smirks()
{

  class_<Helium::Smirks>("Smirks")
    .def("setFixMass", &Helium::Smirks::setFixMass)
    .def("setFixHydrogens", &Helium::Smirks::setFixHydrogens)
    .def("init", &init_1)
    .def("init", &init_2)
    .def("error", &error)
    .def("requiresRingSet", &Helium::Smirks::requiresRingSet)
    .def("requiresCyclicity", &Helium::Smirks::requiresCyclicity)
    .def("requiresExplicitHydrogens", &Helium::Smirks::requiresExplicitHydrogens)
    .def("apply", &apply_1)
    .def("apply", &apply_2)
    .def("react", &react_1)
    .def("react", &react_2)
    .def("react", &react_3)
    .def("react", &react_4)
    .def("react", &react_5)
    .def("react", &react_6)
    ;
  
}
