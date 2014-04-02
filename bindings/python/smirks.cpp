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

void export_smirks()
{

  class_<Helium::Smirks>("Smirks")
    .def("setFixMass", &Helium::Smirks::setFixMass)
    .def("setFixHydrogens", &Helium::Smirks::setFixHydrogens)
    .def("init", &init_1)
    .def("init", &init_2)
    .def("error", &error)
    .def("requiresCycles", &Helium::Smirks::requiresCycles)
    .def("requiresExplicitHydrogens", &Helium::Smirks::requiresExplicitHydrogens)
    .def("apply", &apply_1)
    .def("apply", &apply_2)
    ;
  
}
