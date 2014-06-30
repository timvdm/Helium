#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/kekulize.h"
#include "../../src/algorithms/aromatize.h"
#include "../common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;


bool kekulize_1(Molecule &mol, const Helium::RingSet<Molecule> &rings)
{
  return Helium::kekulize(mol, rings);
}

bool kekulize_2(Molecule &mol)
{
  return Helium::kekulize(mol);
}

bool aromatize_1(Molecule &mol, const Helium::RingSet<Molecule> &rings)
{
  return Helium::aromatize(mol, rings);
}

bool aromatize_2(Molecule &mol)
{
  return Helium::aromatize(mol);
}

void export_aromaticity()
{

  def("kekulize", &kekulize_1);
  def("kekulize", &kekulize_2);
  
  def("aromatize", &aromatize_1);
  def("aromatize", &aromatize_2);

}
