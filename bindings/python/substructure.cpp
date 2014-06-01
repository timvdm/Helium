#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/substructure.h"
#include "common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

Molecule* create_substructure_1(const Molecule &mol, const list &atoms_, const list &bonds_, bool adjustHydrogens)
{
  std::vector<bool> atoms = vector_from_list<bool>(atoms_);
  std::vector<bool> bonds = vector_from_list<bool>(bonds_);

  if (atoms.size() != mol.numAtoms())
    throw std::runtime_error("atoms parameter does not have the correct size");
  if (bonds.size() != mol.numBonds())
    throw std::runtime_error("bonds parameter does not have the correct size");

  std::map<Helium::Index, Helium::Index> indexMap;

  Molecule *sub = new Molecule;

  for (std::size_t i = 0; i < mol.numAtoms(); ++i) {
    if (!atoms[i])
      continue;

    Molecule::atom_type atom = mol.atom(i);
    Molecule::atom_type subAtom = sub->addAtom();

    // count the number of neighbers that are not part of the substructure
    int deltaH = 0;
    if (adjustHydrogens)
      FOREACH_NBR (nbr, atom, mol, Molecule)
        if (!atoms[(*nbr).index()])
          deltaH++;

    subAtom.setAromatic(atom.isAromatic());
    subAtom.setElement(atom.element());
    subAtom.setMass(atom.mass());
    subAtom.setHydrogens(atom.hydrogens() + deltaH);
    subAtom.setCharge(atom.charge());

    indexMap[atom.index()] = subAtom.index();
  }

  for (std::size_t i = 0; i < mol.numBonds(); ++i) {
    if (!bonds[i])
      continue;

    Molecule::bond_type bond = mol.bond(i);
    Molecule::atom_type source = bond.source();
    Molecule::atom_type target = bond.target();

    if (!atoms[source.index()] || !atoms[target.index()])
      continue;

    Molecule::bond_type subBond = sub->addBond(sub->atom(indexMap[source.index()]), sub->atom(indexMap[target.index()]));

    subBond.setAromatic(bond.isAromatic());
    subBond.setOrder(bond.order());
  }

  return sub;
}

Molecule* create_substructure_2(const Molecule &mol, const list &atoms_, const list &bonds_)
{
  return create_substructure_1(mol, atoms_, bonds_, true);
}


void export_substructure()
{

  def("Substructure", &create_substructure_1, return_value_policy<manage_new_object>());
  def("Substructure", &create_substructure_2, return_value_policy<manage_new_object>());

}
