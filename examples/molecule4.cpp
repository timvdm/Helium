// examples/molecule4.cpp
#include <Helium/molecule.h>
#include <Helium/hemol.h>
#include <Helium/smiles.h>
#include <Helium/element.h>

#include <iostream>

using namespace Helium;

Smiles SMILES;

int main()
{
  typedef molecule_traits<HeMol>::atom_type atom_type;
  typedef molecule_traits<HeMol>::bond_type bond_type;

  // create acetone molecule
  HeMol mol;

  // carbon atom 1
  atom_type C1 = add_atom(mol);
  set_element(mol, C1, 6);
  set_mass(mol, C1, Element::averageMass(6));

  // carbon atom 2
  atom_type C2 = add_atom(mol);
  set_element(mol, C2, 6);
  set_mass(mol, C2, Element::averageMass(6));

  // carbon atom 3
  atom_type C3 = add_atom(mol);
  set_element(mol, C3, 6);
  set_mass(mol, C3, Element::averageMass(6));

  // oxygen atom
  atom_type O = add_atom(mol);
  set_element(mol, O, 8);
  set_mass(mol, O, Element::averageMass(8));

  // C1-C2 bond
  add_bond(mol, C1, C2);

  // C2-C3 bond
  add_bond(mol, C2, C3);

  // C2=O bond
  bond_type bond = add_bond(mol, C2, O);
  set_order(mol, bond, 2);

  std::cout << "SMILES: " << SMILES.write(mol) << std::endl;
}
