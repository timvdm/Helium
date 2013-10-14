#include "../src/smiles.h"
#include "../src/fileio/molecules.h"

#include "test.h"

using namespace Helium;

int main()
{
  HeMol mol;

  hemol_from_smiles("CC", mol);
  COMPARE(2, num_atoms(mol));
  COMPARE(1, num_bonds(mol));
  COMPARE(6, get_element(mol, get_atom(mol, 0)));
  COMPARE(6, get_element(mol, get_atom(mol, 1)));
  COMPARE(1, get_order(mol, get_bond(mol, 0)));

  hemol_from_smiles("C=C", mol);
  COMPARE(2, num_atoms(mol));
  COMPARE(1, num_bonds(mol));
  COMPARE(6, get_element(mol, get_atom(mol, 0)));
  COMPARE(6, get_element(mol, get_atom(mol, 1)));
  COMPARE(2, get_order(mol, get_bond(mol, 0)));
  COMPARE(2, get_bond(mol, 0).order());

  hemol_from_smiles("c1ccccc1", mol);
  COMPARE(6, mol.numAtoms());
  COMPARE(6, mol.numBonds());
}
