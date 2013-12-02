#include "../src/smiles.h"
#include "../src/fileio/molecules.h"
#include "../src/algorithms/components.h"

#include "test.h"

using namespace Helium;

void test_parse_smiles()
{
  HeMol mol;

  hemol_from_smiles("C", mol);
  COMPARE(1, num_atoms(mol));
  COMPARE(0, num_bonds(mol));
  COMPARE(6, get_element(mol, get_atom(mol, 0)));
  COMPARE(12, get_mass(mol, get_atom(mol, 0)));
  COMPARE(0, get_charge(mol, get_atom(mol, 0)));
  COMPARE(4, num_hydrogens(mol, get_atom(mol, 0)));

  hemol_from_smiles("N", mol);
  COMPARE(1, num_atoms(mol));
  COMPARE(0, num_bonds(mol));
  COMPARE(7, get_element(mol, get_atom(mol, 0)));
  COMPARE(14, get_mass(mol, get_atom(mol, 0)));
  COMPARE(0, get_charge(mol, get_atom(mol, 0)));
  COMPARE(3, num_hydrogens(mol, get_atom(mol, 0)));

  hemol_from_smiles("[13CH4]", mol);
  COMPARE(1, num_atoms(mol));
  COMPARE(0, num_bonds(mol));
  COMPARE(6, get_element(mol, get_atom(mol, 0)));
  COMPARE(13, get_mass(mol, get_atom(mol, 0)));
  COMPARE(0, get_charge(mol, get_atom(mol, 0)));
  COMPARE(4, num_hydrogens(mol, get_atom(mol, 0)));

  hemol_from_smiles("[NH4+]", mol);
  COMPARE(1, num_atoms(mol));
  COMPARE(0, num_bonds(mol));
  COMPARE(7, get_element(mol, get_atom(mol, 0)));
  COMPARE(14, get_mass(mol, get_atom(mol, 0)));
  COMPARE(1, get_charge(mol, get_atom(mol, 0)));
  COMPARE(4, num_hydrogens(mol, get_atom(mol, 0)));

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

  hemol_from_smiles("CC.CC", mol);
  COMPARE(4, mol.numAtoms());
  COMPARE(2, mol.numBonds());
  COMPARE(2, num_connected_components(mol));
}

void test_write_smiles(const std::string &smiles)
{
  HeMol mol = hemol_from_smiles(smiles);
  COMPARE(smiles, write_smiles(mol));
}

int main()
{
  test_parse_smiles();

  // simple chain
  test_write_smiles("CCCCC");
  // branches
  test_write_smiles("CCC(CC)CC");
  test_write_smiles("CCC(CC)(CC)CC");
  test_write_smiles("CCC(CC(CC)CC)CC");
  // rings
  test_write_smiles("C1CCCC1");
  test_write_smiles("C1CCCC1C2CCCC2");
  // aromatic
  test_write_smiles("n1ccccc1");

  // isotope
  test_write_smiles("[13C]");
  // charge
  test_write_smiles("[C+]");
  test_write_smiles("[C-]");

  test_write_smiles("C=C");
  test_write_smiles("C#C");
  test_write_smiles("C$C");

  test_write_smiles("CCC(=O)C");
  test_write_smiles("CC(N(=O)=O)CC");
  test_write_smiles("CC([N+](=O)[O-])CC");
}
