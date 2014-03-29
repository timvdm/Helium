#include <Helium/algorithms/dfs.h>
#include <Helium/hemol.h>

#include "test.h"

using namespace Helium;

void test_dfs1()
{
  HeMol mol = hemol_from_smiles("C(CC)(CC)CC");

  std::stringstream output;
  DFSDebugVisitor<HeMol> visitor(output);
  depth_first_search(mol, visitor);
  compare_file(datadir() + "dfs1.log", output);
}

void test_dfs2()
{
  HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");

  std::stringstream output;
  DFSDebugVisitor<HeMol> visitor(output);
  depth_first_search(mol, visitor);
  compare_file(datadir() + "dfs2.log", output);
}

void test_dfs3()
{
  HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");

  std::stringstream output;
  DFSDebugVisitor<HeMol> visitor(output);
  exhaustive_depth_first_search(mol, get_atom(mol, 0), visitor);
  compare_file(datadir() + "dfs3.log", output);
}

// basic dfs (all components)
void test_dfs4()
{
  HeMol mol = hemol_from_smiles("C1CCCCC1CC.CC");

  std::stringstream output;
  DFSDebugVisitor<HeMol> visitor(output);
  depth_first_search(mol, visitor);
  compare_file(datadir() + "dfs4.log", output);
}

// basic dfs with atom mask (all components)
void test_dfs5()
{
  HeMol mol = hemol_from_smiles("C1CCCCC1CC.CC");

  std::vector<bool> atomMask(num_atoms(mol), true);
  atomMask[6] = false;
  atomMask[7] = false;

  ASSERT(is_valid_atom_mask(mol, atomMask));

  std::stringstream output;
  DFSDebugVisitor<HeMol> visitor(output);
  depth_first_search_mask(mol, visitor, atomMask);
  compare_file(datadir() + "dfs5.log", output);
}

// basic dfs with atom and bond mask (all components)
void test_dfs6()
{
  HeMol mol = hemol_from_smiles("C1CCCCC1CC.CC");

  std::vector<bool> atomMask(num_atoms(mol), true);
  atomMask[6] = false;
  atomMask[7] = false;
  std::vector<bool> bondMask(num_bonds(mol), true);
  bondMask[6] = false;
  bondMask[7] = false;
  bondMask[8] = false;

  ASSERT(is_valid_bond_mask(mol, atomMask, bondMask));

  std::stringstream output;
  DFSDebugVisitor<HeMol> visitor(output);
  depth_first_search_mask(mol, visitor, atomMask, bondMask);
  compare_file(datadir() + "dfs6.log", output);
}

// basic dfs starting at specific atom (single component)
void test_dfs7()
{
  HeMol mol = hemol_from_smiles("C1CCCCC1CC.CC");

  std::stringstream output;
  DFSDebugVisitor<HeMol> visitor(output);
  depth_first_search(mol, get_atom(mol, 0), visitor);
  compare_file(datadir() + "dfs7.log", output);
}

// basic dfs starting at specific atom with atom mask (single component)
void test_dfs8()
{
  HeMol mol = hemol_from_smiles("C1CCCCC1CC.CC");

  std::vector<bool> atomMask(num_atoms(mol), true);
  atomMask[6] = false;
  atomMask[7] = false;

  ASSERT(is_valid_atom_mask(mol, atomMask));

  std::stringstream output;
  DFSDebugVisitor<HeMol> visitor(output);
  depth_first_search_mask(mol, get_atom(mol, 0), visitor, atomMask);
  compare_file(datadir() + "dfs8.log", output);
}

// basic dfs starting at specific atom with atom and bond mask (single component)
void test_dfs9()
{
  HeMol mol = hemol_from_smiles("C1CCCCC1CC.CC");

  std::vector<bool> atomMask(num_atoms(mol), true);
  atomMask[6] = false;
  atomMask[7] = false;
  std::vector<bool> bondMask(num_bonds(mol), true);
  bondMask[5] = false;
  bondMask[6] = false;
  bondMask[7] = false;

  ASSERT(is_valid_bond_mask(mol, atomMask, bondMask));

  std::stringstream output;
  DFSDebugVisitor<HeMol> visitor(output);
  depth_first_search_mask(mol, get_atom(mol, 0), visitor, atomMask, bondMask);
  compare_file(datadir() + "dfs9.log", output);
}

// dfs with specified order (all component)

// exhaustive dfs starting at specific atom (single component)

int main()
{
  test_dfs1();
  test_dfs2();
  test_dfs3();
  test_dfs4();
  test_dfs5();
  test_dfs6();
  test_dfs7();
  test_dfs8();
  test_dfs9();

}
