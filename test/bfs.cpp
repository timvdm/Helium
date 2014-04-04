#include <Helium/algorithms/bfs.h>
#include <Helium/hemol.h>

#include "test.h"

using namespace Helium;

void test_bfs1()
{
  HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");

  std::stringstream output;
  BFSDebugVisitor<HeMol> visitor(output);
  breadth_first_search(mol, visitor);
  compare_file(datadir() + "bfs1.log", output);
}

void test_bfs2()
{
  HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");

  std::vector<bool> atomMask(num_atoms(mol), true);
  atomMask[4] = false;
  atomMask[5] = false;

  std::stringstream output;
  BFSDebugVisitor<HeMol> visitor(output);
  breadth_first_search_mask(mol, visitor, atomMask);
  compare_file(datadir() + "bfs2.log", output);
}

void test_bfs3()
{
  HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");

  std::vector<bool> atomMask(num_atoms(mol), true);
  atomMask[4] = false;
  atomMask[5] = false;
  std::vector<bool> bondMask(num_bonds(mol), true);
  bondMask[3] = false;
  bondMask[4] = false;
  bondMask[5] = false;

  std::stringstream output;
  BFSDebugVisitor<HeMol> visitor(output);
  breadth_first_search_mask(mol, visitor, atomMask, bondMask);
  compare_file(datadir() + "bfs3.log", output);
}

void test_bfs4()
{
  HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");

  std::stringstream output;
  BFSDebugVisitor<HeMol> visitor(output);
  breadth_first_search(mol, mol.atom(5), visitor);
  compare_file(datadir() + "bfs4.log", output);
}

void test_bfs5()
{
  HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");

  std::vector<bool> atomMask(num_atoms(mol), true);
  atomMask[0] = false;

  std::stringstream output;
  BFSDebugVisitor<HeMol> visitor(output);
  breadth_first_search_mask(mol, mol.atom(5), visitor, atomMask);
  compare_file(datadir() + "bfs5.log", output);
}

void test_bfs6()
{
  HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");

  std::vector<bool> atomMask(num_atoms(mol), true);
  atomMask[0] = false;
  std::vector<bool> bondMask(num_bonds(mol), true);
  bondMask[0] = false;
  bondMask[6] = false;
  bondMask[7] = false;

  std::stringstream output;
  BFSDebugVisitor<HeMol> visitor(output);
  breadth_first_search_mask(mol, mol.atom(5), visitor, atomMask, bondMask);
  compare_file(datadir() + "bfs6.log", output);
}

int main()
{
  test_bfs1();
  test_bfs2();
  test_bfs3();
  test_bfs4();
  test_bfs5();
  test_bfs6();
}
