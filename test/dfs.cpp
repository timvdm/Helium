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

int main()
{
  test_dfs1();
  test_dfs2();
  test_dfs3();
}
