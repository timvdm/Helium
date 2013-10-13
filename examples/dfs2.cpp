// examples/dfs2.cpp
#include <Helium/algorithms/dfs.h>
#include <Helium/hemol.h>

using namespace Helium;

int main()
{
  HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");

  DFSDebugVisitor<HeMol> visitor;
  exhaustive_depth_first_search(mol, get_atom(mol, 0), visitor);
}
