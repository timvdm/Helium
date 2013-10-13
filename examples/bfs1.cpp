// examples/bfs1.cpp
#include <Helium/algorithms/bfs.h>
#include <Helium/hemol.h>

using namespace Helium;

int main()
{
  HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");

  BFSDebugVisitor<HeMol> visitor;
  breadth_first_search(mol, visitor);
}
