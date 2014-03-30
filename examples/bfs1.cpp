// examples/bfs1.cpp
#include <Helium/algorithms/bfs.h>
#include <Helium/hemol.h>
#include <Helium/smiles.h>

using namespace Helium;

int main()
{
  HeMol mol;
  Smiles SMILES;
  if (!SMILES.read("C1CCC(CC)CC1", mol)) {
    std::cerr << SMILES.error().what();
    return -1;
  }

  BFSDebugVisitor<HeMol> visitor;
  breadth_first_search(mol, visitor);
}
