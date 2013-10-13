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

int main()
{
  test_bfs1();
}
