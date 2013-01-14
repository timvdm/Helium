#include "../src/isomorphism.h"
#include "../src/fileio.h"

#include "test.h"

using namespace Helium;

void test_isomorphisms(const std::string &smiles)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol;
  read_smiles(smiles, mol);
  isomorphism_search<DefaultAtomMatcher, DefaultBondMatcher, HeMol, HeMol>(mol, mol);
}

int main()
{
  test_isomorphisms("C");
  test_isomorphisms("CCC");
  test_isomorphisms("C1CCCC1");
  test_isomorphisms("CC1C(C)CC(N)C1");
}
