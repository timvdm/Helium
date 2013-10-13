#include <Helium/algorithms/isomorphism.h>
#include <Helium/smiles.h>

#include "test.h"

using namespace Helium;

void test_isomorphisms(const std::string &smiles)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol;
  try {
    parse_smiles(smiles, mol);
  }
  catch (Smiley::Exception &e) {
    std::cerr << e.what();
  }

  DefaultAtomMatcher<HeMol, HeMol> atomMatcher;
  DefaultBondMatcher<HeMol, HeMol> bondMatcher;
  isomorphism_search(mol, mol, atomMatcher, bondMatcher);
}

int main()
{
  test_isomorphisms("C");
  test_isomorphisms("CCC");
  test_isomorphisms("C1CCCC1");
  test_isomorphisms("CC1C(C)CC(N)C1");
}
