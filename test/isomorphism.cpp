#include <Helium/algorithms/isomorphism.h>
#include <Helium/hemol.h>
#include <Helium/smiles.h>

#include "test.h"

using namespace Helium;

void test_isomorphisms(const std::string &smiles)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol;
  SMILES.read(smiles, mol);

  DefaultAtomMatcher<HeMol, HeMol> atomMatcher;
  DefaultBondMatcher<HeMol, HeMol> bondMatcher;
  ASSERT(isomorphism_search(mol, mol, atomMatcher, bondMatcher));
}

bool test_stereo_isomorphism(const std::string &molSmiles, const std::string &querySmiles)
{
  std::cout << "Testing: " << querySmiles << " in " << molSmiles << std::endl;

  HeMol mol, query;
  Stereochemistry molStereo, queryStereo;
  ASSERT(SMILES.read(molSmiles, mol, molStereo));
  ASSERT(SMILES.read(querySmiles, query, queryStereo));

  DefaultAtomMatcher<HeMol, HeMol> atomMatcher;
  DefaultBondMatcher<HeMol, HeMol> bondMatcher;
  return isomorphism_search(mol, query, molStereo, queryStereo, atomMatcher, bondMatcher);
}

int main()
{
  test_isomorphisms("C");
  test_isomorphisms("CCC");
  test_isomorphisms("C1CCCC1");
  test_isomorphisms("CC1C(C)CC(N)C1");

  //
  // tetrahedral stereochemistry
  //
  // no hydrogens
  ASSERT(test_stereo_isomorphism("C[C@](N)(O)S", "C[C@](N)(O)S"));
  ASSERT(test_stereo_isomorphism("CCC[C@](N)(O)S", "C[C@](N)(O)S"));
  ASSERT(!test_stereo_isomorphism("C[C@@](N)(O)S", "C[C@](N)(O)S"));
  ASSERT(!test_stereo_isomorphism("CCC[C@@](N)(O)S", "C[C@](N)(O)S"));
  // implicit hydrogen in query and queried
  ASSERT(test_stereo_isomorphism("C[C@H](N)O", "C[C@H](N)O"));
  ASSERT(test_stereo_isomorphism("C[C@@H](N)O", "C[C@H](O)N"));
  ASSERT(test_stereo_isomorphism("C[C@H](N)O", "C[C@@H](O)N"));
  // implicit hydrogen in query, explicit in queried
  ASSERT(test_stereo_isomorphism("C[C@]([H])(N)O", "C[C@H](N)O"));
  ASSERT(test_stereo_isomorphism("C[C@@]([H])(N)O", "C[C@@H](N)O"));
  ASSERT(!test_stereo_isomorphism("C[C@]([H])(N)O", "C[C@@H](N)O"));
  ASSERT(!test_stereo_isomorphism("C[C@@]([H])(N)O", "C[C@H](N)O"));
  // explicit hydrogen in query and queried
  ASSERT(test_stereo_isomorphism("C[C@]([H])(N)O", "C[C@]([H])(N)O"));
  ASSERT(test_stereo_isomorphism("C[C@@]([H])(N)O", "C[C@@]([H])(N)O"));
  ASSERT(!test_stereo_isomorphism("C[C@]([H])(N)O", "C[C@@]([H])(N)O"));
  ASSERT(!test_stereo_isomorphism("C[C@@]([H])(N)O", "C[C@]([H])(N)O"));

  // allene
  ASSERT(test_stereo_isomorphism("CC(N)=[C@]=C(O)S", "CC(N)=[C@]=C(O)S"));
  ASSERT(!test_stereo_isomorphism("CC(N)=[C@]=C(O)S", "CC(N)=[C@@]=C(O)S"));
  ASSERT(test_stereo_isomorphism("CC(N)=[C@]=C(O)S", "CC(N)=[C@@]=C(S)O"));
  ASSERT(test_stereo_isomorphism("NC(C)=[C@]=C(O)S", "CC(N)=[C@@]=C(O)S"));
  ASSERT(test_stereo_isomorphism("NC(C)=[C@@]=C(O)S", "CC(N)=[C@]=C(O)S"));

  // square-planar stereochemistry
  ASSERT(test_stereo_isomorphism("[Pt@SP1](C)(N)(O)S", "[Pt@SP1](C)(N)(O)S"));
  ASSERT(!test_stereo_isomorphism("[Pt@SP1](C)(N)(O)S", "[Pt@SP2](C)(N)(O)S"));
  ASSERT(!test_stereo_isomorphism("[Pt@SP1](C)(N)(O)S", "[Pt@SP3](C)(N)(O)S"));
  ASSERT(test_stereo_isomorphism("[Pt@SP1](C)(N)(O)S", "[Pt@SP2](C)(O)(N)S"));
  ASSERT(test_stereo_isomorphism("[Pt@SP1](C)(N)(O)S", "[Pt@SP2](S)(N)(O)C"));
  ASSERT(test_stereo_isomorphism("[Pt@SP1](C)(N)(O)S", "[Pt@SP3](N)(C)(O)S"));
  ASSERT(test_stereo_isomorphism("[Pt@SP1](C)(N)(O)S", "[Pt@SP3](C)(N)(S)O"));

  // trigonal-bipyramidal
  ASSERT(test_stereo_isomorphism("[As@](C)(N)(O)(S)P", "[As@](C)(N)(O)(S)P"));
  ASSERT(!test_stereo_isomorphism("[As@](C)(N)(O)(S)P", "[As@@](C)(N)(O)(S)P"));
  ASSERT(test_stereo_isomorphism("[As@](C)(N)(O)(S)P", "[As@@](C)(S)(O)(N)P"));

  // octahedral
  ASSERT(test_stereo_isomorphism("[Co@](C)(N)(O)(S)(P)F", "[Co@](C)(N)(O)(S)(P)F"));
  ASSERT(!test_stereo_isomorphism("[Co@](C)(N)(O)(S)(P)F", "[Co@@](C)(N)(O)(S)(P)F"));
  ASSERT(test_stereo_isomorphism("[Co@](C)(N)(O)(S)(P)F", "[Co@@](C)(P)(S)(O)(N)F"));

  //
  // cis/trans
  //
  // cis + trans with 2x F
  ASSERT(test_stereo_isomorphism("F/C=C/F", "F/C=C/F"));
  ASSERT(test_stereo_isomorphism("F\\C=C/F", "F\\C=C/F"));
  ASSERT(!test_stereo_isomorphism("F/C=C/F", "F/C=C\\F"));
  ASSERT(!test_stereo_isomorphism("F/C=C/F", "F\\C=C/F"));
  // trans == trans with F & I
  ASSERT(test_stereo_isomorphism("F/C=C/I", "F/C=C/I"));
  ASSERT(test_stereo_isomorphism("F/C=C/I", "F\\C=C\\I"));
  ASSERT(test_stereo_isomorphism("F/C=C/I", "I/C=C/F"));
  ASSERT(test_stereo_isomorphism("F/C=C/I", "I\\C=C\\F"));
  // trans != cis
  ASSERT(!test_stereo_isomorphism("F/C=C/I", "F\\C=C/I"));
  ASSERT(!test_stereo_isomorphism("F/C=C/I", "F/C=C\\I"));
  ASSERT(!test_stereo_isomorphism("F/C=C/I", "I\\C=C/F"));
  ASSERT(!test_stereo_isomorphism("F/C=C/I", "I/C=C\\F"));
  // cis == cis
  ASSERT(test_stereo_isomorphism("F\\C=C/I", "F\\C=C/I"));
  ASSERT(test_stereo_isomorphism("F\\C=C/I", "F/C=C\\I"));
  ASSERT(test_stereo_isomorphism("F\\C=C/I", "I\\C=C/F"));
  ASSERT(test_stereo_isomorphism("F\\C=C/I", "I/C=C\\F"));
}
