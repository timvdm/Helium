#include "../src/molecule.h"
#include "../src/fileio.h"
#include "../src/canonical.h"
#include "../src/extendedconnectivities.h"

#include "test.h"

using namespace Helium;

void testCanonicalPath(const std::string &smiles)
{
  std::cout << "Testing: " << smiles << std::endl;
  Molecule mol;
  read_smiles(smiles, mol);

  std::vector<unsigned int> forward, backward;
  for (unsigned int i = 0; i < num_atoms(&mol); ++i) {
    forward.push_back(i);
    backward.push_back(num_atoms(&mol) - i - 1);
  }

  std::vector<unsigned long> forwardCode = canonicalize_path(&mol, forward).second;
  std::vector<unsigned long> backwardCode = canonicalize_path(&mol, backward).second;

  COMPARE(forwardCode, backwardCode);
}

void testCanonicalize(const std::string &smiles)
{
  std::cout << "Testing: " << smiles << std::endl;
  Molecule mol;
  read_smiles(smiles, mol);

  std::vector<unsigned long> symmetry = extended_connectivities(&mol);
  std::cout << "symmetry: " << symmetry << std::endl;
  canonicalize(&mol, symmetry);
}

int main()
{
  testCanonicalPath("CCCCN");
  testCanonicalPath("CCCCOC");

  testCanonicalize("CCC(C)C");
  testCanonicalize("CCC(C(C)C)C");
  testCanonicalize("c1ccccc1");
}
