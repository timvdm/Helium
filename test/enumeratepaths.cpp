#include <Helium/algorithms/enumeratepaths.h>
#include <Helium/smiles.h>

#include "test.h"

using namespace Helium;

void testEnumeratePaths(const std::string &smiles, int expected)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol;
  std::string errorString;
  try {
    parse_smiles(smiles, mol);
  }
  catch (Smiley::Exception &e) {
    std::cerr << e.what();
  }

  std::vector<std::vector<unsigned int> > paths = enumerate_paths(mol, 7);

  COMPARE(expected, paths.size());

  for (std::size_t i = 0; i < paths.size(); ++i)
    std::cout << paths[i] << std::endl;

}

int main()
{
  testEnumeratePaths("CCCCC", 15);
  testEnumeratePaths("CCC(CC)CC", 28);
}
