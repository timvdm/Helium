#include "../src/molecule.h"
#include "../src/fileio.h"
#include "../src/enumeratepaths.h"

#include "test.h"

using namespace Helium;

void testEnumeratePaths(const std::string &smiles, int expected)
{
  std::cout << "Testing: " << smiles << std::endl;
  Molecule mol;
  read_smiles(smiles, mol);

  std::vector<std::vector<unsigned int> > paths = enumerate_paths(&mol, 7);

  COMPARE(expected, paths.size());

  for (std::size_t i = 0; i < paths.size(); ++i)
    std::cout << paths[i] << std::endl;

}

int main()
{
  testEnumeratePaths("CCCCC", 15);
  testEnumeratePaths("CCC(CC)CC", 28);
}
