#include "../src/smiles.h"
#include "../src/fileio/molecules.h"

#include "test.h"

using namespace Helium;

int main()
{
  HeMol mol;

  try {
    parse_smiles("c1ccccc1", mol);
  }
  catch(Smiley::Exception &e) {
    std::cerr << e.what();
  }
  COMPARE(6, mol.numAtoms());
  COMPARE(6, mol.numBonds());
}
