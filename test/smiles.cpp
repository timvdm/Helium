#include "../src/smiles.h"
#include "../src/fileio.h"

#include "test.h"

using namespace Helium;

int main()
{
  HeMol mol;

  parse_smiles("c1ccccc1", mol);
  COMPARE(6, mol.numAtoms());
  COMPARE(6, mol.numBonds());
}
