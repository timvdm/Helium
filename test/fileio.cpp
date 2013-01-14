#include "../src/fileio.h"

#include "test.h"

using namespace Helium;

int main()
{
  HeMol mol;
  read_smiles("c1ccccc1Cl", mol);

  COMPARE(7, num_atoms(&mol));
  COMPARE(7, num_bonds(&mol));

  COMPARE(6, get_element(&mol, get_atom(&mol, 0)));
}
