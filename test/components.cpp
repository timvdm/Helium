#include "../src/molecule.h"
#include "../src/fileio.h"
#include "../src/components.h"

#include "test.h"

using namespace Helium;

int main()
{
  Molecule mol;
  std::vector<unsigned int> components;

  read_smiles("CC(C)C", mol);
  components = connected_bond_components(&mol);
  COMPARE(3, components.size());
  COMPARE(1, unique_elements(components));

  read_smiles("CC.CC", mol);
  components = connected_bond_components(&mol);
  COMPARE(2, components.size());
  COMPARE(2, unique_elements(components));
  COMPARE(0, components[0]);
  COMPARE(1, components[1]);
}
