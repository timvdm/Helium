#include "../src/components.h"
#include "../src/fileio.h"

#include "test.h"

using namespace Helium;

int main()
{
  HeMol mol;
  std::vector<unsigned int> components;

  //
  // connected_bond_components
  //
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

  read_smiles("C1CC1.CC", mol);
  components = connected_bond_components(&mol);
  COMPARE(4, components.size());
  COMPARE(2, unique_elements(components));
  COMPARE(0, components[0]);
  COMPARE(0, components[1]);
  COMPARE(0, components[2]);
  COMPARE(1, components[3]);

  //
  // connected_atom_components
  //
  read_smiles("CC(C)C", mol);
  components = connected_atom_components(&mol);
  COMPARE(4, components.size());
  COMPARE(1, unique_elements(components));

  read_smiles("CC.CC", mol);
  components = connected_atom_components(&mol);
  COMPARE(4, components.size());
  COMPARE(2, unique_elements(components));
  COMPARE(0, components[0]);
  COMPARE(0, components[1]);
  COMPARE(1, components[2]);
  COMPARE(1, components[3]);

  read_smiles("C1CC1.CC", mol);
  components = connected_atom_components(&mol);
  COMPARE(5, components.size());
  COMPARE(2, unique_elements(components));
  COMPARE(0, components[0]);
  COMPARE(0, components[1]);
  COMPARE(0, components[2]);
  COMPARE(1, components[3]);
  COMPARE(1, components[4]);

  //
  // num_connected_components
  //
  read_smiles("CC(C)C", mol);
  COMPARE(1, num_connected_components(&mol));

  read_smiles("CC.CC", mol);
  COMPARE(2, num_connected_components(&mol));

  read_smiles("C1CC1.CC", mol);
  COMPARE(2, num_connected_components(&mol));
}

