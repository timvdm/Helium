#include <Helium/algorithms/components.h>
#include <Helium/hemol.h>
#include <Helium/smiles.h>

#include "test.h"

using namespace Helium;

int main()
{
  HeMol mol;
  std::vector<unsigned int> components;

  //
  // connected_bond_components
  //
  SMILES.read("CC(C)C", mol);
  components = connected_bond_components(mol);
  COMPARE(3, components.size());
  COMPARE(1, unique_elements(components));

  SMILES.read("CC.CC", mol);
  components = connected_bond_components(mol);
  COMPARE(2, components.size());
  COMPARE(2, unique_elements(components));
  COMPARE(0, components[0]);
  COMPARE(1, components[1]);

  SMILES.read("C1CC1.CC", mol);
  components = connected_bond_components(mol);
  COMPARE(4, components.size());
  COMPARE(2, unique_elements(components));
  COMPARE(0, components[0]);
  COMPARE(0, components[1]);
  COMPARE(0, components[2]);
  COMPARE(1, components[3]);

  //
  // connected_atom_components
  //
  SMILES.read("CC(C)C", mol);
  components = connected_atom_components(mol);
  COMPARE(4, components.size());
  COMPARE(1, unique_elements(components));

  SMILES.read("CC.CC", mol);
  components = connected_atom_components(mol);
  COMPARE(4, components.size());
  COMPARE(2, unique_elements(components));
  COMPARE(0, components[0]);
  COMPARE(0, components[1]);
  COMPARE(1, components[2]);
  COMPARE(1, components[3]);

  SMILES.read("C1CC1.CC", mol);
  components = connected_atom_components(mol);
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
  SMILES.read("CC(C)C", mol);
  COMPARE(1, num_connected_components(mol));

  SMILES.read("CC.CC", mol);
  COMPARE(2, num_connected_components(mol));

  SMILES.read("C1CC1.CC", mol);
  COMPARE(2, num_connected_components(mol));
}

