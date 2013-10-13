// examples/components.cpp
#include <Helium/algorithms/components.h>
#include <Helium/hemol.h>

using namespace Helium;

int main()
{
  HeMol mol = hemol_from_smiles("CCC.CC.C");

  Size numComponents = num_connected_components(mol);

  std::cout << "# components: " << numComponents << std::endl;

  std::vector<unsigned int> atomComponentIndices = connected_atom_components(mol);
  std::vector<unsigned int> bondComponentIndices = connected_bond_components(mol);

  for (Size c = 0; c < numComponents; ++c) {
    std::cout << "component " << c << ":" << std::endl;

    // print component atoms
    std::cout << "  atom indices: ";
    for (Index i = 0; i < num_atoms(mol); ++i)
      if (atomComponentIndices[i] == c)
        std::cout << i << " ";
    std::cout << std::endl;

    // print component bonds
    std::cout << "  bond indices: ";
    for (Index i = 0; i < num_bonds(mol); ++i)
      if (bondComponentIndices[i] == c)
        std::cout << i << " ";
    std::cout << std::endl;
  }
}
