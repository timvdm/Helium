// examples/moleculefile.cpp
#include <Helium/fileio/moleculefile.h>
#include <Helium/hemol.h>
#include <Helium/smiles.h>

#include <iostream>

using namespace Helium;

Smiles SMILES;

int main(int argc, char **argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <molecule_file>" << std::endl;
    return -1;
  }

  // perpare molecule file
  MoleculeFile file;
  try {
    file.load(argv[1]);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return -1;
  }

  // read the molecules and print their SMILES
  HeMol mol;
  for (std::size_t i = 0; i < file.numMolecules(); ++i) {
    file.readMolecule(mol);

    std::cout << SMILES.write(mol) << std::endl;
  }
}
