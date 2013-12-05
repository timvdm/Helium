// examples/moleculefile.cpp
#include <Helium/fileio/moleculefile.h> 
#include <Helium/hemol.h> // for HeMol and hemol_from_smiles()
#include <Helium/smiles.h> // for write_smiles()

#include <iostream>

using namespace Helium;

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

    std::cout << write_smiles(mol) << std::endl;
  }
}
