#include "../src/canonical.h"
#include "../src/fileio.h"
#include "../src/extendedconnectivities.h"
#include "../src/components.h"


#include "test.h"

using namespace Helium;

void test_canonicalize(const std::string &smiles)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol;
  read_smiles(smiles, mol);

  std::vector<unsigned long> symmetry = extended_connectivities(mol);
  std::cout << "symmetry: " << symmetry << std::endl;
  canonicalize(mol, symmetry);
}

bool shuffle_test_mol(HeMol &mol)
{
  bool pass = true;
  std::vector<unsigned int> atoms;
  for (unsigned int i = 0; i < num_atoms(mol); ++i)
    atoms.push_back(i);

  std::vector<unsigned long> ref_code = canonicalize(mol, extended_connectivities(mol)).second;

  for (int i = 0; i < 10; ++i) {
    std::random_shuffle(atoms.begin(), atoms.end());
    mol.renumberAtoms(atoms);
    std::vector<unsigned long> code = canonicalize(mol, extended_connectivities(mol)).second;
    COMPARE(ref_code, code);
    if (ref_code != code)
      pass = false;
  }

  return pass;
}

void shuffle_test_smiles(const std::string &smiles)
{
  std::cout << "Testing " << smiles << "..." << std::endl;
  HeMol mol;
  read_smiles(smiles, mol);
  shuffle_test_mol(mol);
}

void shuffle_test(const std::string &filename)
{
  std::cout << "Shuffle test..." << std::endl;
  MoleculeFile file(filename);

  unsigned int idx, numAtoms = -1;

  HeMol mol;
  while (file.read_molecule(mol)) {
    if (unique_elements(connected_bond_components(mol)) > 1)
      continue;
    if ((file.current() % 1000) == 0)
      std::cout << "Testing molecule " << file.current() << "..." << std::endl;
    if (!shuffle_test_mol(mol) && num_atoms(mol) < numAtoms) {
      numAtoms = num_atoms(mol);
      idx = file.current();    
    }
  }
  std::cout << "smallest: " << idx << std::endl;
}

int main()
{
  shuffle_test_smiles("Clc1ccc2c(CCN2C(=O)C)c1");

  test_canonicalize("CCC(C)C");
  test_canonicalize("CCC(C(C)C)C");
  test_canonicalize("c1ccccc1");


  shuffle_test_smiles("[Cl-].OC(=O)C(CS)[NH3+]");
  shuffle_test_smiles("Cl.NCc1ncc(Br)c(C)c1");
  shuffle_test_smiles("C=CCc1c(O)nc(C)nc1O");
  shuffle_test_smiles("Clc1ccc(cc1)Cc1c(C)nc(N)[nH]c1=O");
  shuffle_test_smiles("COc1cc(C)nc(Cl)n1");
  shuffle_test_smiles("C1CN1");
  
  shuffle_test(datadir() + "1K.hem");
}
