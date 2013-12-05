#include <Helium/algorithms/canonical.h>
#include <Helium/fileio/molecules.h>
#include <Helium/algorithms/extendedconnectivities.h>
#include <Helium/algorithms/components.h>
#include <Helium/smiles.h>


#include "test.h"

using namespace Helium;

void test_canonicalize(const std::string &smiles)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol;
  try {
    parse_smiles(smiles, mol);
  }
  catch (Smiley::Exception &e) {
    std::cerr << e.what();
  }
  std::vector<unsigned long> symmetry = extended_connectivities(mol, AtomInvariant(AtomInvariant::Element));
  std::cout << "symmetry: " << symmetry << std::endl;
  canonicalize(mol, symmetry, AtomInvariant(AtomInvariant::Element), BondInvariant(BondInvariant::Order));
}

bool shuffle_test_mol(HeMol &mol)
{
  bool pass = true;
  std::vector<Index> atoms;
  for (std::size_t i = 0; i < num_atoms(mol); ++i)
    atoms.push_back(i);


  std::pair<std::vector<Index>, std::vector<unsigned long> > ref_canon = canonicalize(mol,
      extended_connectivities(mol, AtomInvariant(AtomInvariant::Element)),
      AtomInvariant(AtomInvariant::Element), BondInvariant(BondInvariant::Order));
  const std::vector<unsigned long> &ref_code = ref_canon.second;
  std::string ref_smiles = write_smiles(mol, ref_canon.first, WriteSmiles::None);

  for (int i = 0; i < 10; ++i) {
    std::random_shuffle(atoms.begin(), atoms.end());
    //std::cout << "P: " << atoms << std::endl;
    //std::cout << mol << std::endl;
    mol.renumberAtoms(atoms);
    //std::cout << mol << std::endl;

    std::pair<std::vector<Index>, std::vector<unsigned long> > canon = canonicalize(mol,
        extended_connectivities(mol, AtomInvariant(AtomInvariant::Element)),
        AtomInvariant(AtomInvariant::Element), BondInvariant(BondInvariant::Order));
    const std::vector<unsigned long> &code = canon.second;
    std::string smiles = write_smiles(mol, canon.first, WriteSmiles::None);

    COMPARE(ref_code, code);
    if (ref_code != code)
      pass = false;

    COMPARE(ref_smiles, smiles);
  }

  return pass;
}

void shuffle_test_smiles(const std::string &smiles)
{
  std::cout << "Testing " << smiles << "..." << std::endl;
  HeMol mol;
  parse_smiles(smiles, mol);
  shuffle_test_mol(mol);
}

void shuffle_test(const std::string &filename)
{
  std::cout << "Shuffle test..." << std::endl;
  MoleculeFile file(filename);

  unsigned int idx, numAtoms = -1;

  HeMol mol;
  for (unsigned int i = 0; i < file.numMolecules(); ++i) {
    file.readMolecule(mol);
    if (unique_elements(connected_bond_components(mol)) > 1)
      continue;
    if ((i % 1000) == 0)
      std::cout << "Testing molecule " << i << "..." << std::endl;
    if (!shuffle_test_mol(mol) && num_atoms(mol) < numAtoms) {
      numAtoms = num_atoms(mol);
      idx = i;
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

  shuffle_test_smiles("CCOC(C)OCC");

  shuffle_test(datadir() + "1K.hel");
}
