#include <Helium/algorithms/canonical.h>
#include <Helium/fileio/moleculefile.h>
#include <Helium/algorithms/extendedconnectivities.h>
#include <Helium/algorithms/components.h>
#include <Helium/smiles.h>


#include "test.h"


using namespace Helium;

#define N 10

void test_canonicalize(const std::string &smiles)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol = hemol_from_smiles(smiles);

  std::vector<unsigned long> symmetry = extended_connectivities(mol, DefaultAtomInvariant());
  std::cout << "symmetry: " << symmetry << std::endl;

  std::vector<unsigned int> atoms = connected_atom_components(mol);
  std::vector<unsigned int> bonds = connected_bond_components(mol);

  canonicalize(mol, symmetry, DefaultAtomInvariant(), DefaultBondInvariant(), atoms, bonds);
}

bool shuffle_test_mol(HeMol &mol)
{
  bool pass = true;
  std::vector<Index> atoms;
  for (std::size_t i = 0; i < num_atoms(mol); ++i)
    atoms.push_back(i);


  std::pair<std::vector<Index>, std::vector<unsigned long> > ref_canon = canonicalize(mol,
        extended_connectivities(mol, DefaultAtomInvariant()), DefaultAtomInvariant(),
        DefaultBondInvariant(), connected_atom_components(mol), connected_bond_components(mol));

  const std::vector<unsigned long> &ref_code = ref_canon.second;
  std::string ref_smiles = SMILES.write(mol, ref_canon.first, Smiles::All);

  // make sure the canonical SMILES is valid and can be read
  try {
  } catch (const Smiley::Exception &e) {
    std::cout << "invalid canonical SMILES:" << std::endl;
    std::cout << e.what() << std::endl;
    ASSERT(0);
    return false;
  }

  for (int i = 0; i < N; ++i) {
    std::cout << "permutation: " << atoms << std::endl;
    std::random_shuffle(atoms.begin(), atoms.end());
    //std::cout << "P: " << atoms << std::endl;
    //std::cout << mol << std::endl;
    mol.renumberAtoms(atoms);
    //std::cout << mol << std::endl;

    std::pair<std::vector<Index>, std::vector<unsigned long> > canon = canonicalize(mol,
        extended_connectivities(mol, DefaultAtomInvariant()), DefaultAtomInvariant(),
        DefaultBondInvariant(), connected_atom_components(mol), connected_bond_components(mol));

    const std::vector<unsigned long> &code = canon.second;
    std::string smiles = SMILES.write(mol, canon.first, Smiles::All);

    std::cout << smiles << std::endl;

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
  HeMol mol = hemol_from_smiles(smiles);
  shuffle_test_mol(mol);
}

void shuffle_test(const std::string &filename)
{
  std::cout << "Shuffle test: " << filename << std::endl;
  MoleculeFile file(filename);

  bool failed = false;
  unsigned int idx = -1, numAtoms = -1;

  HeMol mol;
  for (unsigned int i = 0; i < file.numMolecules(); ++i) {
    file.readMolecule(mol);
    //if (unique_elements(connected_bond_components(mol)) > 1)
    //  continue;
    //if ((i % 1000) == 0)
    //  std::cout << "Testing molecule " << i << "..." << std::endl;

    std::cout << "  testing(" << i << "): " << SMILES.write(mol) << std::endl;

    if (!shuffle_test_mol(mol) && num_atoms(mol) < numAtoms) {
      numAtoms = num_atoms(mol);
      idx = i;
      failed = true;
    }
  }

  if (failed)
    std::cout << "index of smallest molecule that failed: " << idx << std::endl;
}

void test_permutation(const std::string &smiles, const std::vector<Index> &perm)
{
  std::cout << "Testing permuttion: " << smiles << "..." << std::endl;
  std::cout << "    permutation: " << perm << std::endl;
  HeMol mol = hemol_from_smiles(smiles);

  std::pair<std::vector<Index>, std::vector<unsigned long> > ref_canon = canonicalize(mol,
        extended_connectivities(mol, DefaultAtomInvariant()), DefaultAtomInvariant(),
        DefaultBondInvariant(), connected_atom_components(mol), connected_bond_components(mol));

  const std::vector<unsigned long> &ref_code = ref_canon.second;
  std::string ref_smiles = SMILES.write(mol, ref_canon.first, Smiles::All);

  mol.renumberAtoms(perm);

  std::pair<std::vector<Index>, std::vector<unsigned long> > canon = canonicalize(mol,
      extended_connectivities(mol, DefaultAtomInvariant()), DefaultAtomInvariant(),
      DefaultBondInvariant(), connected_atom_components(mol), connected_bond_components(mol));

  const std::vector<unsigned long> &code = canon.second;
  std::string smi = SMILES.write(mol, canon.first, Smiles::All);

  std::cout << smi << std::endl;

  COMPARE(ref_code, code);
  COMPARE(ref_smiles, smi);
}

void test_permutation(const std::string &smiles, const std::string &perm)
{
  std::vector<Index> p;

  std::stringstream ss(perm);
  Index index;
  while ((ss >> index))
    p.push_back(index);

  test_permutation(smiles, p);
}

int main(int argc, char **argv)
{
  bool validate = false;
  if (argc == 2)
    validate = (std::string("-validate") == argv[1]);

  if (argc == 2 && argv[1] == std::string("-bench")) {
    HeMol mol;
    Smiles SMILES;
    SMILES.read("c1(ccc(cc1)/C=N/N(P(=S)(Oc1ccc(cc1)C(N(CP(=O)(O)[O-])CCCC)P(=O)([O-])O)Oc1ccc(cc1)C(P(=O)(O)[O-])N(CP(=O)([O-])O)CCCC)C)OP1(=NP(=NP(=N1)(Oc1ccc(cc1)/C=N/N(P(=S)(Oc1ccc(cc1)C(N(CP(=O)(O)[O-])CCCC)P(=O)([O-])O)Oc1ccc(cc1)C(P(=O)(O)[O-])N(CP(=O)([O-])O)CCCC)C)Oc1ccc(cc1)/C=N/N(P(=S)(Oc1ccc(cc1)C(N(CP(=O)(O)[O-])CCCC)P(=O)([O-])O)Oc1ccc(cc1)C(P(=O)(O)[O-])N(CP(=O)([O-])O)CCCC)C)(Oc1ccc(cc1)/C=N/N(P(=S)(Oc1ccc(cc1)C(N(CP(=O)(O)[O-])CCCC)P(=O)([O-])O)Oc1ccc(cc1)C(P(=O)(O)[O-])N(CP(=O)([O-])O)CCCC)C)Oc1ccc(cc1)/C=N/N(P(=S)(Oc1ccc(cc1)C(N(CP(=O)(O)[O-])CCCC)P(=O)([O-])O)Oc1ccc(cc1)C(P(=O)(O)[O-])N(CP(=O)([O-])O)CCCC)C)Oc1ccc(cc1)/C=N/N(P(=S)(Oc1ccc(cc1)C(N(CP(=O)(O)[O-])CCCC)P(=O)([O-])O)Oc1ccc(cc1)C(P(=O)(O)[O-])N(CP(=O)([O-])O)CCCC)C", mol);
    SMILES.writeCanonical(mol);
    return 0;
  }

  shuffle_test_smiles("CN(P(=N[P+](N=P(N(C)C)(N(C)C)N(C)C)(N=P(N(C)C)(N(C)C)N(C)C)N=P(N(C)C)(N(C)C)N(C)C)(N(C)C)N(C)C)C");


  test_permutation("C[C-]1=[C-][C-]=[C-]([C-]=[C-]1)N.[O+]#[C-].[O+]#[C-].[O+]#[C-].[Cr+9]", "0 8 4 3 13 14 2 7 1 9 10 11 5 12 6");

  shuffle_test_smiles("C[C-]1=[C-][C-]=[C-]([C-]=[C-]1)N.[O+]#[C-].[O+]#[C-].[O+]#[C-].[Cr+9]");



  shuffle_test_smiles("c1ccc(cc1)C[N+]2(CCC3(CC2)OC3)Cc4ccccc4");
  //shuffle_test_smiles("c1ccc(cc1)C[N+]2(CCC3(CC2)OC3)Cc4ccccc4.F[P+](F)(F)(F)(F)F");
  shuffle_test_smiles("[nH]1cccc1");
  shuffle_test_smiles("n1ccccc1");

  shuffle_test_smiles("FC(C1=C[C-]=CC(=C1)C(F)(F)F)(F)F.[Mg+2].[Br-]");

  shuffle_test_smiles("CC[N+]1(CC)CC[N+](CC1)(CC)CC.[Br-].[Br-]");







  shuffle_test_smiles("[Br-].[Br-].[N+]13(CC[N+]2(CC1)CCOCC2)CCOCC3");
  shuffle_test_smiles("[Br-].[Br-].[N+]23(CC[N+]1(CCOCC1)CC2)CCOCC3");
  shuffle_test_smiles("n1ccccc1");
  shuffle_test_smiles("[O-]C(C)=O");

  shuffle_test_smiles("[O-]C(=O)C.[O-]C(=O)C.N.N.N.N.O.O.[Pd+2]");

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

  if (validate) {
    shuffle_test(datadir() + "canonmulti.hel");
    //shuffle_test(datadir() + "1M.hel");
  }
}
