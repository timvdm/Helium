#include <Helium/smiles.h>
#include <Helium/hemol.h>
#include <Helium/fileio/moleculefile.h>
#include <Helium/algorithms/components.h>

#include "test.h"

using namespace Helium;

void test_parse_smiles()
{
  HeMol mol;

  hemol_from_smiles("C", mol);
  COMPARE(1, num_atoms(mol));
  COMPARE(0, num_bonds(mol));
  COMPARE(6, get_element(mol, get_atom(mol, 0)));
  COMPARE(12, get_mass(mol, get_atom(mol, 0)));
  COMPARE(0, get_charge(mol, get_atom(mol, 0)));
  COMPARE(4, get_hydrogens(mol, get_atom(mol, 0)));

  hemol_from_smiles("N", mol);
  COMPARE(1, num_atoms(mol));
  COMPARE(0, num_bonds(mol));
  COMPARE(7, get_element(mol, get_atom(mol, 0)));
  COMPARE(14, get_mass(mol, get_atom(mol, 0)));
  COMPARE(0, get_charge(mol, get_atom(mol, 0)));
  COMPARE(3, get_hydrogens(mol, get_atom(mol, 0)));

  hemol_from_smiles("[13CH4]", mol);
  COMPARE(1, num_atoms(mol));
  COMPARE(0, num_bonds(mol));
  COMPARE(6, get_element(mol, get_atom(mol, 0)));
  COMPARE(13, get_mass(mol, get_atom(mol, 0)));
  COMPARE(0, get_charge(mol, get_atom(mol, 0)));
  COMPARE(4, get_hydrogens(mol, get_atom(mol, 0)));

  hemol_from_smiles("[NH4+]", mol);
  COMPARE(1, num_atoms(mol));
  COMPARE(0, num_bonds(mol));
  COMPARE(7, get_element(mol, get_atom(mol, 0)));
  COMPARE(14, get_mass(mol, get_atom(mol, 0)));
  COMPARE(1, get_charge(mol, get_atom(mol, 0)));
  COMPARE(4, get_hydrogens(mol, get_atom(mol, 0)));

  hemol_from_smiles("CC", mol);
  COMPARE(2, num_atoms(mol));
  COMPARE(1, num_bonds(mol));
  COMPARE(6, get_element(mol, get_atom(mol, 0)));
  COMPARE(6, get_element(mol, get_atom(mol, 1)));
  COMPARE(1, get_order(mol, get_bond(mol, 0)));

  hemol_from_smiles("C=C", mol);
  COMPARE(2, num_atoms(mol));
  COMPARE(1, num_bonds(mol));
  COMPARE(6, get_element(mol, get_atom(mol, 0)));
  COMPARE(6, get_element(mol, get_atom(mol, 1)));
  COMPARE(2, get_order(mol, get_bond(mol, 0)));
  COMPARE(2, get_bond(mol, 0).order());

  hemol_from_smiles("c1ccccc1", mol);
  COMPARE(6, mol.numAtoms());
  COMPARE(6, mol.numBonds());

  hemol_from_smiles("CC.CC", mol);
  COMPARE(4, mol.numAtoms());
  COMPARE(2, mol.numBonds());
  COMPARE(2, num_connected_components(mol));
}

void test_write_smiles(const std::string &smiles)
{
  HeMol mol = hemol_from_smiles(smiles);
  COMPARE(smiles, SMILES.write(mol));
}

void test_tetrahedral(const std::string &smiles, Stereo::Ref center,
    Stereo::Ref ref1, Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4)
{
  std::cout << "Test tetrahedral: " << smiles << std::endl;
  HeMol mol;
  Stereochemistry stereo;
  Smiles SMILES;

  ASSERT(SMILES.read(smiles, mol, stereo));
  COMPARE(stereo.numStereo(), 1);
  COMPARE(stereo.numTetrahedral(), 1);
  ASSERT(stereo.isTetrahedral(center));
  ASSERT(stereo.tetrahedral(center).isValid());
  COMPARE(stereo.tetrahedral(center).type(), Stereo::Tetrahedral);
  COMPARE(stereo.tetrahedral(center).center(), center);
  COMPARE(stereo.tetrahedral(center).ref(0), ref1);
  COMPARE(stereo.tetrahedral(center).ref(1), ref2);
  COMPARE(stereo.tetrahedral(center).ref(2), ref3);
  COMPARE(stereo.tetrahedral(center).ref(3), ref4);
}

void test_tetrahedral()
{
  HeMol mol;
  Stereochemistry stereo;
  Smiles SMILES;

  // - @ @@ @TH1 @TH2
  // - from atom before or after chiral
  // - implicit H
  // - ring closures

  // from atom before chiral
  test_tetrahedral("CC[C@](N)(O)S", 2, 1, 3, 4, 5);
  test_tetrahedral("CC[C@TH1](N)(O)S", 2, 1, 3, 4, 5);
  test_tetrahedral("CC[C@@](N)(O)S", 2, 1, 5, 4, 3);
  test_tetrahedral("CC[C@TH2](N)(O)S", 2, 1, 5, 4, 3);

  // from atom after chiral
  test_tetrahedral("[C@](N)(O)(S)C", 0, 1, 2, 3, 4);
  test_tetrahedral("[C@TH1](N)(O)(S)C", 0, 1, 2, 3, 4);
  test_tetrahedral("[C@@](N)(O)(S)C", 0, 1, 4, 3, 2);
  test_tetrahedral("[C@TH2](N)(O)(S)C", 0, 1, 4, 3, 2);

  // from atom before chiral + implicit H
  test_tetrahedral("CC[C@H](N)O", 2, 1, Stereo::implRef(), 3, 4);
  test_tetrahedral("CC[C@TH1H](N)O", 2, 1, Stereo::implRef(), 3, 4);
  test_tetrahedral("CC[C@@H](N)O", 2, 1, 4, 3, Stereo::implRef());
  test_tetrahedral("CC[C@TH2H](N)O", 2, 1, 4, 3, Stereo::implRef());

  // from atom after chiral + implicit H
  test_tetrahedral("[C@H](N)(O)S", 0, Stereo::implRef(), 1, 2, 3);
  test_tetrahedral("[C@TH1H](N)(O)S", 0, Stereo::implRef(), 1, 2, 3);
  test_tetrahedral("[C@@H](N)(O)S", 0, Stereo::implRef(), 3, 2, 1);
  test_tetrahedral("[C@TH2H](N)(O)S", 0, Stereo::implRef(), 3, 2, 1);

  // ring closure
  test_tetrahedral("[C@]1(Br)(Cl)CCCC(F)C1", 0, 8, 1, 2, 3);
  test_tetrahedral("[C@TH1]1(Br)(Cl)CCCC(F)C1", 0, 8, 1, 2, 3);
  test_tetrahedral("[C@@]1(Br)(Cl)CCCC(F)C1", 0, 8, 3, 2, 1);
  test_tetrahedral("[C@TH2]1(Br)(Cl)CCCC(F)C1", 0, 8, 3, 2, 1);
  test_tetrahedral("[C@H]1(Br)CCCC(F)C1", 0, Stereo::implRef(), 7, 1, 2);
  test_tetrahedral("[C@]1234.C1.N2.O3.S4", 0, 1, 2, 3, 4);
  test_tetrahedral("[C@]1234.C1.N3.O2.S4", 0, 1, 3, 2, 4);
  test_tetrahedral("[C@]1234.C4.N2.O3.S1", 0, 4, 2, 3, 1);
}

void test_roundtrip(const std::string &inputSmiles, std::string expectedSmiles = std::string())
{
  std::cout << "Rountrip: " << inputSmiles << std::endl;
  HeMol mol;
  Stereochemistry stereo;
  Smiles SMILES;

  if (expectedSmiles.empty())
    expectedSmiles = inputSmiles;

  if (!SMILES.read(inputSmiles, mol, stereo)) {
    std::cout << SMILES.error().what() << std::endl;
    bool could_not_read_smiles = false;
    ASSERT(could_not_read_smiles);
  }
  std::string outputSmiles = SMILES.write(mol, stereo);
  COMPARE(expectedSmiles, outputSmiles);
}

void test_roundtrip()
{
  // adding of hydrogens to organic subset
  test_roundtrip("B");
  test_roundtrip("C");
  test_roundtrip("N");
  test_roundtrip("O");
  test_roundtrip("P");
  test_roundtrip("S");
  test_roundtrip("F");
  test_roundtrip("Cl");
  test_roundtrip("Br");
  test_roundtrip("I");
  test_roundtrip("CB");
  test_roundtrip("CC");
  test_roundtrip("CN");
  test_roundtrip("CO");
  test_roundtrip("CP");
  test_roundtrip("CS");
  test_roundtrip("CF");
  test_roundtrip("CCl");
  test_roundtrip("CBr");
  test_roundtrip("CI");

  // simple chains
  test_roundtrip("CCCC");
  test_roundtrip("CC(C)C");
  test_roundtrip("CC(C)C");
  test_roundtrip("CCc1cc(C)ccc1C");
  test_roundtrip("CCC(CC(CC)C)CC");

  // rings
  test_roundtrip("C1CCCCC1");
  test_roundtrip("C1CCCCC1C2CCC2");
  test_roundtrip("C1CCCC2CCCC1CCC2");
  test_roundtrip("c1ccccc1");
  test_roundtrip("[nH]1cccc1");
  test_roundtrip("n1cccc1");

  // bonds
  test_roundtrip("C-C", "CC");
  test_roundtrip("C=C", "C=C");
  test_roundtrip("C#C", "C#C");
  test_roundtrip("C$C", "C$C");
  test_roundtrip("c1:c:c:c:c:c1", "c1ccccc1");
  test_roundtrip("c1ccccc1-c2ccccc2");

  // isotope
  test_roundtrip("[15N]");
  test_roundtrip("[2H]");
  test_roundtrip("[13C]");

  // charge
  test_roundtrip("[C-]");
  test_roundtrip("[N+]");

  // hydrogens
  test_roundtrip("[15NH4+]");
  test_roundtrip("[CH3]");
  test_roundtrip("[CH4]", "C");
  test_roundtrip("[NH3]", "N");

  // various
  test_roundtrip("CCC(=O)C");
  test_roundtrip("CC(N(=O)=O)CC");
  test_roundtrip("CC([N+](=O)[O-])CC");
  test_roundtrip("C(C)(C)(C)C");

  test_roundtrip("F/C=C/I", "FC=CI");
  test_roundtrip("F/C=C/?I", "FC=CI");
  test_roundtrip("F/C=C\\?I", "FC=CI");

  // tetrahedral stereochemistry
  test_roundtrip("C[C@](N)(O)S");
  test_roundtrip("[C@](C)(N)(O)S");
  test_roundtrip("C[C@H](N)O");
  test_roundtrip("[C@H](C)(N)O");
  test_roundtrip("C[C@]123.N1.O2.S3", "C[C@](N)(O)S");
  test_roundtrip("C[C@]123.N2.O1.S3", "C[C@@](N)(O)S");
  test_roundtrip("F[C@]12CCC(C(Cl)C1)C(Br)C2");
  test_roundtrip("F[C@]21CCC(C(Cl)C1)C(Br)C2", "F[C@@]12CCC(C(Cl)C1)C(Br)C2");
  test_roundtrip("F[C@]12CCC(C(Cl)C2)C(Br)C1", "F[C@@]12CCC(C(Cl)C1)C(Br)C2");
  test_roundtrip("F[C@]21CCC(C(Cl)C2)C(Br)C1", "F[C@]12CCC(C(Cl)C1)C(Br)C2");
  test_roundtrip("[C@H]12CCC(C(Cl)C1)C(Br)C2");
  test_roundtrip("[C@H]21CCC(C(Cl)C1)C(Br)C2", "[C@@H]12CCC(C(Cl)C1)C(Br)C2");
  test_roundtrip("[C@H]12CCC(C(Cl)C2)C(Br)C1", "[C@@H]12CCC(C(Cl)C1)C(Br)C2");
  test_roundtrip("[C@H]21CCC(C(Cl)C2)C(Br)C1", "[C@H]12CCC(C(Cl)C1)C(Br)C2");


  //test_roundtrip("");
}

int main()
{
  test_parse_smiles();
  test_tetrahedral();
  test_roundtrip();
}
