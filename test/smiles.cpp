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

void test_allene(const std::string &smiles, Stereo::Ref center,
    Stereo::Ref ref1, Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4)
{
  std::cout << "Test allene: " << smiles << std::endl;
  HeMol mol;
  Stereochemistry stereo;
  Smiles SMILES;

  ASSERT(SMILES.read(smiles, mol, stereo));
  COMPARE(stereo.numStereo(), 1);
  COMPARE(stereo.numAllene(), 1);
  ASSERT(stereo.isAllene(center));
  ASSERT(stereo.allene(center).isValid());
  COMPARE(stereo.allene(center).type(), Stereo::Allene);
  COMPARE(stereo.allene(center).center(), center);
  COMPARE(stereo.allene(center).ref(0), ref1);
  COMPARE(stereo.allene(center).ref(1), ref2);
  COMPARE(stereo.allene(center).ref(2), ref3);
  COMPARE(stereo.allene(center).ref(3), ref4);
}

void test_allene()
{
  test_allene("CC(N)=[C@]=C(O)S", 3, 0, 2, 5, 6);
  test_allene("CC(N)=[C@@]=C(O)S", 3, 0, 6, 5, 2);
  test_allene("CC(N)=[C@AL1]=C(O)S", 3, 0, 2, 5, 6);
  test_allene("CC(N)=[C@AL2]=C(O)S", 3, 0, 6, 5, 2);
  test_allene("C(C)(N)=[C@AL1]=C(O)S", 3, 1, 2, 5, 6);
  test_allene("C(C)(N)=[C@AL2]=C(O)S", 3, 1, 6, 5, 2);
}

void test_cistrans(const std::string &smiles, Stereo::Ref center,
    Stereo::Ref ref1, Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4)
{
  std::cout << "Test cis/trans: " << smiles << std::endl;
  HeMol mol;
  Stereochemistry stereo;
  Smiles SMILES;

  ASSERT(SMILES.read(smiles, mol, stereo));
  COMPARE(stereo.numStereo(), 1);
  COMPARE(stereo.numCisTrans(), 1);
  ASSERT(stereo.isCisTrans(center));
  ASSERT(stereo.cisTrans(center).isValid());
  COMPARE(stereo.cisTrans(center).type(), Stereo::CisTrans);
  COMPARE(stereo.cisTrans(center).center(), center);
  COMPARE(stereo.cisTrans(center).ref(0), ref1);
  COMPARE(stereo.cisTrans(center).ref(1), ref2);
  COMPARE(stereo.cisTrans(center).ref(2), ref3);
  COMPARE(stereo.cisTrans(center).ref(3), ref4);
}

void test_cistrans()
{
  test_cistrans("F/C=C/F", 1, Stereo::implRef(), 0, Stereo::implRef(), 3);
  test_cistrans("F\\C=C/F", 1, 0, Stereo::implRef(), Stereo::implRef(), 3);
  test_cistrans("F/C=C\\F", 1, Stereo::implRef(), 0, 3, Stereo::implRef());
  test_cistrans("F\\C=C\\F", 1, 0, Stereo::implRef(), 3, Stereo::implRef());

  test_cistrans("C(\\F)=C/F", 1, Stereo::implRef(), 1, Stereo::implRef(), 3);
  test_cistrans("C(/F)=C/F", 1, 1, Stereo::implRef(), Stereo::implRef(), 3);
  test_cistrans("C(\\F)=C\\F", 1, Stereo::implRef(), 1, 3, Stereo::implRef());
  test_cistrans("C(/F)=C\\F", 1, 1, Stereo::implRef(), 3, Stereo::implRef());

  test_cistrans("F/C(/I)=C(/F)\\I", 2, 2, 0, 5, 4);
  test_cistrans("F\\C(\\I)=C(/F)\\I", 2, 0, 2, 5, 4);
  test_cistrans("F/C(/I)=C(\\F)/I", 2, 2, 0, 4, 5);
  test_cistrans("F\\C(\\I)=C(\\F)/I", 2, 0, 2, 4, 5);

  test_cistrans("F/C=C/1.F1", 1, Stereo::implRef(), 0, Stereo::implRef(), 3);
  test_cistrans("F/C=C1.F/1", 1, Stereo::implRef(), 0, Stereo::implRef(), 3);
  test_cistrans("F1.C/1=C/F", 1, 0, Stereo::implRef(), Stereo::implRef(), 3);
  test_cistrans("F1.C/1=C/1.F1", 1, 0, Stereo::implRef(), Stereo::implRef(), 3);
}

void test_squareplanar(const std::string &smiles, Stereo::Ref center,
    Stereo::Ref ref1, Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4)
{
  std::cout << "Test square planar: " << smiles << std::endl;
  HeMol mol;
  Stereochemistry stereo;
  Smiles SMILES;

  ASSERT(SMILES.read(smiles, mol, stereo));
  COMPARE(stereo.numStereo(), 1);
  COMPARE(stereo.numSquarePlanar(), 1);
  ASSERT(stereo.isSquarePlanar(center));
  ASSERT(stereo.squarePlanar(center).isValid());
  COMPARE(stereo.squarePlanar(center).type(), Stereo::SquarePlanar);
  COMPARE(stereo.squarePlanar(center).center(), center);
  COMPARE(stereo.squarePlanar(center).ref(0), ref1);
  COMPARE(stereo.squarePlanar(center).ref(1), ref2);
  COMPARE(stereo.squarePlanar(center).ref(2), ref3);
  COMPARE(stereo.squarePlanar(center).ref(3), ref4);
}

void test_squareplanar()
{
  test_squareplanar("C[Pt@SP1](N)(O)S", 1, 0, 2, 3, 4);
  test_squareplanar("[Pt@SP1](C)(N)(O)S", 0, 1, 2, 3, 4);
  test_squareplanar("C[Pt@SP1H](N)O", 1, 0, Stereo::implRef(), 2, 3);
  test_squareplanar("[Pt@SP1H](C)(N)O", 0, Stereo::implRef(), 1, 2, 3);
  test_squareplanar("C[Pt@SP1H2]N", 1, 0, Stereo::implRef(), Stereo::implRef(), 2);
  test_squareplanar("[Pt@SP1H2](C)N", 0, Stereo::implRef(), Stereo::implRef(), 1, 2);

  test_squareplanar("C[Pt@SP2](N)(O)S", 1, 0, 3, 2, 4);
  test_squareplanar("[Pt@SP2](C)(N)(O)S", 0, 1, 3, 2, 4);
  test_squareplanar("C[Pt@SP2H](N)O", 1, 0, 2, Stereo::implRef(), 3);
  test_squareplanar("[Pt@SP2H](C)(N)O", 0, Stereo::implRef(), 2, 1, 3);
  test_squareplanar("C[Pt@SP2H2]N", 1, 0, Stereo::implRef(), Stereo::implRef(), 2);
  test_squareplanar("[Pt@SP2H2](C)N", 0, Stereo::implRef(), 1, Stereo::implRef(), 2);

  test_squareplanar("C[Pt@SP3](N)(O)S", 1, 2, 0, 3, 4);
  test_squareplanar("[Pt@SP3](C)(N)(O)S", 0, 2, 1, 3, 4);
  test_squareplanar("C[Pt@SP3H](N)O", 1, Stereo::implRef(), 0, 2, 3);
  test_squareplanar("[Pt@SP3H](C)(N)O", 0, 1, Stereo::implRef(), 2, 3);
  test_squareplanar("C[Pt@SP3H2]N", 1, Stereo::implRef(), 0, Stereo::implRef(), 2);
  test_squareplanar("[Pt@SP3H2](C)N", 0, Stereo::implRef(), Stereo::implRef(), 1, 2);
}

void test_trigonalbipyramidal(const std::string &smiles, Stereo::Ref center,
    Stereo::Ref ref1, Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4,
    Stereo::Ref ref5)
{
  std::cout << "Test trigonal bipyramidal: " << smiles << std::endl;
  HeMol mol;
  Stereochemistry stereo;
  Smiles SMILES;

  ASSERT(SMILES.read(smiles, mol, stereo));
  COMPARE(stereo.numStereo(), 1);
  COMPARE(stereo.numTrigonalBipyramidal(), 1);
  ASSERT(stereo.isTrigonalBipyramidal(center));
  ASSERT(stereo.trigonalBipyramidal(center).isValid());
  //std::cout << stereo.trigonalBipyramidal(center) << std::endl;
  COMPARE(stereo.trigonalBipyramidal(center).type(), Stereo::TrigonalBipyramidal);
  COMPARE(stereo.trigonalBipyramidal(center).center(), center);
  COMPARE(stereo.trigonalBipyramidal(center).ref(0), ref1);
  COMPARE(stereo.trigonalBipyramidal(center).ref(1), ref2);
  COMPARE(stereo.trigonalBipyramidal(center).ref(2), ref3);
  COMPARE(stereo.trigonalBipyramidal(center).ref(3), ref4);
  COMPARE(stereo.trigonalBipyramidal(center).ref(4), ref5);
}

void test_trigonalbipyramidal()
{
  // TB1/TB2
  test_trigonalbipyramidal("C[As@](N)(O)(S)P", 1, 0, 2, 3, 4, 5);
  test_trigonalbipyramidal("C[As@TB1](N)(O)(S)P", 1, 0, 2, 3, 4, 5);
  test_trigonalbipyramidal("C[As@@](N)(O)(S)P", 1, 5, 2, 3, 4, 0);
  test_trigonalbipyramidal("C[As@TB2](N)(O)(S)P", 1, 5, 2, 3, 4, 0);

  test_trigonalbipyramidal("[As@](C)(N)(O)(S)P", 0, 1, 2, 3, 4, 5);
  test_trigonalbipyramidal("[As@TB1](C)(N)(O)(S)P", 0, 1, 2, 3, 4, 5);
  test_trigonalbipyramidal("[As@@](C)(N)(O)(S)P", 0, 5, 2, 3, 4, 1);
  test_trigonalbipyramidal("[As@TB2](C)(N)(O)(S)P", 0, 5, 2, 3, 4, 1);

  test_trigonalbipyramidal("C[As@H](N)(O)S", 1, 0, Stereo::implRef(), 2, 3, 4);
  test_trigonalbipyramidal("C[As@TB1H](N)(O)S", 1, 0, Stereo::implRef(), 2, 3, 4);
  test_trigonalbipyramidal("C[As@@H](N)(O)S", 1, 4, Stereo::implRef(), 2, 3, 0);
  test_trigonalbipyramidal("C[As@TB2H](N)(O)S", 1, 4, Stereo::implRef(), 2, 3, 0);

  test_trigonalbipyramidal("[As@H](C)(N)(O)S", 0, Stereo::implRef(), 1, 2, 3, 4);
  test_trigonalbipyramidal("[As@TB1H](C)(N)(O)S", 0, Stereo::implRef(), 1, 2, 3, 4);
  test_trigonalbipyramidal("[As@@H](C)(N)(O)S", 0, 4, 1, 2, 3, Stereo::implRef());
  test_trigonalbipyramidal("[As@TB2H](C)(N)(O)S", 0, 4, 1, 2, 3, Stereo::implRef());

  test_trigonalbipyramidal("C[As@H2](N)O", 1, 0, Stereo::implRef(), Stereo::implRef(), 2, 3);
  test_trigonalbipyramidal("C[As@TB1H2](N)O", 1, 0, Stereo::implRef(), Stereo::implRef(), 2, 3);
  test_trigonalbipyramidal("C[As@@H2](N)O", 1, 3, Stereo::implRef(), Stereo::implRef(), 2, 0);
  test_trigonalbipyramidal("C[As@TB2H2](N)O", 1, 3, Stereo::implRef(), Stereo::implRef(), 2, 0);

  test_trigonalbipyramidal("[As@H2](C)(N)O", 0, Stereo::implRef(), Stereo::implRef(), 1, 2, 3);
  test_trigonalbipyramidal("[As@TB1H2](C)(N)O", 0, Stereo::implRef(), Stereo::implRef(), 1, 2, 3);
  test_trigonalbipyramidal("[As@@H2](C)(N)O", 0, 3, Stereo::implRef(), 1, 2, Stereo::implRef());
  test_trigonalbipyramidal("[As@TB2H2](C)(N)O", 0, 3, Stereo::implRef(), 1, 2, Stereo::implRef());

  // TB3/TB4: a-d
  test_trigonalbipyramidal("C[As@TB3](N)(O)(S)P", 1, 0, 2, 3, 5, 4);
  test_trigonalbipyramidal("C[As@TB4](N)(O)(S)P", 1, 4, 2, 3, 5, 0);
  // TB5/TB6: a-c
  test_trigonalbipyramidal("C[As@TB5](N)(O)(S)P", 1, 0, 2, 4, 5, 3);
  test_trigonalbipyramidal("C[As@TB6](N)(O)(S)P", 1, 3, 2, 4, 5, 0);
  // TB7/TB8: a-b
  test_trigonalbipyramidal("C[As@TB7](N)(O)(S)P", 1, 0, 3, 4, 5, 2);
  test_trigonalbipyramidal("C[As@TB8](N)(O)(S)P", 1, 2, 3, 4, 5, 0);
  // TB9/TB11: b-e
  test_trigonalbipyramidal("C[As@TB9](N)(O)(S)P",  1, 2, 0, 3, 4, 5);
  test_trigonalbipyramidal("C[As@TB11](N)(O)(S)P", 1, 5, 0, 3, 4, 2);
  // TB10/TB12: b-d
  test_trigonalbipyramidal("C[As@TB10](N)(O)(S)P", 1, 2, 0, 3, 5, 4);
  test_trigonalbipyramidal("C[As@TB12](N)(O)(S)P", 1, 4, 0, 3, 5, 2);
  // TB13/TB14: b-c
  test_trigonalbipyramidal("C[As@TB13](N)(O)(S)P",  1, 2, 0, 4, 5, 3);
  test_trigonalbipyramidal("C[As@TB14](N)(O)(S)P",  1, 3, 0, 4, 5, 2);
  // TB15/TB20: c-e
  test_trigonalbipyramidal("C[As@TB15](N)(O)(S)P",  1, 3, 0, 2, 4, 5);
  test_trigonalbipyramidal("C[As@TB20](N)(O)(S)P",  1, 5, 0, 2, 4, 3);
  // TB16/TB19: c-d
  test_trigonalbipyramidal("C[As@TB16](N)(O)(S)P",  1, 3, 0, 2, 5, 4);
  test_trigonalbipyramidal("C[As@TB19](N)(O)(S)P",  1, 4, 0, 2, 5, 3);
  // TB17/TB18: d-e
  test_trigonalbipyramidal("C[As@TB17](N)(O)(S)P",  1, 4, 0, 2, 3, 5);
  test_trigonalbipyramidal("C[As@TB18](N)(O)(S)P",  1, 5, 0, 2, 3, 4);
}

void test_octahedral(const std::string &smiles, Stereo::Ref center,
    Stereo::Ref ref1, Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4,
    Stereo::Ref ref5, Stereo::Ref ref6)
{
  std::cout << "Test octahedral: " << smiles << std::endl;
  HeMol mol;
  Stereochemistry stereo;
  Smiles SMILES;

  ASSERT(SMILES.read(smiles, mol, stereo));
  COMPARE(stereo.numStereo(), 1);
  COMPARE(stereo.numOctahedral(), 1);
  ASSERT(stereo.isOctahedral(center));
  ASSERT(stereo.octahedral(center).isValid());
  std::cout << stereo.octahedral(center) << std::endl;
  COMPARE(stereo.octahedral(center).type(), Stereo::Octahedral);
  COMPARE(stereo.octahedral(center).center(), center);
  COMPARE(stereo.octahedral(center).ref(0), ref1);
  COMPARE(stereo.octahedral(center).ref(1), ref2);
  COMPARE(stereo.octahedral(center).ref(2), ref3);
  COMPARE(stereo.octahedral(center).ref(3), ref4);
  COMPARE(stereo.octahedral(center).ref(4), ref5);
  COMPARE(stereo.octahedral(center).ref(5), ref6);
}

void test_octahedral()
{
  // OH1/OH2
  test_octahedral("C[Co@](N)(O)(S)(P)F", 1, 0, 2, 3, 4, 5, 6);
  test_octahedral("C[Co@OH1](N)(O)(S)(P)F", 1, 0, 2, 3, 4, 5, 6);
  test_octahedral("C[Co@@](N)(O)(S)(P)F", 1, 0, 5, 4, 3, 2, 6);
  test_octahedral("C[Co@OH2](N)(O)(S)(P)F", 1, 0, 5, 4, 3, 2, 6);

  test_octahedral("[Co@](C)(N)(O)(S)(P)F", 0, 1, 2, 3, 4, 5, 6);
  test_octahedral("[Co@OH1](C)(N)(O)(S)(P)F", 0, 1, 2, 3, 4, 5, 6);
  test_octahedral("[Co@@](C)(N)(O)(S)(P)F", 0, 1, 5, 4, 3, 2, 6);
  test_octahedral("[Co@OH2](C)(N)(O)(S)(P)F", 0, 1, 5, 4, 3, 2, 6);

  test_octahedral("C[Co@H](N)(O)(S)P", 1, 0, Stereo::implRef(), 2, 3, 4, 5);
  test_octahedral("C[Co@OH1H](N)(O)(S)P", 1, 0, Stereo::implRef(), 2, 3, 4, 5);
  test_octahedral("C[Co@@H](N)(O)(S)P", 1, 0, 4, 3, 2, Stereo::implRef(), 5);
  test_octahedral("C[Co@OH2H](N)(O)(S)P", 1, 0, 4, 3, 2, Stereo::implRef(), 5);

  test_octahedral("[Co@H](C)(N)(O)(S)P", 0, Stereo::implRef(), 1, 2, 3, 4, 5);
  test_octahedral("[Co@OH1H](C)(N)(O)(S)P", 0, Stereo::implRef(), 1, 2, 3, 4, 5);
  test_octahedral("[Co@@H](C)(N)(O)(S)P", 0, Stereo::implRef(), 4, 3, 2, 1, 5);
  test_octahedral("[Co@OH2H](C)(N)(O)(S)P", 0, Stereo::implRef(), 4, 3, 2, 1, 5);

  test_octahedral("C[Co@H2](N)(O)S", 1, 0, Stereo::implRef(), Stereo::implRef(), 2, 3, 4);
  test_octahedral("C[Co@OH1H2](N)(O)S", 1, 0, Stereo::implRef(), Stereo::implRef(), 2, 3, 4);
  test_octahedral("C[Co@@H2](N)(O)S", 1, 0, 3, 2, Stereo::implRef(), Stereo::implRef(), 4);
  test_octahedral("C[Co@OH2H2](N)(O)S", 1, 0, 3, 2, Stereo::implRef(), Stereo::implRef(), 4);

  test_octahedral("[Co@H2](C)(N)(O)S", 0, Stereo::implRef(), Stereo::implRef(), 1, 2, 3, 4);
  test_octahedral("[Co@OH1H2](C)(N)(O)S", 0, Stereo::implRef(), Stereo::implRef(), 1, 2, 3, 4);
  test_octahedral("[Co@@H2](C)(N)(O)S", 0, Stereo::implRef(), 3, 2, 1, Stereo::implRef(), 4);
  test_octahedral("[Co@OH2H2](C)(N)(O)S", 0, Stereo::implRef(), 3, 2, 1, Stereo::implRef(), 4);

  // OH3/16: a-e U-shape
  test_octahedral("[Co@OH3](C)(N)(O)(S)(P)F",  0, 1, 2, 3, 4, 6, 5);
  test_octahedral("[Co@OH16](C)(N)(O)(S)(P)F", 0, 5, 2, 3, 4, 6, 1);
  // OH6/OH18: a-d U-shape
  test_octahedral("[Co@OH6](C)(N)(O)(S)(P)F",  0, 1, 2, 3, 5, 6, 4);
  test_octahedral("[Co@OH18](C)(N)(O)(S)(P)F", 0, 4, 2, 3, 5, 6, 1);
  // OH19/OH24: a-c U-shape
  test_octahedral("[Co@OH19](C)(N)(O)(S)(P)F", 0, 1, 2, 4, 5, 6, 3);
  test_octahedral("[Co@OH24](C)(N)(O)(S)(P)F", 0, 3, 2, 4, 5, 6, 1);
  // OH25/OH30: a-b U-shape
  test_octahedral("[Co@OH25](C)(N)(O)(S)(P)F", 0, 1, 3, 4, 5, 6, 2);
  test_octahedral("[Co@OH30](C)(N)(O)(S)(P)F", 0, 2, 3, 4, 5, 6, 1);
  // OH4/OH14: a-f Z-shape
  test_octahedral("[Co@OH4](C)(N)(O)(S)(P)F",  0, 1, 2, 3, 5, 4, 6);
  test_octahedral("[Co@OH14](C)(N)(O)(S)(P)F", 0, 6, 2, 3, 5, 4, 1);
  // OH5/OH15: a-e Z-shape
  test_octahedral("[Co@OH5](C)(N)(O)(S)(P)F",  0, 1, 2, 3, 6, 4, 5);
  test_octahedral("[Co@OH15](C)(N)(O)(S)(P)F", 0, 5, 2, 3, 6, 4, 1);
  // OH7/OH17: a-d Z-shape
  test_octahedral("[Co@OH7](C)(N)(O)(S)(P)F",  0, 1, 2, 3, 6, 5, 4);
  test_octahedral("[Co@OH17](C)(N)(O)(S)(P)F", 0, 4, 2, 3, 6, 5, 1);
  // OH20/OH23: a-c Z-shape
  test_octahedral("[Co@OH20](C)(N)(O)(S)(P)F", 0, 1, 2, 4, 6, 5, 3);
  test_octahedral("[Co@OH23](C)(N)(O)(S)(P)F", 0, 3, 2, 4, 6, 5, 1);
  // OH26/OH29: a-b Z-shape
  test_octahedral("[Co@OH26](C)(N)(O)(S)(P)F", 0, 1, 3, 4, 6, 5, 2);
  test_octahedral("[Co@OH29](C)(N)(O)(S)(P)F", 0, 2, 3, 4, 6, 5, 1);
  // OH10/OH8: a-f 4-shape
  test_octahedral("[Co@OH10](C)(N)(O)(S)(P)F", 0, 1, 5, 3, 4, 2, 6);
  test_octahedral("[Co@OH8](C)(N)(O)(S)(P)F",  0, 6, 5, 3, 4, 2, 1);
  // OH11/OH9: a-e 4-shape
  test_octahedral("[Co@OH11](C)(N)(O)(S)(P)F", 0, 1, 6, 3, 4, 2, 5);
  test_octahedral("[Co@OH9](C)(N)(O)(S)(P)F",  0, 5, 6, 3, 4, 2, 1);
  // OH13/OH12: a-d 4-shape
  test_octahedral("[Co@OH13](C)(N)(O)(S)(P)F", 0, 1, 6, 3, 5, 2, 4);
  test_octahedral("[Co@OH12](C)(N)(O)(S)(P)F", 0, 4, 6, 3, 5, 2, 1);
  // OH22/OH21: a-c 4-shape
  test_octahedral("[Co@OH22](C)(N)(O)(S)(P)F", 0, 1, 6, 4, 5, 2, 3);
  test_octahedral("[Co@OH21](C)(N)(O)(S)(P)F", 0, 3, 6, 4, 5, 2, 1);
  // OH28/OH27: a-b 4-shape
  test_octahedral("[Co@OH28](C)(N)(O)(S)(P)F", 0, 1, 6, 4, 5, 3, 2);
  test_octahedral("[Co@OH27](C)(N)(O)(S)(P)F", 0, 2, 6, 4, 5, 3, 1);
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

  // cis/trans stereochemistry
  test_roundtrip("F\\C=C\\F", "F/C=C/F");
  test_roundtrip("F/C=C/F", "F\\C=C\\F");
  test_roundtrip("F/C=C/?I", "FC=CI");
  test_roundtrip("F/C=C\\?I", "FC=CI");
  test_roundtrip("F/C=C/C=C/F", "F\\C=C\\C=C\\F");

  // square-planar stereochemistry
  test_roundtrip("C[Pt@SP1](N)(O)S");
  test_roundtrip("C[Pt@SP2](N)(O)S");
  test_roundtrip("C[Pt@SP3](N)(O)S");
  test_roundtrip("[Pt@SP1](C)(N)(O)S");
  test_roundtrip("[Pt@SP2](C)(N)(O)S");
  test_roundtrip("[Pt@SP3](C)(N)(O)S");
  test_roundtrip("C[Pt@SP1H](N)O");
  test_roundtrip("C[Pt@SP2H](N)O");
  test_roundtrip("C[Pt@SP3H](N)O");
  test_roundtrip("[Pt@SP1H](C)(N)O");
  test_roundtrip("[Pt@SP2H](C)(N)O");
  test_roundtrip("[Pt@SP3H](C)(N)O");
  test_roundtrip("C[Pt@SP1H2]N");
  test_roundtrip("C[Pt@SP2H2]N", "C[Pt@SP1H2]N");
  test_roundtrip("C[Pt@SP3H2]N");
  test_roundtrip("[Pt@SP1H2](C)N");
  test_roundtrip("[Pt@SP2H2](C)N");
  test_roundtrip("[Pt@SP3H2](C)N", "[Pt@SP1H2](C)N");

  // trigonalbipyramidal stereochemistry
  test_roundtrip("C[As@](N)(O)(S)P", "C[As@TB1](N)(O)(S)P");
  test_roundtrip("C[As@@](N)(O)(S)P", "C[As@TB2](N)(O)(S)P");
  test_roundtrip("C[As@TB1](N)(O)(S)P");
  test_roundtrip("C[As@TB2](N)(O)(S)P");
  test_roundtrip("[As@TB1](C)(N)(O)(S)P");
  test_roundtrip("[As@TB2](C)(N)(O)(S)P");
  test_roundtrip("[As@TB1H](C)(N)(O)S");
  test_roundtrip("[As@TB2H](C)(N)(O)S");
  test_roundtrip("[As@TB1H2](C)(N)O");
  test_roundtrip("[As@TB2H2](C)(N)O");
  test_roundtrip("C[As@TB3](N)(O)(S)P");
  test_roundtrip("C[As@TB4](N)(O)(S)P");
  test_roundtrip("C[As@TB5](N)(O)(S)P");
  test_roundtrip("C[As@TB6](N)(O)(S)P");
  test_roundtrip("C[As@TB7](N)(O)(S)P");
  test_roundtrip("C[As@TB8](N)(O)(S)P");
  test_roundtrip("C[As@TB9](N)(O)(S)P");
  test_roundtrip("C[As@TB10](N)(O)(S)P");
  test_roundtrip("C[As@TB11](N)(O)(S)P");
  test_roundtrip("C[As@TB12](N)(O)(S)P");
  test_roundtrip("C[As@TB13](N)(O)(S)P");
  test_roundtrip("C[As@TB14](N)(O)(S)P");
  test_roundtrip("C[As@TB15](N)(O)(S)P");
  test_roundtrip("C[As@TB16](N)(O)(S)P");
  test_roundtrip("C[As@TB17](N)(O)(S)P");
  test_roundtrip("C[As@TB18](N)(O)(S)P");
  test_roundtrip("C[As@TB19](N)(O)(S)P");
  test_roundtrip("C[As@TB20](N)(O)(S)P");

  // octahedral stereochemistry
  test_roundtrip("C[Co@](N)(O)(S)(P)F", "C[Co@OH1](N)(O)(S)(P)F");
  test_roundtrip("C[Co@@](N)(O)(S)(P)F", "C[Co@OH2](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH1](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH2](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH1H](N)(O)(S)P");
  test_roundtrip("C[Co@OH2H](N)(O)(S)P");
  test_roundtrip("C[Co@OH1H2](N)(O)S");
  test_roundtrip("C[Co@OH2H2](N)(O)S");
  test_roundtrip("C[Co@OH3](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH4](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH5](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH6](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH7](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH8](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH9](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH10](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH11](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH12](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH13](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH14](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH15](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH16](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH17](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH18](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH19](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH20](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH21](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH22](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH23](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH24](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH25](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH26](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH27](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH28](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH29](N)(O)(S)(P)F");
  test_roundtrip("C[Co@OH30](N)(O)(S)(P)F");

  // unspecified stereochemistry
  test_roundtrip("C[C@?](N)(O)S", "CC(N)(O)S");
  test_roundtrip("C[C@TH?](N)(O)S", "CC(N)(O)S");
  test_roundtrip("CC=[C@AL?]=C(N)O", "CC=C=C(N)O");
  test_roundtrip("C[Pt@SP?](N)(O)S", "C[Pt](N)(O)S");
  test_roundtrip("C[As@TB?](N)(O)(S)P", "C[As](N)(O)(S)P");
  test_roundtrip("C[Co@OH?](N)(O)(S)(P)F", "C[Co](N)(O)(S)(P)F");
  test_roundtrip("F/?C=C/F", "FC=CF");
  test_roundtrip("F/C=C/?F", "FC=CF");
  test_roundtrip("F/?C=C/?F", "FC=CF");

  // invalid stereochemistry
  test_roundtrip("C[C@](N)O", "C[C](N)O"); // TH 3 nbrs
  test_roundtrip("C[C@TH1](N)(O)(S)P", "CC(N)(O)(S)P"); // TH 5 nbrs
  test_roundtrip("C[C@SP1](N)O", "C[C](N)O"); // SP 3 nbrs
  test_roundtrip("C[C@SP1](N)(O)(S)P", "CC(N)(O)(S)P"); // SP 5 nbrs
  test_roundtrip("C[C@TB1](N)(O)S", "CC(N)(O)S"); // TB 4 nbrs
  test_roundtrip("C[C@TB1](N)(O)(S)(P)F", "CC(N)(O)(S)(P)F"); // TB 6 nbrs
  test_roundtrip("C[C@OH1](N)(O)(S)P", "CC(N)(O)(S)P"); // OH 5 nbrs
  test_roundtrip("C[C@OH1](N)(O)(S)(P)(F)I", "CC(N)(O)(S)(P)(F)I"); // OH 7 nbrs


  //test_roundtrip("");
}

int main()
{
  test_parse_smiles();
  test_tetrahedral();
  test_allene();
  test_cistrans();
  test_squareplanar();
  test_trigonalbipyramidal();
  test_octahedral();
  test_roundtrip();
}
