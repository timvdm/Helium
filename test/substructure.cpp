#include "../src/substructure.h"
#include "../src/smiles.h"
#include "../src/concepts.h"

#include "test.h"

using namespace Helium;

void test_substructure1()
{
  HeMol mol = hemol_from_smiles("c1ccccc1c2ccccc2");

  std::vector<bool> atoms1(num_atoms(mol));
  std::vector<bool> bonds1(num_bonds(mol));
  std::vector<bool> atoms2(num_atoms(mol));
  std::vector<bool> bonds2(num_bonds(mol));
  for (int i = 0; i < 6; ++i) {
    atoms1[i] = true;
    bonds1[i] = true;
    atoms2[i + 6] = true;
    bonds2[i + 7] = true;
  }


  Substructure<HeMol> substruct1(mol, atoms1, bonds1);
  Substructure<HeMol> substruct2(mol, atoms2, bonds2);

  COMPARE(6, substruct2.oldAtomIndex(get_atom(substruct2, 0)));
  COMPARE(7, substruct2.oldAtomIndex(get_atom(substruct2, 1)));
  COMPARE(8, substruct2.oldAtomIndex(get_atom(substruct2, 2)));
  COMPARE(9, substruct2.oldAtomIndex(get_atom(substruct2, 3)));
  COMPARE(10, substruct2.oldAtomIndex(get_atom(substruct2, 4)));
  COMPARE(11, substruct2.oldAtomIndex(get_atom(substruct2, 5)));

  COMPARE(7, substruct2.oldBondIndex(get_bond(substruct2, 0)));
  COMPARE(8, substruct2.oldBondIndex(get_bond(substruct2, 1)));
  COMPARE(9, substruct2.oldBondIndex(get_bond(substruct2, 2)));
  COMPARE(10, substruct2.oldBondIndex(get_bond(substruct2, 3)));
  COMPARE(11, substruct2.oldBondIndex(get_bond(substruct2, 4)));
  COMPARE(12, substruct2.oldBondIndex(get_bond(substruct2, 5)));
}


int main()
{
  HeMol mol;
  SMILES.read("CNOSP", mol);

  std::vector<bool> atoms(5), bonds(4);
  atoms[1] = true;
  atoms[2] = true;
  atoms[3] = true;
  bonds[1] = true;
  bonds[2] = true;

  Substructure<HeMol> substruct(mol, atoms, bonds);

  ////////////////////////////////////////////////
  //
  // Substructure
  //
  ////////////////////////////////////////////////

  // test num_atoms
  COMPARE(3, num_atoms(substruct));

  // test get_atoms
  unsigned int numAtoms = 0;
  for (auto &atom : get_atoms(substruct))
    ++numAtoms;
  COMPARE(3, numAtoms);

  // test get_atom
  COMPARE(get_atom(mol, 1), get_atom(substruct, 0));
  COMPARE(get_atom(mol, 2), get_atom(substruct, 1));
  COMPARE(get_atom(mol, 3), get_atom(substruct, 2));

  // test num_bonds
  COMPARE(2, num_bonds(substruct));

  // test get_bonds
  unsigned int numBonds = 0;
  for (auto &bond : get_bonds(substruct))
    ++numBonds;
  COMPARE(2, numBonds);

  // test get_bond
  COMPARE(get_bond(mol, 1), get_bond(substruct, 0));
  COMPARE(get_bond(mol, 2), get_bond(substruct, 1));

  COMPARE(get_bond(mol, 1), get_bond(substruct, get_atom(substruct, 0), get_atom(substruct, 1)));
  COMPARE(get_bond(mol, 2), get_bond(substruct, get_atom(substruct, 1), get_atom(substruct, 2)));

  ////////////////////////////////////////////////
  //
  // Atom
  //
  ////////////////////////////////////////////////

  // test get_index
  COMPARE(0, get_index(substruct, get_atom(substruct, 0)));
  COMPARE(1, get_index(substruct, get_atom(substruct, 1)));
  COMPARE(2, get_index(substruct, get_atom(substruct, 2)));

  // test get_bonds
  numBonds = 0;
  for (auto &bond : get_bonds(substruct, get_atom(substruct, 0)))
    ++numBonds;
  COMPARE(1, numBonds);

  numBonds = 0;
  for (auto &bond : get_bonds(substruct, get_atom(substruct, 1)))
    ++numBonds;
  COMPARE(2, numBonds);

  numBonds = 0;
  for (auto &bond : get_bonds(substruct, get_atom(substruct, 2)))
    ++numBonds;
  COMPARE(1, numBonds);

  // test get_nbrs
  unsigned int numNbrs = 0;
  for (auto &nbr : get_nbrs(substruct, get_atom(substruct, 0)))
    ++numNbrs;
  COMPARE(1, numNbrs);

  numNbrs = 0;
  for (auto &nbr : get_nbrs(substruct, get_atom(substruct, 1)))
    ++numNbrs;
  COMPARE(2, numNbrs);

  numNbrs = 0;
  for (auto &nbr : get_nbrs(substruct, get_atom(substruct, 2)))
    ++numNbrs;
  COMPARE(1, numNbrs);

  // test get_degree
  COMPARE(1, get_degree(substruct, get_atom(substruct, 0)));
  COMPARE(2, get_degree(substruct, get_atom(substruct, 1)));
  COMPARE(1, get_degree(substruct, get_atom(substruct, 2)));

  ////////////////////////////////////////////////
  //
  // Atom
  //
  ////////////////////////////////////////////////

  // test get_index
  COMPARE(0, get_index(substruct, get_bond(substruct, 0)));
  COMPARE(1, get_index(substruct, get_bond(substruct, 1)));


  test_substructure1();

  //check_molecule_concept(Substructure<HeMol>(HeMol(), std::vector<bool>(), std::vector<bool>()));
}
