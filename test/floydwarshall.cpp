#include <Helium/algorithms/floydwarshall.h>
#include <Helium/smiles.h>
#include <Helium/fileio/moleculefile.h>

#include "test.h"

using namespace Helium;

int main()
{
  typedef molecule_traits<HeMol>::atom_type atom_type;

  HeMol mol;
  atom_type a1 = mol.addAtom();
  atom_type a2 = mol.addAtom();
  atom_type a3 = mol.addAtom();
  atom_type a4 = mol.addAtom();
  atom_type a5 = mol.addAtom();
  atom_type a6 = mol.addAtom();
  atom_type a7 = mol.addAtom();

  mol.addBond(a1, a2);
  mol.addBond(a2, a3);
  mol.addBond(a3, a4);
  mol.addBond(a4, a5);
  mol.addBond(a5, a6);
  mol.addBond(a6, a7);
  mol.addBond(a7, a1);
  mol.addBond(a2, a6);

  DistanceMatrix D = floyd_warshall(mol);

  COMPARE(0, D(0, 0));
  COMPARE(1, D(0, 1));
  COMPARE(2, D(0, 2));
  COMPARE(3, D(0, 3));
  COMPARE(4, D(0, 4));
  COMPARE(2, D(0, 5));
  COMPARE(3, D(0, 6));

  std::cout << D << std::endl;
}
