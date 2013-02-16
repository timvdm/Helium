#include "../src/substructure.h"
#include "../src/smiles.h"

#include "test.h"

using namespace Helium;

int main()
{
  HeMol mol;
  parse_smiles("CNOSP", mol);

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
  molecule_traits<Substructure<HeMol> >::atom_iter atom, end_atoms;
  tie(atom, end_atoms) = get_atoms(substruct);
  for (; atom != end_atoms; ++atom)
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
  molecule_traits<Substructure<HeMol> >::bond_iter bond, end_bonds;
  tie(bond, end_bonds) = get_bonds(substruct);
  for (; bond != end_bonds; ++bond)
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
  molecule_traits<Substructure<HeMol> >::incident_iter atom_bond, end_atom_bonds;
  tie(atom_bond, end_atom_bonds) = get_bonds(substruct, get_atom(substruct, 0));
  for (; atom_bond != end_atom_bonds; ++atom_bond)
    ++numBonds;
  COMPARE(1, numBonds);

  numBonds = 0;
  tie(atom_bond, end_atom_bonds) = get_bonds(substruct, get_atom(substruct, 1));
  for (; atom_bond != end_atom_bonds; ++atom_bond)
    ++numBonds;
  COMPARE(2, numBonds);

  numBonds = 0;
  tie(atom_bond, end_atom_bonds) = get_bonds(substruct, get_atom(substruct, 2));
  for (; atom_bond != end_atom_bonds; ++atom_bond)
    ++numBonds;
  COMPARE(1, numBonds);

  // test get_nbrs
  unsigned int numNbrs = 0;
  molecule_traits<Substructure<HeMol> >::nbr_iter nbr, end_nbrs;
  tie(nbr, end_nbrs) = get_nbrs(substruct, get_atom(substruct, 0));
  for (; nbr != end_nbrs; ++nbr)
    ++numNbrs;
  COMPARE(1, numNbrs);

  numNbrs = 0;
  tie(nbr, end_nbrs) = get_nbrs(substruct, get_atom(substruct, 1));
  for (; nbr != end_nbrs; ++nbr)
    ++numNbrs;
  COMPARE(2, numNbrs);

  numNbrs = 0;
  tie(nbr, end_nbrs) = get_nbrs(substruct, get_atom(substruct, 2));
  for (; nbr != end_nbrs; ++nbr)
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

}
