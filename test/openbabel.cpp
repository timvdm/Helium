#include <Helium/toolkits/openbabel.h>

#include "test.h"

using namespace Helium;
using namespace OpenBabel;

template<typename EditableMoleculeType>
void test_atom_functions(EditableMoleculeType &mol)
{
  typedef typename molecule_traits<EditableMoleculeType>::atom_type atom_type;
  typedef typename molecule_traits<EditableMoleculeType>::bond_type bond_type;

  COMPARE(0, num_atoms(mol));
  COMPARE(0, num_bonds(mol));

  atom_type atom1 = add_atom(mol);

  COMPARE(1, num_atoms(mol));

  // atom index
  COMPARE(0, get_index(mol, atom1));
  // atom aromaticity
  COMPARE(false, is_aromatic(mol, atom1));
  set_aromatic(mol, atom1, true);
  COMPARE(true, is_aromatic(mol, atom1));
  // atom element
  set_element(mol, atom1, 7);
  COMPARE(7, get_element(mol, atom1));
  set_element(mol, atom1, 9);
  COMPARE(9, get_element(mol, atom1));
  // atom mass
  set_mass(mol, atom1, 13);
  COMPARE(13, get_mass(mol, atom1));
  set_mass(mol, atom1, 2);
  COMPARE(2, get_mass(mol, atom1));
  // atom number of implicit hydrogens
  set_hydrogens(mol, atom1, 4);
  COMPARE(4, get_hydrogens(mol, atom1));
  set_hydrogens(mol, atom1, 3);
  COMPARE(3, get_hydrogens(mol, atom1));
  // atom charge
  set_charge(mol, atom1, -1);
  COMPARE(-1, get_charge(mol, atom1));
  set_charge(mol, atom1, 2);
  COMPARE(2, get_charge(mol, atom1));
  // atom degree
  COMPARE(0, get_degree(mol, atom1));

  //
  // add more atoms
  //
  atom_type atom2 = add_atom(mol);
  COMPARE(2, num_atoms(mol));
  COMPARE(1, get_index(mol, atom2));

  atom_type atom3 = add_atom(mol);
  COMPARE(3, num_atoms(mol));
  COMPARE(2, get_index(mol, atom3));

  atom_type atom4 = add_atom(mol);
  COMPARE(4, num_atoms(mol));
  COMPARE(3, get_index(mol, atom4));

  // add bond
  COMPARE(0, num_bonds(mol));
  bond_type bond1 = add_bond(mol, atom1, atom2);
  COMPARE(1, num_bonds(mol));
  COMPARE(0, get_index(mol, bond1));
  COMPARE(atom1, get_source(mol, bond1));
  COMPARE(atom2, get_target(mol, bond1));
  COMPARE(atom2, get_other(mol, bond1, atom1));
  COMPARE(atom1, get_other(mol, bond1, atom2));
  COMPARE(false, is_aromatic(mol, bond1));


  // bond aromaticity
  set_aromatic(mol, bond1, true);
  COMPARE(true, is_aromatic(mol, bond1));
  set_aromatic(mol, bond1, false);
  COMPARE(false, is_aromatic(mol, bond1));
  // bond order
  set_order(mol, bond1, 2);
  COMPARE(2, get_order(mol, bond1));
  set_order(mol, bond1, 1);
  COMPARE(1, get_order(mol, bond1));

  // atom degree
  COMPARE(1, get_degree(mol, atom1));
  COMPARE(1, get_degree(mol, atom2));
  COMPARE(0, get_degree(mol, atom3));
  COMPARE(0, get_degree(mol, atom4));

  bond_type bond2 = add_bond(mol, atom2, atom3);
  COMPARE(1, get_index(mol, bond2));

  // atom degree
  COMPARE(1, get_degree(mol, atom1));
  COMPARE(2, get_degree(mol, atom2));
  COMPARE(1, get_degree(mol, atom3));
  COMPARE(0, get_degree(mol, atom4));

  bond_type bond3 = add_bond(mol, atom2, atom4);
  COMPARE(2, get_index(mol, bond3));

  // atom degree
  COMPARE(1, get_degree(mol, atom1));
  COMPARE(3, get_degree(mol, atom2));
  COMPARE(1, get_degree(mol, atom3));
  COMPARE(1, get_degree(mol, atom4));

  //
  // iterators
  //
  OpenBabel::OBAtom* atoms[4] = { atom1, atom2, atom3, atom4 };

  int i = 0;
  FOREACH_ATOM_T (atom, mol, EditableMoleculeType) {
    COMPARE(atoms[i], *atom);
    ++i;
  }
  COMPARE(i, 4);

  OpenBabel::OBAtom* sources[3] = { atom1, atom2, atom2 };
  OpenBabel::OBAtom* targets[3] = { atom2, atom3, atom4 };

  i = 0;
  FOREACH_BOND_T (bond, mol, EditableMoleculeType) {
    COMPARE(sources[i], get_source(mol, *bond));
    COMPARE(targets[i], get_target(mol, *bond));
    ++i;
  }
  COMPARE(i, 3);

  OpenBabel::OBAtom* nbrs[3] = { atom1, atom3, atom4 };

  i = 0;
  FOREACH_NBR_T (nbr, atom2, mol, EditableMoleculeType) {
    COMPARE(nbrs[i], *nbr);
    ++i;
  }
  COMPARE(i, 3);

  OpenBabel::OBBond* bonds[3] = { bond1, bond2, bond3 };

  i = 0;
  FOREACH_INCIDENT_T (bond, atom2, mol, EditableMoleculeType) {
    COMPARE(bonds[i], *bond);
    ++i;
  }
  COMPARE(i, 3);

}

int main()
{
  OBMol mol;
  test_atom_functions(mol);
}
