#include <iterator>
#include <Helium/hemol.h>
#include <Helium/smartmol.h>
#include <Helium/concepts.h>


#include "test.h"

using namespace Helium;

void test_element_predicate()
{
  // ElementPredicate
  HeMol mol = hemol_from_smiles("C");
  ASSERT( molecule_has_atom(mol, element_eq_predicate(mol, 6)));
  ASSERT(!molecule_has_atom(mol, element_eq_predicate(mol, 7)));
  ASSERT( molecule_has_atom(mol, element_lt_predicate(mol, 7)));
  ASSERT(!molecule_has_atom(mol, element_lt_predicate(mol, 6)));
  ASSERT( molecule_has_atom(mol, element_gt_predicate(mol, 5)));
  ASSERT(!molecule_has_atom(mol, element_gt_predicate(mol, 6)));
  ASSERT( molecule_has_atom(mol, element_leq_predicate(mol, 6)));
  ASSERT( molecule_has_atom(mol, element_leq_predicate(mol, 7)));
  ASSERT(!molecule_has_atom(mol, element_leq_predicate(mol, 5)));
  ASSERT( molecule_has_atom(mol, element_geq_predicate(mol, 5)));
  ASSERT( molecule_has_atom(mol, element_geq_predicate(mol, 6)));
  ASSERT(!molecule_has_atom(mol, element_geq_predicate(mol, 7)));
}

void test_mass_predicate()
{
  // MassPredicate
  HeMol mol = hemol_from_smiles("[13C]");
  ASSERT( molecule_has_atom(mol, mass_eq_predicate(mol, 13)));
  ASSERT(!molecule_has_atom(mol, mass_eq_predicate(mol, 12)));
  ASSERT( molecule_has_atom(mol, mass_lt_predicate(mol, 14)));
  ASSERT(!molecule_has_atom(mol, mass_lt_predicate(mol, 13)));
  ASSERT( molecule_has_atom(mol, mass_gt_predicate(mol, 12)));
  ASSERT(!molecule_has_atom(mol, mass_gt_predicate(mol, 13)));
  ASSERT( molecule_has_atom(mol, mass_leq_predicate(mol, 13)));
  ASSERT(!molecule_has_atom(mol, mass_leq_predicate(mol, 12)));
  ASSERT( molecule_has_atom(mol, mass_geq_predicate(mol, 13)));
  ASSERT(!molecule_has_atom(mol, mass_geq_predicate(mol, 14)));
}

void test_charge_predicate()
{
  // ChargePredicate
  HeMol mol = hemol_from_smiles("C[NH3+]");
  ASSERT( molecule_has_atom(mol, charge_eq_predicate(mol, 1)));
  ASSERT( molecule_has_atom(mol, charge_eq_predicate(mol, 0)));
  ASSERT(!molecule_has_atom(mol, charge_eq_predicate(mol, -1)));
  ASSERT( molecule_has_atom(mol, charge_lt_predicate(mol, 2)));
  ASSERT( molecule_has_atom(mol, charge_lt_predicate(mol, 1)));
  ASSERT(!molecule_has_atom(mol, charge_lt_predicate(mol, 0)));
  ASSERT( molecule_has_atom(mol, charge_gt_predicate(mol, 0)));
  ASSERT(!molecule_has_atom(mol, charge_gt_predicate(mol, 1)));
  ASSERT( molecule_has_atom(mol, charge_leq_predicate(mol, 1)));
  ASSERT( molecule_has_atom(mol, charge_leq_predicate(mol, 0)));
  ASSERT(!molecule_has_atom(mol, charge_leq_predicate(mol, -1)));
  ASSERT( molecule_has_atom(mol, charge_geq_predicate(mol, 1)));
  ASSERT(!molecule_has_atom(mol, charge_geq_predicate(mol, 2)));
}

void test_degree_predicate()
{
  // DegreePredicate
  HeMol mol = hemol_from_smiles("CC(O)C");
  ASSERT( molecule_has_atom(mol, degree_eq_predicate(mol, 1)));
  ASSERT(!molecule_has_atom(mol, degree_eq_predicate(mol, 2)));
  ASSERT( molecule_has_atom(mol, degree_eq_predicate(mol, 3)));
  ASSERT( molecule_has_atom(mol, degree_lt_predicate(mol, 2)));
  ASSERT(!molecule_has_atom(mol, degree_lt_predicate(mol, 1)));
  ASSERT( molecule_has_atom(mol, degree_gt_predicate(mol, 2)));
  ASSERT(!molecule_has_atom(mol, degree_gt_predicate(mol, 3)));
  ASSERT( molecule_has_atom(mol, degree_leq_predicate(mol, 1)));
  ASSERT(!molecule_has_atom(mol, degree_leq_predicate(mol, 0)));
  ASSERT( molecule_has_atom(mol, degree_geq_predicate(mol, 3)));
  ASSERT(!molecule_has_atom(mol, degree_geq_predicate(mol, 4)));
}

void test_aromatic_atom_predicate()
{
  // AromaticAtomPredicate
  HeMol mol1 = hemol_from_smiles("c1ccccc1");
  ASSERT( molecule_has_atom(mol1, aromatic_atom_predicate(mol1)));
  ASSERT(!molecule_has_atom(mol1, not_aromatic_atom_predicate(mol1)));

  HeMol mol2 = hemol_from_smiles("C1CCCCC1");
  ASSERT(!molecule_has_atom(mol2, aromatic_atom_predicate(mol2)));
  ASSERT( molecule_has_atom(mol2, not_aromatic_atom_predicate(mol2)));
}

void test_atom_and_predicate()
{
  HeMol mol = hemol_from_smiles("[13C]");

  ASSERT( atom_and_predicates(mol, element_eq_predicate(mol, 6), mass_eq_predicate(mol, 13))(mol, get_atom(mol, 0)));
  ASSERT(!atom_and_predicates(mol, element_eq_predicate(mol, 7), mass_eq_predicate(mol, 13))(mol, get_atom(mol, 0)));
  ASSERT(!atom_and_predicates(mol, element_eq_predicate(mol, 6), mass_eq_predicate(mol, 14))(mol, get_atom(mol, 0)));
}

void test_atom_or_predicate()
{
  HeMol mol = hemol_from_smiles("[13C]");

  ASSERT( atom_or_predicates(mol, element_eq_predicate(mol, 6), mass_eq_predicate(mol, 13))(mol, get_atom(mol, 0)));
  ASSERT( atom_or_predicates(mol, element_eq_predicate(mol, 7), mass_eq_predicate(mol, 13))(mol, get_atom(mol, 0)));
  ASSERT( atom_or_predicates(mol, element_eq_predicate(mol, 6), mass_eq_predicate(mol, 14))(mol, get_atom(mol, 0)));
  ASSERT(!atom_or_predicates(mol, element_eq_predicate(mol, 7), mass_eq_predicate(mol, 14))(mol, get_atom(mol, 0)));
}

void test_atom_predicates()
{
  test_element_predicate();
  test_mass_predicate();
  test_charge_predicate();
  test_degree_predicate();
  test_aromatic_atom_predicate();
  test_atom_and_predicate();
  test_atom_or_predicate();
}

void test_molecule_has_atom()
{
  HeMol methane = hemol_from_smiles("C");

  ASSERT(molecule_has_atom(methane, ElementPredicate<HeMol>(6)));
  ASSERT(!molecule_has_atom(methane, ElementPredicate<HeMol>(7)));
}

template<typename MoleculeType>
void test_editable_molecule()
{
  typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
  typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
  MoleculeType mol;

  // create 4 atoms
  atom_type atom1 = add_atom(mol);
  COMPARE(1, num_atoms(mol));
  atom_type atom2 = add_atom(mol);
  COMPARE(2, num_atoms(mol));
  atom_type atom3 = add_atom(mol);
  COMPARE(3, num_atoms(mol));
  atom_type atom4 = add_atom(mol);
  COMPARE(4, num_atoms(mol));

  // set some atom properties
  set_aromatic(mol, atom1, true);
  set_aromatic(mol, atom2, true);
  set_aromatic(mol, atom3, false);
  set_aromatic(mol, atom4, false);
  COMPARE(true, is_aromatic(mol, atom1));
  COMPARE(true, is_aromatic(mol, atom2));
  COMPARE(false, is_aromatic(mol, atom3));
  COMPARE(false, is_aromatic(mol, atom4));

  set_element(mol, atom1, 6);
  set_element(mol, atom2, 7);
  set_element(mol, atom3, 8);
  set_element(mol, atom4, 9);
  COMPARE(6, get_element(mol, atom1));
  COMPARE(7, get_element(mol, atom2));
  COMPARE(8, get_element(mol, atom3));
  COMPARE(9, get_element(mol, atom4));

  set_mass(mol, atom1, 60);
  set_mass(mol, atom2, 70);
  set_mass(mol, atom3, 80);
  set_mass(mol, atom4, 90);
  COMPARE(60, get_mass(mol, atom1));
  COMPARE(70, get_mass(mol, atom2));
  COMPARE(80, get_mass(mol, atom3));
  COMPARE(90, get_mass(mol, atom4));

  set_hydrogens(mol, atom1, 1);
  set_hydrogens(mol, atom2, 2);
  set_hydrogens(mol, atom3, 3);
  set_hydrogens(mol, atom4, 4);
  COMPARE(1, get_hydrogens(mol, atom1));
  COMPARE(2, get_hydrogens(mol, atom2));
  COMPARE(3, get_hydrogens(mol, atom3));
  COMPARE(4, get_hydrogens(mol, atom4));

  set_charge(mol, atom1, -2);
  set_charge(mol, atom2, -1);
  set_charge(mol, atom3, 0);
  set_charge(mol, atom4, 2);
  COMPARE(-2, get_charge(mol, atom1));
  COMPARE(-1, get_charge(mol, atom2));
  COMPARE(0, get_charge(mol, atom3));
  COMPARE(2, get_charge(mol, atom4));

  // add 3 bonds
  bond_type bond1 = add_bond(mol, atom1, atom2);
  COMPARE(1, num_bonds(mol));
  bond_type bond2 = add_bond(mol, atom2, atom3);
  COMPARE(2, num_bonds(mol));
  bond_type bond3 = add_bond(mol, atom3, atom4);
  COMPARE(3, num_bonds(mol));

  // set some bond properties
  set_aromatic(mol, bond1, false);
  set_aromatic(mol, bond2, true);
  set_aromatic(mol, bond3, false);
  COMPARE(false, is_aromatic(mol, bond1));
  COMPARE(true, is_aromatic(mol, bond2));
  COMPARE(false, is_aromatic(mol, bond3));

  set_order(mol, bond1, 1);
  set_order(mol, bond2, 2);
  set_order(mol, bond3, 3);
  COMPARE(1, get_order(mol, bond1));
  COMPARE(2, get_order(mol, bond2));
  COMPARE(3, get_order(mol, bond3));

  // remove a bond
  remove_bond(mol, bond3);
  COMPARE(4, num_atoms(mol));
  COMPARE(2, num_bonds(mol));
  COMPARE(1, get_order(mol, get_bond(mol, 0)));
  COMPARE(2, get_order(mol, get_bond(mol, 1)));

  // remove an atom
  remove_atom(mol, atom2);
  COMPARE(3, num_atoms(mol));
  COMPARE(0, num_bonds(mol));
  COMPARE(60, get_mass(mol, get_atom(mol, 0)));
  COMPARE(80, get_mass(mol, get_atom(mol, 1)));
  COMPARE(90, get_mass(mol, get_atom(mol, 2)));

  // clear the molecule
  clear_molecule(mol);
  COMPARE(0, num_atoms(mol));
  COMPARE(0, num_bonds(mol));
}

void test_substructure()
{
  HeMol mol = hemol_from_smiles("Oc1cc(CC)ccc1N");

  std::vector<bool> atoms(mol.numAtoms(), true);
  atoms[0] = false;
  atoms[4] = false;
  atoms[5] = false;
  atoms[9] = false;

  std::vector<bool> bonds(mol.numBonds(), true);
  bonds[0] = false;
  bonds[3] = false;
  bonds[4] = false;
  bonds[9] = false;

  HeMol sub;
  make_substructure(sub, mol, atoms, bonds);

  Smiles SMILES;
  COMPARE("c1ccccc1", SMILES.write(sub));
}

int main()
{
  test_substructure();

  test_atom_predicates();
  test_molecule_has_atom();

  test_editable_molecule<HeMol>();
  test_editable_molecule<SmartMol>();

  HeMol hemol;
  //check_editable_molecule_concept(hemol);

  /*
  try {
    SmartMol smartmol;
    check_editable_molecule_concept(smartmol);
  } catch (...) {}
  */
}
