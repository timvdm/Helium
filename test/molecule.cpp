#include <Helium/hemol.h>

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

int main()
{
  test_atom_predicates();
  test_molecule_has_atom();
}
