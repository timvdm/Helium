#include <Helium/algorithms/smirks.h>
#include <Helium/hemol.h>

#include "test.h"

using namespace Helium;

void test_atom_classes()
{
  Smirks smirks;
  ASSERT(!smirks.init("CCCC"));
  ASSERT(smirks.init("CC>>CC"));
  ASSERT(smirks.init("[C:1]C>>[C:1]C"));
  ASSERT(!smirks.init("[C:1][C:2]>>[C:1]C"));
  ASSERT(!smirks.init("[C:1]C>>[C:1][C:2]"));
  ASSERT(!smirks.init("[C:1][C:1]>>[C:1][C:2]"));
  ASSERT(!smirks.init("[C:1][C:2]>>[C:1][C:3]"));
  ASSERT(smirks.init("[C:1][C:2]>>[C:1][O:2]"));

  ASSERT(!smirks.init("[C:1]>>[C,N:1]"));
  ASSERT(!smirks.init("[C:1]>>[!C:1]"));
}

void test_simple_atoms()
{
  // C -> O
  {
    HeMol mol = hemol_from_smiles("C");
    Smirks smirks;
    smirks.init("[C:1]>>[O:1]");
    smirks.apply(mol, RingSet<HeMol>(mol));
    COMPARE(8, get_element(mol, get_atom(mol, 0)));
  }
  // C -> O-
  {
    HeMol mol = hemol_from_smiles("C");
    Smirks smirks;
    smirks.init("[C:1]>>[O-:1]");
    smirks.apply(mol, RingSet<HeMol>(mol));
    COMPARE(8, get_element(mol, get_atom(mol, 0)));
    COMPARE(-1, get_charge(mol, get_atom(mol, 0)));
  }
  // C -> 13C
  {
    HeMol mol = hemol_from_smiles("C");
    Smirks smirks;
    smirks.init("[C:1]>>[13C:1]");
    smirks.apply(mol, RingSet<HeMol>(mol));
    COMPARE(6, get_element(mol, get_atom(mol, 0)));
    COMPARE(13, get_mass(mol, get_atom(mol, 0)));
  }
  // C -> N+
  {
    HeMol mol = hemol_from_smiles("C");
    Smirks smirks;
    smirks.init("[C:1]>>[#7+:1]");
    smirks.apply(mol, RingSet<HeMol>(mol));
    COMPARE(7, get_element(mol, get_atom(mol, 0)));
    COMPARE(1, get_charge(mol, get_atom(mol, 0)));
  }
  // Ch4 -> Ch2
  {
    HeMol mol = hemol_from_smiles("C");
    Smirks smirks;
    smirks.init("[C:1]>>[Ch2:1]");
    smirks.apply(mol, RingSet<HeMol>(mol));
    COMPARE(6, get_element(mol, get_atom(mol, 0)));
    COMPARE(2, num_hydrogens(mol, get_atom(mol, 0)));
  }
  // CH4 -> CH3
  {
    HeMol mol = hemol_from_smiles("C");
    Smirks smirks;
    ASSERT(smirks.init("[C:1]>>[CH3:1]"));
    smirks.apply(mol, RingSet<HeMol>(mol));
    COMPARE(6, get_element(mol, get_atom(mol, 0)));
    COMPARE(3, num_hydrogens(mol, get_atom(mol, 0)));

    mol = hemol_from_smiles("C[H]");
    smirks.apply(mol, RingSet<HeMol>(mol));
    COMPARE(6, get_element(mol, get_atom(mol, 0)));
    COMPARE(2, num_hydrogens(mol, get_atom(mol, 0)));
  }
}

int main()
{
  Smirks smirks;
  smirks.init("[C:1][C:2]>>[C:1][O:2]");

  test_atom_classes();

  test_simple_atoms();
}
