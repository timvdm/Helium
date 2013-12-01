#include <Helium/algorithms/smarts.h>
#include <Helium/algorithms/cycles.h>
#include <Helium/hemol.h>

#include "test.h"

using namespace Helium;

void print_smarts(const std::string &smarts)
{
  Smiley::PrintCallback callback;
  Smiley::Parser<Smiley::PrintCallback> parser(callback, Smiley::Parser<Smiley::PrintCallback>::SmartsMode);

  try {
    parser.parse(smarts);
  } catch (Smiley::Exception &e) {
    std::ostringstream errorStream;
    if (e.type() == Smiley::Exception::SyntaxError)
      errorStream << "Syntax";
    else
      errorStream << "Semantics";
    errorStream << "Error: " << e.what() << "." << std::endl;
    errorStream << smarts << std::endl;
    for (std::size_t i = 0; i < e.pos(); ++i)
      errorStream << " ";
    for (std::size_t i = 0; i < e.length(); ++i)
      errorStream << "^";
    errorStream << std::endl;

    std::cerr << errorStream.str();
  }
}

void test_parser_error()
{
  Smarts smarts;
  ASSERT(smarts.init("C"));
  ASSERT(smarts.error().type() == Smiley::Exception::NoError);

  ASSERT(!smarts.init("gsfsf"));
  std::cout << smarts.error().what() << std::endl;
}

void test_organic_subset()
{
  // AE_AliphaticElement
  Smarts smarts1;
  smarts1.init("C");
  COMPARE(Smiley::AE_AliphaticElement, smarts1.trees().atom(0)->type);
  COMPARE(6, smarts1.trees().atom(0)->value);
  Smarts smarts3;
  smarts3.init("O");
  COMPARE(Smiley::AE_AliphaticElement, smarts3.trees().atom(0)->type);
  COMPARE(8, smarts3.trees().atom(0)->value);

  // AE_AromaticElement
  Smarts smarts2;
  smarts2.init("c");
  COMPARE(Smiley::AE_AromaticElement, smarts2.trees().atom(0)->type);
  COMPARE(6, smarts2.trees().atom(0)->value);
  Smarts smarts4;
  smarts4.init("o");
  COMPARE(Smiley::AE_AromaticElement, smarts4.trees().atom(0)->type);
  COMPARE(8, smarts4.trees().atom(0)->value);

  // AE_Aliphatic
  Smarts smarts5;
  smarts5.init("A");
  COMPARE(Smiley::AE_Aliphatic, smarts5.trees().atom(0)->type);

  // AE_Aromatic
  Smarts smarts6;
  smarts6.init("a");
  COMPARE(Smiley::AE_Aromatic, smarts6.trees().atom(0)->type);
}

void test_cyclic()
{
  // AE_Cyclic
  Smarts smarts1;
  smarts1.init("[R]");
  COMPARE(Smiley::AE_Cyclic, smarts1.trees().atom(0)->type);

  Smarts smarts3;
  smarts3.init("[r]");
  COMPARE(Smiley::AE_Cyclic, smarts3.trees().atom(0)->type);

  // AE_Acyclic
  Smarts smarts2;
  smarts2.init("[R0]");
  COMPARE(Smiley::AE_Acyclic, smarts2.trees().atom(0)->type);
}

void test_isotope()
{
  // AE_Isotope
  Smarts smarts1;
  smarts1.init("[13]");
  COMPARE(Smiley::AE_Isotope, smarts1.trees().atom(0)->type);
  COMPARE(13, smarts1.trees().atom(0)->value);
}

void test_atomic_number()
{
  // AE_AtomicNumber
  Smarts smarts1;
  smarts1.init("[#7]");
  COMPARE(Smiley::AE_AtomicNumber, smarts1.trees().atom(0)->type);
  COMPARE(7, smarts1.trees().atom(0)->value);
}

void test_degree()
{
  // AE_Degree
  Smarts smarts1;
  smarts1.init("[D]");
  COMPARE(Smiley::AE_Degree, smarts1.trees().atom(0)->type);
  COMPARE(1, smarts1.trees().atom(0)->value);

  Smarts smarts2;
  smarts2.init("[D3]");
  COMPARE(Smiley::AE_Degree, smarts2.trees().atom(0)->type);
  COMPARE(3, smarts2.trees().atom(0)->value);
}

void test_valence()
{
  // AE_Valence
  Smarts smarts1;
  smarts1.init("[v]");
  COMPARE(Smiley::AE_Valence, smarts1.trees().atom(0)->type);
  COMPARE(1, smarts1.trees().atom(0)->value);

  Smarts smarts2;
  smarts2.init("[v3]");
  COMPARE(Smiley::AE_Valence, smarts2.trees().atom(0)->type);
  COMPARE(3, smarts2.trees().atom(0)->value);
}

void test_connectivity()
{
  // AE_Connectivity
  Smarts smarts1;
  smarts1.init("[X]");
  COMPARE(Smiley::AE_Connectivity, smarts1.trees().atom(0)->type);
  COMPARE(1, smarts1.trees().atom(0)->value);

  Smarts smarts2;
  smarts2.init("[X3]");
  COMPARE(Smiley::AE_Connectivity, smarts2.trees().atom(0)->type);
  COMPARE(3, smarts2.trees().atom(0)->value);
}

void test_total_h()
{
  // AE_TotalH
  Smarts smarts1;
  smarts1.init("[H]");
  COMPARE(Smiley::AE_TotalH, smarts1.trees().atom(0)->type);
  COMPARE(1, smarts1.trees().atom(0)->value);

  Smarts smarts2;
  smarts2.init("[H3]");
  COMPARE(Smiley::AE_TotalH, smarts2.trees().atom(0)->type);
  COMPARE(3, smarts2.trees().atom(0)->value);
}

void test_implicit_h()
{
  // AE_ImplicitH
  Smarts smarts1;
  smarts1.init("[h]");
  COMPARE(Smiley::AE_ImplicitH, smarts1.trees().atom(0)->type);
  COMPARE(-1, smarts1.trees().atom(0)->value); // >= h1

  Smarts smarts2;
  smarts2.init("[h3]");
  COMPARE(Smiley::AE_ImplicitH, smarts2.trees().atom(0)->type);
  COMPARE(3, smarts2.trees().atom(0)->value);
}

void test_ring_membership()
{
  // AE_RingMembership
  Smarts smarts1;
  smarts1.init("[R1]");
  COMPARE(Smiley::AE_RingMembership, smarts1.trees().atom(0)->type);
  COMPARE(1, smarts1.trees().atom(0)->value);

  Smarts smarts2;
  smarts2.init("[R2]");
  COMPARE(Smiley::AE_RingMembership, smarts2.trees().atom(0)->type);
  COMPARE(2, smarts2.trees().atom(0)->value);
}

void test_ring_size()
{
  // AE_RingSize
  Smarts smarts1;
  smarts1.init("[r3]");
  COMPARE(Smiley::AE_RingSize, smarts1.trees().atom(0)->type);
  COMPARE(3, smarts1.trees().atom(0)->value);

  Smarts smarts2;
  smarts2.init("[r6]");
  COMPARE(Smiley::AE_RingSize, smarts2.trees().atom(0)->type);
  COMPARE(6, smarts2.trees().atom(0)->value);
}

void test_ring_connectivity()
{
  // AE_RingConnectivity
  Smarts smarts1;
  smarts1.init("[x]");
  COMPARE(Smiley::AE_RingConnectivity, smarts1.trees().atom(0)->type);
  COMPARE(-1, smarts1.trees().atom(0)->value); // >= x1

  Smarts smarts2;
  smarts2.init("[x3]");
  COMPARE(Smiley::AE_RingConnectivity, smarts2.trees().atom(0)->type);
  COMPARE(3, smarts2.trees().atom(0)->value);
}

void test_charge()
{
  // AE_Charge
  Smarts smarts1;
  smarts1.init("[-]");
  COMPARE(Smiley::AE_Charge, smarts1.trees().atom(0)->type);
  COMPARE(-1, smarts1.trees().atom(0)->value);

  Smarts smarts2;
  smarts2.init("[-3]");
  COMPARE(Smiley::AE_Charge, smarts2.trees().atom(0)->type);
  COMPARE(-3, smarts2.trees().atom(0)->value);

  Smarts smarts3;
  smarts3.init("[+]");
  COMPARE(Smiley::AE_Charge, smarts3.trees().atom(0)->type);
  COMPARE(1, smarts3.trees().atom(0)->value);

  Smarts smarts4;
  smarts4.init("[+2]");
  COMPARE(Smiley::AE_Charge, smarts4.trees().atom(0)->type);
  COMPARE(2, smarts4.trees().atom(0)->value);
}

void test_chirality()
{
  // AE_Chirality
}

void test_atom_class()
{
  // AE_AtomClass
  Smarts smarts1;
  smarts1.init("[:2]");
  COMPARE(Smiley::AE_AtomClass, smarts1.trees().atom(0)->type);
  COMPARE(2, smarts1.trees().atom(0)->value);

  Smarts smarts2;
  smarts2.init("[:5]");
  COMPARE(Smiley::AE_AtomClass, smarts2.trees().atom(0)->type);
  COMPARE(5, smarts2.trees().atom(0)->value);
}

void test_atom_not()
{
  Smarts smarts;
  smarts.init("[!C!N]");

  impl::SmartsAtomExpr *root = smarts.trees().atom(0);

  COMPARE(Smiley::OP_AndHi, root->type);

  COMPARE(Smiley::OP_Not, root->left->type);
  COMPARE(Smiley::OP_Not, root->right->type);

  COMPARE(Smiley::AE_AliphaticElement, root->left->arg->type);
  COMPARE(Smiley::AE_AliphaticElement, root->right->arg->type);
}

void test_atom_and_or_complex()
{
  Smarts smarts;
  smarts.init("[13C,13C;13C,13C]");

  impl::SmartsAtomExpr *root = smarts.trees().atom(0);

  COMPARE(Smiley::OP_AndLo, root->type);

  COMPARE(Smiley::OP_Or, root->left->type);
  COMPARE(Smiley::OP_Or, root->right->type);

  COMPARE(Smiley::OP_AndHi, root->left->left->type);
  COMPARE(Smiley::OP_AndHi, root->left->right->type);
  COMPARE(Smiley::OP_AndHi, root->right->left->type);
  COMPARE(Smiley::OP_AndHi, root->right->right->type);

  COMPARE(Smiley::AE_AliphaticElement, root->left->left->left->type);
  COMPARE(Smiley::AE_Isotope, root->left->left->right->type);
  COMPARE(Smiley::AE_AliphaticElement, root->left->right->left->type);
  COMPARE(Smiley::AE_Isotope, root->left->right->right->type);
  COMPARE(Smiley::AE_AliphaticElement, root->right->left->left->type);
  COMPARE(Smiley::AE_Isotope, root->right->left->right->type);
  COMPARE(Smiley::AE_AliphaticElement, root->right->right->left->type);
  COMPARE(Smiley::AE_Isotope, root->right->right->right->type);
}

void test_implicit_bonds()
{
  // implicit BE_Single
  Smarts smarts1;
  smarts1.init("CC");
  COMPARE(Smiley::BE_Single, smarts1.trees().bond(0)->type);

  // implicit BE_Single or BE_Aromatic
  Smarts smarts2;
  smarts2.init("cc");
  impl::SmartsBondExpr *root = smarts2.trees().bond(0);
  COMPARE(Smiley::OP_Or, root->type);
  COMPARE(Smiley::BE_Single, root->left->type);
  COMPARE(Smiley::BE_Aromatic, root->right->type);
}

void test_bonds()
{
  // BE_Single
  Smarts smarts1;
  smarts1.init("C-C");
  COMPARE(Smiley::BE_Single, smarts1.trees().bond(0)->type);

  // BE_Double
  Smarts smarts2;
  smarts2.init("C=C");
  COMPARE(Smiley::BE_Double, smarts2.trees().bond(0)->type);

  // BE_Triple
  Smarts smarts3;
  smarts3.init("C#C");
  COMPARE(Smiley::BE_Triple, smarts3.trees().bond(0)->type);

  // BE_Quadriple
  Smarts smarts4;
  smarts4.init("C$C");
  COMPARE(Smiley::BE_Quadriple, smarts4.trees().bond(0)->type);

  // BE_Aromatic
  Smarts smarts5;
  smarts5.init("c:c");
  COMPARE(Smiley::BE_Aromatic, smarts5.trees().bond(0)->type);

  // BE_Up
  Smarts smarts6;
  smarts6.init("C/C");
  COMPARE(Smiley::BE_Up, smarts6.trees().bond(0)->type);

  // BE_Down
  Smarts smarts7;
  smarts7.init("C\\C");
  COMPARE(Smiley::BE_Down, smarts7.trees().bond(0)->type);

  // BE_Ring
  Smarts smarts8;
  smarts8.init("C@C");
  COMPARE(Smiley::BE_Ring, smarts8.trees().bond(0)->type);

  // BE_Any
  Smarts smarts9;
  smarts9.init("C~C");
  COMPARE(Smiley::BE_Any, smarts9.trees().bond(0)->type);

  // multiple bonds
  Smarts smarts10;
  smarts10.init("C-C=C#C");
  COMPARE(Smiley::BE_Single, smarts10.trees().bond(0)->type);
  COMPARE(Smiley::BE_Double, smarts10.trees().bond(1)->type);
  COMPARE(Smiley::BE_Triple, smarts10.trees().bond(2)->type);

  Smarts smarts11;
  smarts11.init("[C]-[C]=[C]#[C]");
  COMPARE(Smiley::BE_Single, smarts11.trees().bond(0)->type);
  COMPARE(Smiley::BE_Double, smarts11.trees().bond(1)->type);
  COMPARE(Smiley::BE_Triple, smarts11.trees().bond(2)->type);
}

void test_bond_not()
{
  Smarts smarts;
  smarts.init("C!-!=C");

  impl::SmartsBondExpr *root = smarts.trees().bond(0);

  COMPARE(Smiley::OP_AndHi, root->type);

  COMPARE(Smiley::OP_Not, root->left->type);
  COMPARE(Smiley::OP_Not, root->right->type);

  COMPARE(Smiley::BE_Double, root->left->arg->type);
  COMPARE(Smiley::BE_Single, root->right->arg->type);
}

void test_bond_and_or_complex()
{
  Smarts smarts;
  smarts.init("C-@,=@C");

  impl::SmartsBondExpr *root = smarts.trees().bond(0);

  COMPARE(Smiley::OP_Or, root->type);

  COMPARE(Smiley::OP_AndHi, root->left->type);
  COMPARE(Smiley::OP_AndHi, root->right->type);

  COMPARE(Smiley::BE_Ring, root->left->left->type);
  COMPARE(Smiley::BE_Double, root->left->right->type);
  COMPARE(Smiley::BE_Ring, root->right->left->type);
  COMPARE(Smiley::BE_Single, root->right->right->type);
}

void test_bond_ring()
{
  Smarts smarts1;
  smarts1.init("C1CC-,=1");
  impl::SmartsBondExpr *root = smarts1.trees().bond(2);
  COMPARE(Smiley::OP_Or, root->type);
  COMPARE(Smiley::BE_Double, root->left->type);
  COMPARE(Smiley::BE_Single, root->right->type);

  Smarts smarts2;
  smarts2.init("C-,=1CC1");
  root = smarts2.trees().bond(2);
  COMPARE(Smiley::OP_Or, root->type);
  COMPARE(Smiley::BE_Double, root->left->type);
  COMPARE(Smiley::BE_Single, root->right->type);
}

void test_smarts_match(const std::string &smarts, const std::string &smiles, bool expected = true)
{
  std::cout << "Testing: " << smarts << " in " << smiles << std::endl;
  Smarts s;
  if (!s.init(smarts)) {
    std::cerr << s.error().what();
    return;
  }

  HeMol mol = hemol_from_smiles(smiles);
  RingSet<HeMol> rings = relevant_cycles(mol);

  COMPARE(expected, s.search(mol, rings));
}

void test_simple_atom_match()
{
  // any
  test_smarts_match("*", "C");
  test_smarts_match("*", "N");
  test_smarts_match("*", "O");
  test_smarts_match("[*]", "C");
  test_smarts_match("[*]", "N");
  test_smarts_match("[*]", "O");

  // organic subset: aromatic/aliphatic element
  test_smarts_match("C", "C");
  test_smarts_match("C", "N", false);
  test_smarts_match("c", "c");
  test_smarts_match("c", "n", false);

  // aromatic/aliphatic element
  test_smarts_match("[C]", "[C]");
  test_smarts_match("[C]", "[N]", false);
  test_smarts_match("[c]", "[c]");
  test_smarts_match("[c]", "[n]", false);

  // aromatic
  test_smarts_match("a", "c");
  test_smarts_match("a", "C", false);
  test_smarts_match("[a]", "c");
  test_smarts_match("[a]", "C", false);

  // aliphatic
  test_smarts_match("A", "C");
  test_smarts_match("A", "c", false);
  test_smarts_match("[A]", "C");
  test_smarts_match("[A]", "c", false);

  // cyclic
  test_smarts_match("[R]", "C1CC1");
  test_smarts_match("[R]", "CCC", false);
  test_smarts_match("[r]", "C1CC1");
  test_smarts_match("[r]", "CCC", false);

  // acyclic
  test_smarts_match("[R0]", "CCC");
  test_smarts_match("[R0]", "C1CC1", false);

  // isotope
  test_smarts_match("[13]", "[13C]");
  test_smarts_match("[13]", "C", false);
  test_smarts_match("[12]", "C");

  // atomic number
  test_smarts_match("[#6]", "C");
  test_smarts_match("[#6]", "c");
  test_smarts_match("[#6]", "N", false);

  // degree
  test_smarts_match("[D]", "CC"); //default: exactly 1
  test_smarts_match("[D]", "C", false);
  test_smarts_match("[D0]", "C");
  test_smarts_match("[D1]", "C", false);
  test_smarts_match("[D1]", "CC");
  test_smarts_match("[D2]", "CC", false);
  test_smarts_match("[D2]", "CCC");
  test_smarts_match("[D3]", "CCC", false);
  test_smarts_match("[D3]", "CC(C)C");
  test_smarts_match("[D4]", "CC(C)C", false);

  // valence
  test_smarts_match("[v]", "CC"); // default: exactly 1
  test_smarts_match("[v]", "C", false);
  test_smarts_match("[v0]", "C");
  test_smarts_match("[v1]", "CC");
  test_smarts_match("[v1]", "C", false);
  test_smarts_match("[v2]", "CC", false);
  test_smarts_match("[v2]", "C=C");
  test_smarts_match("[v2]", "CCC");
  test_smarts_match("[v3]", "CCC", false);
  test_smarts_match("[v3]", "CC=C");
  test_smarts_match("[v3]", "CC(C)C");
  test_smarts_match("[v4]", "CC(C)C", false);
  test_smarts_match("[v4]", "CC=C=C");

  // connectivity
  test_smarts_match("[X]", "Cl"); // default: exactly 1
  test_smarts_match("[X]", "C", false);
  test_smarts_match("[X4]", "C");
  test_smarts_match("[X3]", "C", false);
  test_smarts_match("[X4]", "C=C");
  test_smarts_match("[X4]", "C([H])[H]");
  test_smarts_match("[X3]", "C([H])[H]", false);
  test_smarts_match("[X1]", "C([H])[H]");
  test_smarts_match("[X3]", "N");
  test_smarts_match("[X2]", "N", false);

  // total H
  test_smarts_match("[H]", "[CH]"); // default: exactly 1
  test_smarts_match("[H]", "[CH2]", false); // default: exactly 1
  test_smarts_match("[H4]", "C"); // all implicit
  test_smarts_match("[H3]", "C", false);
  test_smarts_match("[H4]", "C([H])([H])([H])[H]"); // all explicit
  test_smarts_match("[H3]", "C([H])([H])([H])[H]", false);
  test_smarts_match("[H4]", "C([H])[H]"); // implicit + explicit
  test_smarts_match("[H3]", "C([H])[H]", false);
  test_smarts_match("[H4]", "[CH2]([H])[H]"); // implicit + explicit
  test_smarts_match("[H5]", "[CH2]([H])[H]", false);
  test_smarts_match("[H3]", "[CH3-]"); // implicit
  test_smarts_match("[H2]", "[CH3-]", false);
  test_smarts_match("[H2]", "C=C"); // implicit
  test_smarts_match("[H1]", "C=C", false);
  test_smarts_match("[H1]", "C#C"); // implicit
  test_smarts_match("[H0]", "C#C", false);

  // implicit H
  test_smarts_match("[h]", "[CH]"); // default: at least 1
  test_smarts_match("[h]", "[CH2]");
  test_smarts_match("[h]", "[CH3]");
  test_smarts_match("[h]", "[C]", false);
  test_smarts_match("[h4]", "C"); // all implicit
  test_smarts_match("[h0]", "C([H])([H])([H])[H]"); // all explicit
  test_smarts_match("[h2]", "C([H])[H]"); // implicit + explicit
  test_smarts_match("[h2]", "[CH2]([H])[H]"); // implicit + explicit
  test_smarts_match("[h3]", "[CH3-]"); // implicit
  test_smarts_match("[h2]", "C=C"); // implicit
  test_smarts_match("[h1]", "C#C"); // implicit

  // ring membership
  test_smarts_match("[R1]", "C1CC1");
  test_smarts_match("[R1]", "CCC", false);
  test_smarts_match("[R2]", "C12CC1CC2");
  test_smarts_match("[R2]", "C1CC1", false);
  test_smarts_match("[R3]", "C123CCC1CCC2CCC3");
  test_smarts_match("[R3]", "C12CC1CC2", false);

  // ring size
  test_smarts_match("[r3]", "C1CC1");
  test_smarts_match("[r4]", "C1CCC1");
  test_smarts_match("[r5]", "C1CCCC1");
  test_smarts_match("[r6]", "C1CCCCC1");
  test_smarts_match("[r4]", "C1CC1", false);
  test_smarts_match("[r5]", "C1CCC1", false);
  test_smarts_match("[r6]", "C1CCCC1", false);
  test_smarts_match("[r7]", "C1CCCCC1", false);

  // ring connectivity
  test_smarts_match("[x]", "C1CC1"); // default: at least 1
  test_smarts_match("[x]", "C12CC1CC2");
  test_smarts_match("[x]", "CCC", false);
  test_smarts_match("[x2]", "C1CC1");
  test_smarts_match("[x1]", "C1CC1", false);
  test_smarts_match("[x3]", "C1CC1", false);
  test_smarts_match("[x1]", "C12CC1CC2", false);
  test_smarts_match("[x2]", "C12CC1CC2");
  test_smarts_match("[x3]", "C12CC1CC2");

  // charge
  test_smarts_match("[-]", "[C-]");
  test_smarts_match("[-]", "C", false);
  test_smarts_match("[-]", "[C-2]", false);
  test_smarts_match("[+]", "[C+]");
  test_smarts_match("[+]", "C", false);
  test_smarts_match("[+]", "[C+2]", false);
  test_smarts_match("[-2]", "[C-2]");
  test_smarts_match("[-2]", "C", false);
  test_smarts_match("[-2]", "[C-]", false);
  test_smarts_match("[+2]", "[C+2]");
  test_smarts_match("[+2]", "C", false);
  test_smarts_match("[+2]", "[C+]", false);
}

void test_atom_operators()
{
  // implicit OP_AndHi
  test_smarts_match("[CD0]", "C");
  test_smarts_match("[ND0]", "C", false);
  test_smarts_match("[CD1]", "C", false);

  // explicit OP_AndHi
  test_smarts_match("[C&D0]", "C");
  test_smarts_match("[N&D0]", "C", false);
  test_smarts_match("[C&D1]", "C", false);

  // OP_AndLo
  test_smarts_match("[C;D0]", "C");
  test_smarts_match("[N;D0]", "C", false);
  test_smarts_match("[C;D1]", "C", false);

  // OP_Or
  test_smarts_match("[C,N]", "C");
  test_smarts_match("[C,N]", "N");
  test_smarts_match("[C,N]", "O", false);

  // OP_Not
  test_smarts_match("[!C]", "N");
  test_smarts_match("[!C]", "O");
  test_smarts_match("[!C]", "C", false);

  // combination
  test_smarts_match("[C,N;+]", "[C+]");
  test_smarts_match("[C,N;+]", "[N+]");
  test_smarts_match("[C,N;+]", "[C]", false);
  test_smarts_match("[C,N;+]", "[N]", false);
  test_smarts_match("[C,N;+]", "[O+]", false);
  test_smarts_match("[C,N;+]", "[O+]", false);

  test_smarts_match("[!C!N]", "O");
  test_smarts_match("[!C!N]", "C", false);
  test_smarts_match("[!C!N]", "N", false);
}

void test_simple_bond_match()
{
  // implicit single
  test_smarts_match("CC", "CC");
  test_smarts_match("CC", "cc", false);
  test_smarts_match("CC", "C=C", false);

  // implicit aromatic
  test_smarts_match("cc", "cc");
  test_smarts_match("cc", "CC", false);
  test_smarts_match("cc", "C=C", false);

  // any
  test_smarts_match("*~*", "CC");
  test_smarts_match("*~*", "C-C");
  test_smarts_match("*~*", "cc");
  test_smarts_match("*~*", "c:c");
  test_smarts_match("*~*", "C=C");
  test_smarts_match("*~*", "C#C");

  // single
  test_smarts_match("C-C", "CC");
  test_smarts_match("C-C", "C-C");
  test_smarts_match("C-C", "C=C", false);
  test_smarts_match("C-C", "cc", false);

  // double
  test_smarts_match("C=C", "C=C");
  test_smarts_match("C=C", "CC", false);
  test_smarts_match("C=C", "cc", false);

  // triple
  test_smarts_match("C#C", "C#C");
  test_smarts_match("C#C", "CC", false);
  test_smarts_match("C#C", "cc", false);

  // quadriple

  // aromatic
  test_smarts_match("*:*", "cc");
  test_smarts_match("*:*", "cC", false);
  test_smarts_match("*:*", "Cc", false);

  // up
  // down

  // ring
  test_smarts_match("*@*", "C1CC1");
  test_smarts_match("*@*", "CCC", false);
}

void test_bond_operators()
{
  // OP_Not
  test_smarts_match("*!-*", "C=C");
  test_smarts_match("*!-*", "C-C", false);

  // implicit OP_AndHi
  test_smarts_match("*-@*", "C1CC1");
  test_smarts_match("*-@*", "CCC", false);

  // OP_AndHi
  test_smarts_match("*-&@*", "C1CC1");
  test_smarts_match("*-&@*", "CCC", false);

  // OP_AndLo
  test_smarts_match("*-;@*", "C1CC1");
  test_smarts_match("*-;@*", "CCC", false);

  // OP_Or
  test_smarts_match("*-,=*", "C-C");
  test_smarts_match("*-,=*", "C=C");
  test_smarts_match("*-,=*", "C#C", false);
}

void test_daylight_examples()
{
  // cc           any pair of attached aromatic carbons
  test_smarts_match("cc", "cc");
  test_smarts_match("cc", "c-c");
  test_smarts_match("cc", "c:c");

  // c:c         aromatic carbons joined by an aromatic bond
  test_smarts_match("c:c", "cc");
  test_smarts_match("c:c", "c-c", false);
  test_smarts_match("c:c", "c:c");

  // c-c         aromatic carbons joined by a single bond (e.g. biphenyl).
  test_smarts_match("c-c", "cc", false);
  test_smarts_match("c-c", "c-c");
  test_smarts_match("c-c", "c:c", false);

  // O           any aliphatic oxygen
  test_smarts_match("O", "O");
  test_smarts_match("O", "[O]");
  test_smarts_match("O", "CO");
  test_smarts_match("O", "C[OH]");
  test_smarts_match("O", "C[O-]");
  test_smarts_match("O", "N", false);

  // [O;H1]      simple hydroxy oxygen
  test_smarts_match("[O;H1]", "O", false);
  test_smarts_match("[O;H1]", "[O]", false);
  test_smarts_match("[O;H1]", "CO");
  test_smarts_match("[O;H1]", "C[OH]");
  test_smarts_match("[O;H1]", "C[O-]", false);
  test_smarts_match("[O;H1]", "N", false);

  // [O;D1]      1-connected (hydroxy or hydroxide) oxygen
  test_smarts_match("[O;D1]", "O", false);
  test_smarts_match("[O;D1]", "[O]", false);
  test_smarts_match("[O;D1]", "CO");
  test_smarts_match("[O;D1]", "C[OH]");
  test_smarts_match("[O;D1]", "C[O-]");
  test_smarts_match("[O;D1]", "N", false);

  // [O;D2]      2-connected (etheric) oxygen
  test_smarts_match("[O;D2]", "COC");
  test_smarts_match("[O;D2]", "O", false);
  test_smarts_match("[O;D2]", "[O]", false);
  test_smarts_match("[O;D2]", "CO", false);
  test_smarts_match("[O;D2]", "C[OH]", false);
  test_smarts_match("[O;D2]", "C[O-]", false);
  test_smarts_match("[O;D2]", "N", false);

  // [C,c]       any carbon
  test_smarts_match("[C,c]", "C");
  test_smarts_match("[C,c]", "c");
  test_smarts_match("[C,c]", "N", false);
  test_smarts_match("[C,c]", "n", false);

  // [F,Cl,Br,I]  the 1st four halogens.
  test_smarts_match("[F,Cl,Br,I]", "F");
  test_smarts_match("[F,Cl,Br,I]", "Cl");
  test_smarts_match("[F,Cl,Br,I]", "Br");
  test_smarts_match("[F,Cl,Br,I]", "I");
  test_smarts_match("[F,Cl,Br,I]", "[F]");
  test_smarts_match("[F,Cl,Br,I]", "[Cl]");
  test_smarts_match("[F,Cl,Br,I]", "[Br]");
  test_smarts_match("[F,Cl,Br,I]", "[I]");
  test_smarts_match("[F,Cl,Br,I]", "CCF");
  test_smarts_match("[F,Cl,Br,I]", "CC", false);

  // [N;R]       must be aliphatic nitrogen AND in a ring
  test_smarts_match("[N;R]", "N1CCCCC1");
  test_smarts_match("[N;R]", "C1CCCCC1", false);
  test_smarts_match("[N;R]", "NCC", false);
  test_smarts_match("[N;R]", "n1ccccc1", false);

  // [!C;R]      ( NOT aliphatic carbon ) AND in a ring
  test_smarts_match("[!C;R]", "N1CC1");
  test_smarts_match("[!C;R]", "NCC", false);
  test_smarts_match("[!C;R]", "CCC", false);
  test_smarts_match("[!C;R]", "C1CC1", false);

  // [n;H1]      H-pyrrole nitrogen
  test_smarts_match("[n;H1]", "[nH]1cccc1");
  test_smarts_match("[n;H1]", "c1cccc1", false);
  test_smarts_match("[n;H1]", "o1cccc1", false);

  // [n&H1]      same as above
  test_smarts_match("[n;H1]", "[nH]1cccc1");
  test_smarts_match("[n;H1]", "c1cccc1", false);
  test_smarts_match("[n;H1]", "o1cccc1", false);

  // [c,n&H1]    any arom carbon OR H-pyrrole nitrogen
  test_smarts_match("[c,n&H1]", "[nH]1cccc1");
  test_smarts_match("[c,n&H1]", "c1cccc1");
  test_smarts_match("[c,n&H1]", "o1cccc1");

  // [c,n;H1]    (arom carbon OR arom nitrogen) and exactly one H
  test_smarts_match("[c,n;H1]", "[nH]1cccc1");
  test_smarts_match("[c,n;H1]", "c1cccc1");
  test_smarts_match("[c,n;H1]", "o1cccc1");
  test_smarts_match("[c,n;H1]", "c1(F)c(F)c(F)c(F)c1F", false);

  // *!@*        two atoms connected by a non-ringbond
  test_smarts_match("*!@*", "CC");
  test_smarts_match("*!@*", "CO");
  test_smarts_match("*!@*", "C1CC1", false);

  // *@;!:*      two atoms connected by a non-aromatic ringbond
  test_smarts_match("*@;!:*", "CC", false);
  test_smarts_match("*@;!:*", "CO", false);
  test_smarts_match("*@;!:*", "C1CC1");
  test_smarts_match("*@;!:*", "c1ccccc1", false);

  // [C,c]=,#[C,c]       two carbons connected by a double or triple bond
  test_smarts_match("[C,c]=,#[C,c]", "C=C");
  test_smarts_match("[C,c]=,#[C,c]", "C#C");
  test_smarts_match("[C,c]=,#[C,c]", "c=c");
  test_smarts_match("[C,c]=,#[C,c]", "c#c");
  test_smarts_match("[C,c]=,#[C,c]", "C-C", false);
  test_smarts_match("[C,c]=,#[C,c]", "C=O", false);
  test_smarts_match("[C,c]=,#[C,c]", "C-c", false);
  test_smarts_match("[C,c]=,#[C,c]", "C=O", false);
}

void test_complex()
{
  test_smarts_match("C(=O)N", "c1cc(CC(=O)NC)ccc1");
}

void test_disconnected()
{
  test_smarts_match("C.O", "CO");
  test_smarts_match("C.O", "CC", false);
}

void test_disconnected_mapping1()
{
  // index:   01 23 45
  // SMARTS:  CO.CN.CF
  // SMILES:  CN.CF.CO

  std::cout << "Testing: CO.CN.CF in CN.CF.CO" << std::endl;
  Smarts s;
  s.init("CO.CN.CF");
  HeMol mol = hemol_from_smiles("CN.CF.CO");
  RingSet<HeMol> rings = relevant_cycles(mol);

  // perform match
  MappingList mappings;
  ASSERT(s.search(mol, mappings, rings));

  // check mapping
  COMPARE(1, mappings.maps.size());
  const IsomorphismMapping &map = mappings.maps[0];

  COMPARE(4, map[0]);
  COMPARE(5, map[1]);
  COMPARE(0, map[2]);
  COMPARE(1, map[3]);
  COMPARE(2, map[4]);
  COMPARE(3, map[5]);
}

void test_disconnected_mapping2()
{
  // index:   01 23 45
  // SMARTS:  CO.CO.CC
  // SMILES:  OCCCCCO
  // index:   0123456

  std::cout << "Testing: CO.CO.CC in OCCCCO" << std::endl;
  Smarts s;
  s.init("CO.CO.CC");
  HeMol mol = hemol_from_smiles("OCCCCCO");
  RingSet<HeMol> rings = relevant_cycles(mol);

  // perform match
  MappingList mappings;
  ASSERT(s.search(mol, mappings, rings));

  // check mapping
  COMPARE(4, mappings.maps.size());

  const IsomorphismMapping &map1 = mappings.maps[0];
  COMPARE(1, map1[0]);
  COMPARE(0, map1[1]);
  COMPARE(5, map1[2]);
  COMPARE(6, map1[3]);
  COMPARE(2, map1[4]);
  COMPARE(3, map1[5]);

  const IsomorphismMapping &map2 = mappings.maps[1];
  COMPARE(1, map2[0]);
  COMPARE(0, map2[1]);
  COMPARE(5, map2[2]);
  COMPARE(6, map2[3]);
  COMPARE(3, map2[4]);
  COMPARE(4, map2[5]);

  const IsomorphismMapping &map3 = mappings.maps[2];
  COMPARE(5, map3[0]);
  COMPARE(6, map3[1]);
  COMPARE(1, map3[2]);
  COMPARE(0, map3[3]);
  COMPARE(2, map3[4]);
  COMPARE(3, map3[5]);

  const IsomorphismMapping &map4 = mappings.maps[3];
  COMPARE(5, map4[0]);
  COMPARE(6, map4[1]);
  COMPARE(1, map4[2]);
  COMPARE(0, map4[3]);
  COMPARE(3, map4[4]);
  COMPARE(4, map4[5]);
}

void test_recursive_parsing()
{
  Smarts smarts;
  smarts.init("[R,$(CCO),N]");

  COMPARE(1, num_atoms(smarts.query()));
  COMPARE(1, smarts.trees().atoms().size());
  COMPARE(1, smarts.recursiveMols().size());
  COMPARE(1, smarts.recursiveTrees().size());

  impl::SmartsAtomExpr *root = smarts.trees().atom(0);
  COMPARE(Smiley::OP_Or, root->type);
  COMPARE(Smiley::AE_Cyclic, root->right->type);
  COMPARE(Smiley::OP_Or, root->left->type);
  COMPARE(Smiley::AE_AliphaticElement, root->left->left->type);
  COMPARE(Smiley::AE_Recursive, root->left->right->type);
  COMPARE(0, root->left->right->value);

  COMPARE(3, num_atoms(smarts.recursiveMols()[0]));
  COMPARE(3, smarts.recursiveTrees()[0].atoms().size());
}

void test_recursive()
{
  test_smarts_match("[$(*O);$(*N)]", "OCN");
  test_smarts_match("[$(*O);$(*N)]", "OCO", false);
  test_smarts_match("[$(*O);$(*N)]", "NCN", false);

  // CaaO  C ortho to O
  test_smarts_match("CaaO", "Cc1c(O)cccc1"); // ortho
  test_smarts_match("CaaO", "Cc1cc(O)ccc1", false); // meta
  test_smarts_match("CaaO", "Cc1ccc(O)cc1", false); // para

  // CaaaN       C meta to N
  test_smarts_match("CaaaN", "Cc1c(N)cccc1", false); // ortho
  test_smarts_match("CaaaN", "Cc1cc(N)ccc1"); // meta
  test_smarts_match("CaaaN", "Cc1ccc(N)cc1", false); // para

  // Caa(O)aN    C ortho to O and meta to N (but 2O,3N only)
  test_smarts_match("Caa(O)aN", "Cc1c(O)c(N)ccc1"); // 2O,3N
  test_smarts_match("Caa(O)aN", "Cc1c(O)ccc(N)c1", false); // 2O,5N

  // Ca(aO)aaN   C ortho to O and meta to N (but 2O,5N only)
  test_smarts_match("Ca(aO)aaN", "Cc1c(O)c(N)ccc1", false); // 2O,3N
  test_smarts_match("Ca(aO)aaN", "Cc1c(O)ccc(N)c1"); // 2O,5N

  // C[$(aaO);$(aaaN)]   C ortho to O and meta to N (all cases)
  test_smarts_match("C[$(aaO);$(aaaN)]", "Cc1c(O)c(N)ccc1"); // 2O,3N
  test_smarts_match("C[$(aaO);$(aaaN)]", "Cc1c(O)ccc(N)c1"); // 2O,5N
}

int main()
{
  //
  // Parsing
  //

  test_parser_error();

  // atoms
  test_organic_subset();
  test_cyclic();
  test_isotope();
  test_atomic_number();
  test_degree();
  test_valence();
  test_connectivity();
  test_total_h();
  test_implicit_h();
  test_ring_membership();
  test_ring_size();
  test_ring_connectivity();
  test_charge();
  test_chirality();
  test_atom_class();
  test_atom_not();
  test_atom_and_or_complex();

  // bonds
  test_implicit_bonds();
  test_bonds();
  test_bond_not();
  test_bond_and_or_complex();
  test_bond_ring();

  //
  // Matching
  //

  // atoms
  test_simple_atom_match();
  test_atom_operators();

  // bonds
  test_simple_bond_match();
  test_bond_operators();


  //
  // Complex matching
  //
  test_daylight_examples();
  test_complex();

  //
  // Disconnected SMARTS
  //
  test_disconnected();
  test_disconnected_mapping1();
  test_disconnected_mapping2();

  //
  // Recursive SMARTS
  //
  test_recursive_parsing();
  test_recursive();

  //print_smarts("[A$(CCC)R]");
}
