#include <Helium/algorithms/cycles.h>
#include <Helium/smiles.h>
#include <Helium/fileio/moleculefile.h>

#include "cycles.h"

using namespace Helium;

void test_cyclomatic_number(const std::string &smiles, unsigned int expected)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol = hemol_from_smiles(smiles);
  COMPARE(expected, cyclomatic_number(mol));
}

void test_cycle_membership(const HeMol &mol)
{
  std::vector<bool> cyclic_atoms, cyclic_bonds;
  cycle_membership(mol, cyclic_atoms, cyclic_bonds);
  RingSet<HeMol> rings = relevant_cycles(mol);

  //std::cout << write_smiles(mol, WriteSmiles::Order) << std::endl;

  HeMol::atom_iter atom, end_atom;
  TIE(atom, end_atom) = get_atoms(mol);
  for (; atom != end_atom; ++atom) {
    if (rings.isAtomInRing(*atom))
      COMPARE(true, cyclic_atoms[get_index(mol, *atom)]);
    else
      COMPARE(false, cyclic_atoms[get_index(mol, *atom)]);
  }

  HeMol::bond_iter bond, end_bond;
  TIE(bond, end_bond) = get_bonds(mol);
  for (; bond != end_bond; ++bond) {
    if (rings.isBondInRing(*bond))
      COMPARE(true, cyclic_bonds[get_index(mol, *bond)]);
    else
      COMPARE(false, cyclic_bonds[get_index(mol, *bond)]);
  }
}

void test_cycle_membership(const std::string &filename)
{
  std::cout << "Testing cycle_membership()..." << std::endl;
  MoleculeFile file(filename);

  HeMol mol;
  for (unsigned int i = 0; i < file.numMolecules(); ++i) {
    file.readMolecule(mol);
    test_cycle_membership(mol);
  }
}

void test_relevant_cycles(const std::string &smiles, std::vector<std::pair<unsigned int, unsigned int> > &expected)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol = hemol_from_smiles(smiles);

  RingSet<HeMol> cycles = relevant_cycles(mol);
  std::map<unsigned int, unsigned int> cycleSizeCounts;
  for (std::size_t i = 0; i < cycles.size(); ++i)
    cycleSizeCounts[cycles.ring(i).size()]++;

  for (std::size_t i = 0; i < expected.size(); ++i)
    COMPARE(expected[i].second, cycleSizeCounts[expected[i].first]);
}

struct IsomorphismCycleAlgorithm
{
  template<typename MoleculeType>
  std::vector<TestCycle> operator()(const MoleculeType &mol) const
  {
    RingSet<HeMol> cycles = relevant_cycles(mol);

    std::vector<TestCycle> result;
    for (std::size_t i = 0; i < cycles.size(); ++i) {
      TestCycle cycle(num_bonds(mol));
      for (std::size_t j = 0; j < cycles.ring(i).size(); ++j)
        cycle.cycle()[get_index(mol, cycles.ring(i).bond(j))] = true;
      result.push_back(cycle);
    }

    return result;
  }
};

void test_cycle_bit_matrix1()
{
  impl::CycleBitMatrix m(4);
  m.addRow();
  m.addRow();
  m.addRow();

  m.set(0, 0, true);
  m.set(0, 1, true);
  m.set(1, 2, true);
  m.set(1, 3, true);
  m.set(2, 0, true);
  m.set(2, 1, true);
  m.set(2, 2, true);
  m.set(2, 3, true);

  std::cout << m << std::endl;

  COMPARE(2, m.eliminate());

  std::cout << m << std::endl;
}

void test_cycle_bit_matrix2()
{
  impl::CycleBitMatrix m(5);
  m.addRow();
  m.addRow();
  m.addRow();

  m.set(0, 2, true);
  m.set(0, 3, true);
  m.set(0, 4, true);
  m.set(1, 0, true);
  m.set(1, 1, true);
  m.set(1, 2, true);
  m.set(2, 0, true);
  m.set(2, 1, true);
  m.set(2, 3, true);
  m.set(2, 4, true);

  std::cout << m << std::endl;

  COMPARE(2, m.eliminate());

  std::cout << m << std::endl;
}

void test_cycle_bit_matrix3()
{
  impl::CycleBitMatrix m(6);
  m.addRow();
  m.addRow();
  m.addRow();

  m.set(0, 0, true);
  m.set(0, 1, true);
  m.set(0, 2, true);
  m.set(0, 3, true);

  m.set(1, 0, true);
  m.set(1, 1, true);
  m.set(1, 4, true);
  m.set(1, 5, true);

  m.set(2, 2, true);
  m.set(2, 3, true);
  m.set(2, 4, true);
  m.set(2, 5, true);

  std::cout << m << std::endl;

  COMPARE(2, m.eliminate());

  std::cout << m << std::endl;
}

int main()
{
  // counter example to relevant_cycles
  HeMol mol = hemol_from_smiles("C(NC1CN2CCC1CC2)CN3CCC5(CCC3)NCc4ccccc4O5");
  test_cycle_membership(mol);


  test_cycle_bit_matrix1();
  test_cycle_bit_matrix2();
  test_cycle_bit_matrix3();


  // test cyclomatc number
  test_cyclomatic_number("CCC", 0);
  test_cyclomatic_number("C1CC1", 1);
  test_cyclomatic_number("C1CC1C1CC1", 2);
  test_cyclomatic_number("C1CC1.C1CC1", 2);

  test_cycle_membership(datadir() + "100K.hel");

  std::vector<std::pair<unsigned int, unsigned int> > cycles;

  // test relevant_cycles
  cycles.clear();
  test_relevant_cycles("CCC", cycles);

  cycles.clear();
  cycles.push_back(std::make_pair(3, 1));
  test_relevant_cycles("C1CC1", cycles);

  cycles.clear();
  cycles.push_back(std::make_pair(3, 2));
  test_relevant_cycles("C1CC1.C1CC1", cycles);

  cycles.clear();
  cycles.push_back(std::make_pair(3, 1));
  cycles.push_back(std::make_pair(4, 1));
  test_relevant_cycles("C1CC1.C1CCC1", cycles);

  cycles.clear();
  cycles.push_back(std::make_pair(3, 1));
  cycles.push_back(std::make_pair(6, 1));
  test_relevant_cycles("C1C(N)C1Cc1ccc(O)cc1", cycles);

  test_cycle_perception(IsomorphismCycleAlgorithm());
}
