#include <Helium/algorithms/cycles.h>
#include <Helium/smiles.h>
#include <Helium/fileio/molecules.h>

#include "test.h"

using namespace Helium;

void test_cyclomatic_number(const std::string &smiles, unsigned int expected)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol;
  parse_smiles(smiles, mol);
  COMPARE(expected, cyclomatic_number(mol));
}

void test_cycle_membership(const std::string &filename)
{
  std::cout << "Testing cycle_membership()..." << std::endl;
  MoleculeFile file(filename);

  HeMol mol;
  for (unsigned int i = 0; i < file.numMolecules(); ++i) {
    file.read_molecule(mol);
    std::vector<bool> cyclic_atoms, cyclic_bonds;
    cycle_membership(mol, cyclic_atoms, cyclic_bonds);

    HeMol::atom_iter atom, end_atom;
    tie(atom, end_atom) = get_atoms(mol);
    for (; atom != end_atom; ++atom) {
      if (is_cyclic(mol, *atom))
        COMPARE(true, cyclic_atoms[get_index(mol, *atom)]);
      else
        COMPARE(false, cyclic_atoms[get_index(mol, *atom)]);
    }

    HeMol::bond_iter bond, end_bond;
    tie(bond, end_bond) = get_bonds(mol);
    for (; bond != end_bond; ++bond) {
      if (is_cyclic(mol, *bond))
        COMPARE(true, cyclic_bonds[get_index(mol, *bond)]);
      else
        COMPARE(false, cyclic_bonds[get_index(mol, *bond)]);
    }
  }
}
 
void test_relevant_cycles(const std::string &smiles, std::vector<std::pair<unsigned int, unsigned int> > &expected)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol;
  parse_smiles(smiles, mol);


  std::vector<std::vector<Index> > cycles = relevant_cycles(mol);
  std::map<unsigned int, unsigned int> cycleSizeCounts;
  for (std::size_t i = 0; i < cycles.size(); ++i)
    cycleSizeCounts[cycles[i].size()]++;

  for (std::size_t i = 0; i < expected.size(); ++i)
    COMPARE(expected[i].second, cycleSizeCounts[expected[i].first]);
}

int main()
{
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

}
