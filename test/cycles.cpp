#include "../src/cycles.h"
#include "../src/fileio.h"

#include "test.h"

using namespace Helium;

void test_cyclomatic_number(const std::string &smiles, unsigned int expected)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol;
  read_smiles(smiles, mol);
  COMPARE(expected, cyclomatic_number(&mol));
}

void test_relevant_cycles(const std::string &smiles, std::vector<std::pair<unsigned int, unsigned int> > &expected)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol;
  read_smiles(smiles, mol);


  std::vector<std::vector<Index> > cycles = relevant_cycles(&mol);
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
