#include <Helium/algorithms/enumeratesubgraphs.h>
#include <Helium/smiles.h>

#include "test.h"

using namespace Helium;

std::ostream& operator<<(std::ostream &os, const Subgraph &subgraph)
{
  os << "        Subgraph:" << std::endl;
  os << "            atoms: ";
  for (std::size_t i = 0; i < subgraph.atoms.size(); ++i)
    os << subgraph.atoms[i] << " ";
  os << std::endl;
  os << "            bonds: ";
  for (std::size_t i = 0; i < subgraph.bonds.size(); ++i)
    os << subgraph.bonds[i] << " ";
  return os;
}

struct EnumerateSubgraphsCallback
{
  EnumerateSubgraphsCallback() : count(0)
  {
  }

  void operator()(const Subgraph &subgraph)
  {
    ++count;
    subgraphs.insert(std::make_pair(subgraph.atoms, subgraph.bonds));
    std::cout << subgraph << std::endl;
  }

  int count;
  std::set<std::pair<std::vector<bool>, std::vector<bool> > > subgraphs;
};

void test_impl_all_combinations()
{
  std::vector<char> v;
  v.push_back('A');
  v.push_back('B');
  v.push_back('C');

  std::vector<std::vector<char> > result;
  impl::_all_combinations(v, v.size() - 1, 0, result);
  std::cout << result << std::endl;

  COMPARE(8, result.size());

}

void testEnumerateSubgraphs(const std::string &smiles, int size = 7, bool trees = false)
{
  std::cout << "Testing: " << smiles << std::endl;
  HeMol mol;
  SMILES.read(smiles, mol);

  EnumerateSubgraphsCallback callback_correct, callback_slow, callback_fast;

  std::cout << "    Correct algorithm:" << std::endl;
  enumerate_subgraphs_correct(mol, callback_correct, size, trees);
  //std::cout << "    Slow algorithm:" << std::endl;
  //enumerate_subgraphs_slow(mol, callback_slow, size);
  std::cout << "    Fast algorithm:" << std::endl;
  enumerate_subgraphs(mol, callback_fast, size, trees);

  std::cout << mol;

  /*
  std::set<std::pair<std::vector<bool>, std::vector<bool> > >::iterator i;
  for (i = callback_correct.subgraphs.begin(); i != callback_correct.subgraphs.end(); ++i)
    if (callback_slow.subgraphs.find(*i) == callback_slow.subgraphs.end())
      std::cout << "NOT FOUND " << *i << std::endl;
  */

  COMPARE(callback_correct.count, callback_fast.count);

  // make sure there are no duplicates...
  COMPARE(callback_correct.subgraphs.size(), callback_correct.count);
  COMPARE(callback_fast.subgraphs.size(), callback_fast.count);
}

int main()
{
  test_impl_all_combinations();

  // subgraphs
  testEnumerateSubgraphs("S1C=N1");
  testEnumerateSubgraphs("C(=O)CC");
  testEnumerateSubgraphs("CCC12C(C2C)C1");
  testEnumerateSubgraphs("C2C1C(C2)C1C");
  testEnumerateSubgraphs("C2C1C(C2)C1C");
  testEnumerateSubgraphs("NCc1ccc2c(cnn2C)c1");

  testEnumerateSubgraphs("S1C=N1", 7, true);
  testEnumerateSubgraphs("C(=O)CC", 7, true);
  testEnumerateSubgraphs("CCC12C(C2C)C1", 7, true);
  testEnumerateSubgraphs("NCc1ccc2c(cnn2C)c1", 7, true);

  testEnumerateSubgraphs("Clc1ccccc1", 7, true);

  testEnumerateSubgraphs("ClC1CC1", 7, true);
}
