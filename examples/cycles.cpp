// examples/cycles.cpp
#include <Helium/algorithms/cycles.h>
#include <Helium/hemol.h>

using namespace Helium;

int main()
{
  HeMol mol = hemol_from_smiles("c1cccc2c1cc[nH]2CCCC3CC3");

  std::cout << "cyclomatic number: " << cyclomatic_number(mol) << std::endl;

  // determine ring membership
  std::vector<bool> cyclicAtoms, cyclicBonds;
  cycle_membership(mol, cyclicAtoms, cyclicBonds);

  std::cout << "cyclic atoms: ";
  for (std::size_t i = 0; i < cyclicAtoms.size(); ++i)
    if (cyclicAtoms[i])
      std::cout << i << " ";
  std::cout << std::endl;

  std::cout << "cyclic bonds: ";
  for (std::size_t i = 0; i < cyclicBonds.size(); ++i)
    if (cyclicBonds[i])
      std::cout << i << " ";
  std::cout << std::endl;

  // find cycles
  RingSet<HeMol> rings = relevant_cycles(mol);

  std::cout << "# cycles: " << rings.size() << std::endl;
  for (std::size_t i = 0; i < rings.size(); ++i) {
    std::cout << "rings " << i << ": ";
    for (std::size_t j = 0; j < rings.ring(i).size(); ++j) {
      std::cout << get_index(mol, rings.ring(i).atom(j)) << " ";
    }
    std::cout << std::endl;
  }
}
