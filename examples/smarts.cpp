// examples/dfs1.cpp
#include <Helium/smarts.h>
#include <Helium/hemol.h>
#include <Helium/algorithms/cycles.h>


using namespace Helium;

int main()
{
  std::string SMARTS = "C(O)=O";
  std::string SMILES = "CCCCC(=O)[O-]";

  HeMol mol = hemol_from_smiles(SMILES);

  Smarts smarts;
  if (!smarts.init(SMARTS)) {
    std::cerr << "Error: " << smarts.error().what() << std::endl;
    return -1;
  }

  bool match;
  SingleMapping mapping;

  // only compute rings when reqlly needed!
  if (smarts.requiresCycles())
    match = smarts.search(mol, mapping, relevant_cycles(mol));
  else
    match = smarts.search(mol, mapping, RingSet<HeMol>(mol));

  if (match) {
    std::cout << SMARTS << "\t->\t" << SMILES << std::endl;
    for (std::size_t i = 0; i < mapping.map.size(); ++i)
      std::cout << i << "\t->\t" << mapping.map[i] << std::endl;
  } else {
    std::cout << "SMARTS not found in molecule" << std::endl;
  }
}
