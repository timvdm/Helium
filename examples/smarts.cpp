// examples/dfs1.cpp
#include <Helium/smarts.h>
#include <Helium/hemol.h>
#include <Helium/smiles.h>
#include <Helium/algorithms/cycles.h>


using namespace Helium;

int main()
{
  std::string smartsString = "C(O)=O";
  std::string smilesString = "CCCCC(=O)[O-]";

  HeMol mol;
  Smiles SMILES;
  if (!SMILES.read(smilesString, mol)) {
    std::cerr << SMILES.error().what();
    return -1;
  }

  Smarts smarts;
  if (!smarts.init(smartsString)) {
    std::cerr << "Error: " << smarts.error().what() << std::endl;
    return -1;
  }

  bool match;
  SingleMapping mapping;

  // only compute rings when reqlly needed!
  if (smarts.requiresCycles())
    match = smarts.findMapping(mol, relevant_cycles(mol), mapping);
  else
    match = smarts.findMapping(mol, RingSet<HeMol>(mol), mapping);

  if (match) {
    std::cout << smartsString << "\t->\t" << smilesString << std::endl;
    for (std::size_t i = 0; i < mapping.map.size(); ++i)
      std::cout << i << "\t->\t" << mapping.map[i] << std::endl;
  } else {
    std::cout << "SMARTS not found in molecule" << std::endl;
  }
}
