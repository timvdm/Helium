#include <Helium/fileio/molecules.h>
#include <Helium/algorithms/enumeratesubgraphs.h>
#include <Helium/algorithms/extendedconnectivities.h>
#include <Helium/algorithms/canonical.h>
#include <Helium/smiles.h>

using namespace Helium;

template<typename MoleculeType>
struct SubgraphsCallback
{
  SubgraphsCallback(MoleculeType &mol_) //FIXME make const
    : mol(mol_)
  {
  }

  void operator()(const Subgraph &subgraph)
  {
    // create the subgraph molecule
    Substructure<MoleculeType> substruct(mol, subgraph.atoms, subgraph.bonds);
    // compute symmetry classes
    std::vector<unsigned long> symmetry = extended_connectivities(substruct);
    // canonicalize the subgraph
    std::vector<Index> canon = canonicalize(substruct, symmetry, AtomElementAttribute(), BondOrderAttribute()).first;
    // write subgraph SMILES
    std::string smiles = write_smiles(substruct, canon, 0);
    // add to features
    std::map<std::string, int>::iterator feature = features.find(smiles);
    if (feature == features.end())
      features[smiles] = 1;
    else
      feature->second++;
  }

  MoleculeType &mol; //!< The molecule
  std::map<std::string, int> features;
};


int main(int argc, char **argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
    return -1;
  }

  MemoryMappedMoleculeFile moleculeFile;
  try {
    moleculeFile.load(argv[1]);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }

  HeMol mol;
  for (unsigned int i = 0; i < moleculeFile.numMolecules(); ++i) {
    if ((i % 100) == 0)
      std::cerr << i << std::endl;

    moleculeFile.read_molecule(i, mol);

    SubgraphsCallback<HeMol> callback(mol);
    enumerate_subgraphs(mol, callback, 7, false);

    for (std::map<std::string, int>::iterator i = callback.features.begin(); i != callback.features.end(); ++i)
      std::cout << i->first << " " << i->second << " ";
    std::cout << std::endl;


  }

}
