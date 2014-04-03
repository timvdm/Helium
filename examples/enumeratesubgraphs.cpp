#include <Helium/algorithms/enumeratesubgraphs.h>
#include <Helium/algorithms/extendedconnectivities.h>
#include <Helium/algorithms/canonical.h>
#include <Helium/smiles.h>

using namespace Helium;

Smiles SMILES;

template<typename MoleculeType>
struct SubgraphsCallback
{
  SubgraphsCallback(const MoleculeType &mol_)
    : mol(mol_)
  {
  }

  void operator()(const Subgraph &subgraph)
  {
    // create the subgraph molecule
    Substructure<MoleculeType> substruct(mol, subgraph.atoms, subgraph.bonds);
    // compute symmetry classes
    std::vector<unsigned long> symmetry = extended_connectivities(substruct,
        DefaultAtomInvariant(DefaultAtomInvariant::Element));
    // canonicalize the subgraph
    std::vector<Index> canon = canonicalize_component(substruct, symmetry,
        DefaultAtomInvariant(DefaultAtomInvariant::Element),
        DefaultBondInvariant(DefaultBondInvariant::Order)).first;
    // write subgraph SMILES
    std::string smiles = SMILES.write(substruct, canon);
    // add to features
    std::map<std::string, int>::iterator feature = features.find(smiles);
    if (feature == features.end())
      features[smiles] = 1;
    else
      feature->second++;
  }

  const MoleculeType &mol; //!< The molecule
  std::map<std::string, int> features;
};


int main(int argc, char **argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <SMILES>" << std::endl;
    return -1;
  }

  HeMol mol;
  Smiles SMILES;
  if (!SMILES.read(argv[1], mol)) {
    std::cerr << SMILES.error().what();
    return -1;
  }

  SubgraphsCallback<HeMol> callback(mol);
  enumerate_subgraphs(mol, callback, 7, false);

  for (std::map<std::string, int>::iterator i = callback.features.begin(); i != callback.features.end(); ++i)
    std::cout << i->first << " " << i->second << " ";
  std::cout << std::endl;
}
