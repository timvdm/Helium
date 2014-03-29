#include <Helium/algorithms/smartcomponents.h>
#include <Helium/smartmol.h>
#include <Helium/smiles.h>
#include <Helium/fileio/moleculefile.h>

#include "benchmark.h"

using namespace Helium;

class BenchmarkAtomComponents : public SmartAttribute
{
  public:
    BenchmarkAtomComponents(const SmartMol &mol) : SmartAttribute(mol)
    {
    }

    std::string name() const
    {
      return "Helium::benchmark_connected_atom_components";
    }

    void addAtom(Index index)
    {
      m_atomComponents = connected_atom_components(molecule());
    }

    void addBond(Index index)
    {
      m_atomComponents = connected_atom_components(molecule());
    }

    void removeAtom(Index index)
    {
      m_atomComponents = connected_atom_components(molecule());
    }

    void removeBond(Index index)
    {
      m_atomComponents = connected_atom_components(molecule());
    }

    void clear()
    {
      m_atomComponents = connected_atom_components(molecule());
    }

  private:
    std::vector<unsigned int> m_atomComponents;
};


template<typename Attr>
void benchmark_connected_atom_components()
{
  MoleculeFile file(datadir() + "100K.hel");

  for (unsigned int i = 0; i < file.numMolecules(); ++i) {
    SmartMol mol;
    SmartAttribute *attr = new Attr(mol);
    mol.addAttribute(attr);
    file.readMolecule(mol);

    while (num_atoms(mol))
      remove_atom(mol, get_atom(mol, 0));
  }
}


int main(int argc, char **argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <n>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "    1     connected_atom_components (called for each add_[atom/bond] and remove_[atom/bond]" << std::endl;
    std::cerr << "    2     DynamicAtomComponents" << std::endl;
    std::cerr << std::endl;
    return -1;
  }

  int n = atoi(argv[1]);

  switch (n) {
    case 1:
      benchmark_connected_atom_components<BenchmarkAtomComponents>();
      break;
    case 2:
      benchmark_connected_atom_components<DynamicAtomComponents>();
      break;
    default:
      std::cerr << "Invalid benchmark: " << n << std::endl;
      break;
  }

}

