#include <Helium/algorithms/dijkstra.h>
#include <Helium/hemol.h>
#include <Helium/smiles.h>
#include <Helium/fileio/moleculefile.h>

#include "test.h"

using namespace Helium;

void test_all_atoms()
{
  HeMol mol;
  try {
    parse_smiles("C1CCCC2C1CCC2", mol);
  } catch(Smiley::Exception &e) {
    std::cerr << e.what();
  }

  Dijkstra<HeMol> d(mol, get_atom(mol, 0));

  COMPARE(0, d.distance(get_atom(mol, 0)));
  COMPARE(1, d.distance(get_atom(mol, 1)));
  COMPARE(2, d.distance(get_atom(mol, 2)));
  COMPARE(3, d.distance(get_atom(mol, 3)));
  COMPARE(2, d.distance(get_atom(mol, 4)));
  COMPARE(1, d.distance(get_atom(mol, 5)));
  COMPARE(2, d.distance(get_atom(mol, 6)));
  COMPARE(3, d.distance(get_atom(mol, 7)));
  COMPARE(3, d.distance(get_atom(mol, 8)));

  std::vector<HeAtom> path;

  path = d.path(get_atom(mol, 0));
  COMPARE(1, path.size());
  COMPARE(0, get_index(mol, path[0]));

  path = d.path(get_atom(mol, 1));
  COMPARE(2, path.size());
  COMPARE(0, get_index(mol, path[0]));
  COMPARE(1, get_index(mol, path[1]));

  path = d.path(get_atom(mol, 2));
  COMPARE(3, path.size());
  COMPARE(0, get_index(mol, path[0]));
  COMPARE(1, get_index(mol, path[1]));
  COMPARE(2, get_index(mol, path[2]));

  path = d.path(get_atom(mol, 3));
  COMPARE(4, path.size());
  COMPARE(0, get_index(mol, path[0]));
  COMPARE(1, get_index(mol, path[1]));
  COMPARE(2, get_index(mol, path[2]));
  COMPARE(3, get_index(mol, path[3]));

  path = d.path(get_atom(mol, 4));
  COMPARE(3, path.size());
  COMPARE(0, get_index(mol, path[0]));
  COMPARE(5, get_index(mol, path[1]));
  COMPARE(4, get_index(mol, path[2]));

  path = d.path(get_atom(mol, 5));
  COMPARE(2, path.size());
  COMPARE(0, get_index(mol, path[0]));
  COMPARE(5, get_index(mol, path[1]));

  path = d.path(get_atom(mol, 6));
  COMPARE(3, path.size());
  COMPARE(0, get_index(mol, path[0]));
  COMPARE(5, get_index(mol, path[1]));
  COMPARE(6, get_index(mol, path[2]));

  path = d.path(get_atom(mol, 7));
  COMPARE(4, path.size());
  COMPARE(0, get_index(mol, path[0]));
  COMPARE(5, get_index(mol, path[1]));
  COMPARE(6, get_index(mol, path[2]));
  COMPARE(7, get_index(mol, path[3]));

  path = d.path(get_atom(mol, 8));
  COMPARE(4, path.size());
  COMPARE(0, get_index(mol, path[0]));
  COMPARE(5, get_index(mol, path[1]));
  COMPARE(4, get_index(mol, path[2]));
  COMPARE(8, get_index(mol, path[3]));
}

struct DijkstraExcludeIndex4Atoms
{
  template<typename MoleculeType, typename AtomType>
  bool operator()(const MoleculeType &mol, AtomType atom) const
  {
    return get_index(mol, atom) != 4;
  }
};


int main()
{
  test_all_atoms();
}
