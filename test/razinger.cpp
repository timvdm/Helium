#include <Helium/stereo/razinger.h>
#include <Helium/hemol.h>
#include <Helium/smiles.h>

#include "test.h"

using namespace Helium;

int main()
{
  HeMol mol;

  RingSet<HeMol> rings(mol);
  std::vector<unsigned long> ec;
  std::vector<bool> aromaticBonds;
  find_stereogenic_units(mol, rings, aromaticBonds, ec);
}
