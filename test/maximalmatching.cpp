#include "../src/maximalmatching.h"
#include "../src/fileio.h"

#include "test.h"

using namespace Helium;

int main()
{
  HeMol mol;

  read_smiles("c1ccccc1", mol);

  maximal_matching(&mol);
}
