#include "../src/fileio/sdf.h"
#include "../src/smiles.h"

#include "test.h"

using namespace Helium;

// differences:
// - atom order
// [HH] -> [H][H]
// [2HH] -> [2H][H]
// [3HH] -> [3H][H]
void test_mdlbench()
{
  SDFInputFile file(datadir() + "mdlbench.sdf");
  if (file.error()) {
    std::cerr << file.error().what() << std::endl;
    return;
  }

  std::string title;
  Smiles SMILES;
  HeMol mol;

  std::stringstream output;

  while (file.isGood()) {
    if (!file.read(mol, title)) {
      std::cerr << file.error().what() << std::endl;
      continue;
    }

    output << SMILES.write(mol) << " " << title << std::endl; 
    //break;
  }

  compare_file(datadir() + "mdlbench_helium.smi", output.str());

}

int main()
{
  test_mdlbench();
}
