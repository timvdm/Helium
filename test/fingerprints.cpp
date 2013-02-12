#include "../src/fingerprints.h"
#include "../src/smiles.h"

#include "test.h"

using namespace Helium;


void test_path_fingerprint(const std::string &substructure, const std::string &superstructure)
{
  std::cout << "Testing (path): " << substructure << " < " << superstructure << std::endl;
  HeMol sub, super;
  parse_smiles(substructure, sub);
  parse_smiles(superstructure, super);

  Word subFp[16], superFp[16];
  path_fingerprint(sub, subFp);
  path_fingerprint(super, superFp);

  //print(subFp, 16);
  //print(superFp, 16);

  ASSERT(bitvec_is_subset_superset(subFp, superFp, 16));
}

void test_tree_fingerprint(const std::string &substructure, const std::string &superstructure)
{
  std::cout << "Testing (tree): " << substructure << " < " << superstructure << std::endl;
  HeMol sub, super;
  parse_smiles(substructure, sub);
  parse_smiles(superstructure, super);

  Word subFp[16], superFp[16];
  std::cout << "    generating query fingerprint..." << std::endl;
  tree_fingerprint(sub, subFp);
  std::cout << "    generating queried fingerprint..." << std::endl;
  tree_fingerprint(super, superFp);

  //print(subFp, 16);
  //print(superFp, 16);
  
  for (unsigned int i = 0; i < 1024; ++i)
    if (bitvec_get(i, subFp) && !bitvec_get(i, superFp))
      std::cout << "bit " << i << " is not in queried fingerprint" << std::endl;

  ASSERT(bitvec_is_subset_superset(subFp, superFp, 16));
}

void test_fingerprint(const std::string &substructure, const std::string &superstructure)
{
  test_path_fingerprint(substructure, superstructure);
  test_tree_fingerprint(substructure, superstructure);
}


int main()
{
  test_fingerprint("C", "CC");
  test_fingerprint("CC", "CC(C)C");
  test_fingerprint("CCC", "CC(C)C");
  test_fingerprint("CC(C)C", "CC(C)C");
  test_fingerprint("C", "ClC(Cl)Cl");
  test_fingerprint("CCl", "ClC(Cl)Cl");
  test_fingerprint("ClCCl", "ClC(Cl)Cl");
  
  test_fingerprint("CCC", "C1CCC1");
  test_fingerprint("CCC", "c1ccccc1");
  test_fingerprint("C1CCCCC1", "c1ccccc1");
  test_fingerprint("CCCOCCl", "c1ccccc1COCCl");
  
  
  test_fingerprint("Clc1ccccc1", "CCc1cc(Cl)c(C=O)cc1");

}
