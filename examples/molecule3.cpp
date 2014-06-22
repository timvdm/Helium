// examples/molecule3.cpp
#include <Helium/molecule.h>
#include <Helium/hemol.h>
#include <Helium/smiles.h>

#include <iostream>

using namespace Helium;

template<typename MoleculeType>
void print_stuff(const MoleculeType &mol)
{
  // iterate over the bonds
  FOREACH_BOND (bond, mol) {
    std::cout << "bond " << get_index(mol, *bond) << ":" << std::endl;
    std::cout << "    source atom index: " << get_index(mol, get_source(mol, *bond)) << std::endl;
    std::cout << "    target atom index: " << get_index(mol, get_target(mol, *bond)) << std::endl;
    std::cout << "    source atom index retrieved using get_other: "
              << get_index(mol, get_other(mol, *bond, get_target(mol, *bond))) << std::endl;
    std::cout << "    order: " << get_order(mol, *bond) << std::endl;
  }
}

int main()
{
  HeMol mol;
  Smiles SMILES;
  if (!SMILES.read("CC(=O)C", mol)) {
    std::cerr << SMILES.error().what();
    return -1;
  }

  print_stuff(mol);
}
