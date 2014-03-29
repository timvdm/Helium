// examples/molecule3.cpp
#include <Helium/molecule.h>
#include <Helium/hemol.h> // for HeMol and hemol_from_smiles()

#include <iostream>

using namespace Helium;

template<typename MoleculeType>
void print_stuff(const MoleculeType &mol)
{
  // iterate over the bonds
  FOREACH_BOND_T (bond, mol, MoleculeType) {
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
  HeMol mol = hemol_from_smiles("CC(=O)C");
  print_stuff(mol);
}
