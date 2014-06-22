// examples/molecule1.cpp
#include <Helium/smartmol.h> // for HeMol and hemol_from_smiles()

#include <iostream>

using namespace Helium;

template<typename MoleculeType>
void print_stuff(const MoleculeType &mol)
{
  // num_atoms() and num_bonds()
  std::cout << "# atoms: " << num_atoms(mol) << std::endl;
  std::cout << "# bonds: " << num_bonds(mol) << std::endl;

  // get_atom() and get_bond()
  std::cout << "degree of 2nd atom: " << get_degree(mol, get_atom(mol, 1)) << std::endl;
  std::cout << "order of 1st bond: " << get_order(mol, get_bond(mol, 0)) << std::endl;
  std::cout << "the bond C=O has order: " << get_order(mol, get_bond(mol, get_atom(mol, 1), get_atom(mol, 2))) << std::endl;

  // iterate over the atoms
  FOREACH_ATOM (atom, mol) {
    std::cout << "atom " << get_index(mol, *atom) << " has element " << get_element(mol, *atom) << std::endl;
  }

  // iterate over the bonds
  FOREACH_BOND (bond, mol) {
    std::cout << "bond " << get_index(mol, *bond) << " has order " << get_order(mol, *bond) << std::endl;
  }
}

int main()
{
  SmartMol mol;

  /*
  for (auto atom : mol.atoms()) {
  }
  */
  //print_stuff(mol);
}
