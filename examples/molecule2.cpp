// examples/molecule2.cpp
#include <Helium/molecule.h>
#include <Helium/hemol.h> // for HeMol and hemol_from_smiles()

#include <iostream>

using namespace Helium;

template<typename MoleculeType>
void print_stuff(const MoleculeType &mol)
{
  // iterate over the atoms
  FOREACH_ATOM (atom, mol, MoleculeType) {
    std::cout << "atom " << get_index(mol, *atom) << ":" << std::endl;
    std::cout << "    element: " << get_element(mol, *atom) << std::endl;
    std::cout << "    mass: " << get_mass(mol, *atom) << std::endl;
    std::cout << "    charge: " << get_charge(mol, *atom) << std::endl;
    std::cout << "    degree: " << get_degree(mol, *atom) << std::endl;
    std::cout << "    number of hydrogens: " << num_hydrogens(mol, *atom) << std::endl;
  }

  // print neighbor indices for atom 1
  std::cout << "neighbor indices for atom 1: ";
  FOREACH_NBR (nbr, get_atom(mol, 1), mol, MoleculeType)
    std::cout << get_index(mol, *nbr) << " ";
  std::cout << std::endl;

  // print incident bond indices for atom 1
  std::cout << "incident bond indices for atom 1: ";
  FOREACH_INCIDENT (bond, get_atom(mol, 1), mol, MoleculeType)
    std::cout << get_index(mol, *bond) << " ";
  std::cout << std::endl;
}

int main()
{
  HeMol mol = hemol_from_smiles("CC(=O)C");
  print_stuff(mol);
}
