// examples/molecule2.cpp
#include <Helium/molecule.h>
#include <Helium/hemol.h> // for HeMol and hemol_from_smiles()
#include <Helium/smiles.h>

#include <iostream>

using namespace Helium;

template<typename MoleculeType>
void print_stuff(const MoleculeType &mol)
{
  // iterate over the atoms
  FOREACH_ATOM (atom, mol) {
    std::cout << "atom " << get_index(mol, *atom) << ":" << std::endl;
    std::cout << "    element: " << get_element(mol, *atom) << std::endl;
    std::cout << "    mass: " << get_mass(mol, *atom) << std::endl;
    std::cout << "    charge: " << get_charge(mol, *atom) << std::endl;
    std::cout << "    degree: " << get_degree(mol, *atom) << std::endl;
    std::cout << "    number of hydrogens: " << get_hydrogens(mol, *atom) << std::endl;
  }

  // print neighbor indices for atom 1
  std::cout << "neighbor indices for atom 1: ";
  FOREACH_NBR (nbr, get_atom(mol, 1), mol)
    std::cout << get_index(mol, *nbr) << " ";
  std::cout << std::endl;

  // print incident bond indices for atom 1
  std::cout << "incident bond indices for atom 1: ";
  FOREACH_INCIDENT (bond, get_atom(mol, 1), mol)
    std::cout << get_index(mol, *bond) << " ";
  std::cout << std::endl;
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
