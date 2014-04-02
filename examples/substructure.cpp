// examples/substructure.cpp
#include <Helium/substructure.h>
#include <Helium/algorithms/cycles.h>
#include <Helium/hemol.h>
#include <Helium/smiles.h>
#include <Helium/smarts.h>

using namespace Helium;

template<typename MoleculeType>
void print_stuff(const MoleculeType &mol)
{
  // see examples/smarts.cpp
  Smarts smarts;
  if (!smarts.init("c1ccccc1")) {
    std::cerr << "Error: " << smarts.error().what() << std::endl;
    return;
  }

  // perform a SMARTS search
  SingleMapping mapping;
  if (!smarts.findMapping(mol, relevant_cycles(mol), mapping)) {
    std::cout << "SMARTS did not match SMILES" << std::endl;
    return;
  }

  // create the substructure atom mask
  std::vector<bool> atoms(num_atoms(mol));
  for (std::size_t i = 0; i < mapping.map.size(); ++i)
    atoms[mapping.map[i]] = true;

  std::vector<bool> bonds(num_atoms(mol));
  FOREACH_BOND_T (bond, mol, MoleculeType)
    if (atoms[get_index(mol, get_source(mol, *bond))] &&
        atoms[get_index(mol, get_target(mol, *bond))])
      bonds[get_index(mol, *bond)] = true;

  // create the substructure
  Substructure<HeMol> substruct(mol, atoms, bonds);

  std::cout << "# atoms: " << num_atoms(substruct) << std::endl;
  std::cout << "# bonds: " << num_bonds(substruct) << std::endl;

  // iterate over substruct's atoms
  FOREACH_ATOM_T (atom, substruct, Substructure<MoleculeType>) {
    std::cout << "atom " << get_index(substruct, *atom) << " has element " << get_element(substruct, *atom) << std::endl;
  }

  // iterate over substruct's bonds
  FOREACH_BOND_T (bond, substruct, Substructure<MoleculeType>) {
    std::cout << "bond " << get_index(substruct, get_source(substruct, *bond)) << "-" <<
      get_index(substruct, get_target(substruct, *bond)) << " has order " << get_order(substruct, *bond) << std::endl;
  }
}

int main()
{
  // read a SMILES string
  HeMol mol;
  Smiles SMILES;
  if (!SMILES.read("c1ccccc1C(=O)[O-]", mol)) {
    std::cerr << SMILES.error().what();
    return -1;
  }

  print_stuff(mol);
}
