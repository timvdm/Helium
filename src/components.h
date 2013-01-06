#ifndef HELIUM_COMPONENTS_H
#define HELIUM_COMPONENTS_H

#include "tie.h"

namespace Helium {

  namespace impl {
    
    template<typename MoleculeType>
    void connected_bond_components(MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom,
        unsigned int number, std::vector<unsigned int> &components)
    {
      typedef typename molecule_traits<MoleculeType>::atom_bond_iter atom_bond_iter;

      atom_bond_iter bond, end_bonds;
      tie(bond, end_bonds) = get_bonds(mol, atom);
      for (; bond != end_bonds; ++bond) {
        if (components[get_index(mol, *bond)] != molecule_traits<MoleculeType>::null_index())
          continue;
        components[get_index(mol, *bond)] = number;
        connected_bond_components(mol, get_other(mol, *bond, atom), number, components);      
      }
    }
   
  }

  /**
   * Get a std::vector containing the component number for each bond. This
   * vector is indexed by bond index starting from 0. The components are
   * sequentially numbered starting from 0.
   */
  template<typename MoleculeType>
  std::vector<unsigned int> connected_bond_components(MoleculeType *mol)
  {
    typedef typename molecule_traits<MoleculeType>::mol_bond_iter mol_bond_iter;

    unsigned int number = 0;
    std::vector<unsigned int> components(num_bonds(mol), molecule_traits<MoleculeType>::null_index());

    mol_bond_iter bond, end_bonds;
    tie(bond, end_bonds) = get_bonds(mol);
    for (; bond != end_bonds; ++bond) {
      if (components[get_index(mol, *bond)] != molecule_traits<MoleculeType>::null_index())
        continue;
      impl::connected_bond_components(mol, get_source(mol, *bond), number++, components);
    }

    return components;
  }

}

#endif
