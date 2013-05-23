#ifndef HELIUM_INVARIANTS_H
#define HELIUM_INVARIANTS_H

#include <Helium/tie.h>
#include <Helium/molecule.h>

namespace Helium {

  template<typename MoleculeType>
  unsigned int atom_invariant(MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    unsigned int invariant = 0;
    // bit 0: aromatic
    //invariant |= is_aromatic(mol, atom);
    // bit 1: cyclic
    //invariant |= is_cyclic(mol, atom) << 1;
    // bit 2-8: element
    invariant |= get_element(mol, atom) << 2;
    // bit 9-15: mass
    //invariant |= get_mass(mol, atom) << 9;
    // bit 16-19: degree
    invariant |= get_degree(mol, atom) << 16;
    // bit 20-22: hydrogens
    //invariant |= num_hydrogens(mol, atom) << 20;
    // bit 23-26: charge
    //invariant |= (get_charge(mol, atom) + 7) << 23;

    return invariant;
  }

  template<typename MoleculeType>
  unsigned int bond_invariant(MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond)
  {
    unsigned int invariant = 0;
    // bit 0: aromatic
    //invariant |= is_aromatic(mol, bond);
    // bit 1: cyclic
    //invariant |= is_cyclic(mol, bond) << 1;
    // bit 2-4: element
    invariant |= get_order(mol, bond) << 2;

    return invariant;
  }


}

#endif
