#ifndef HELIUM_EXTENDEDCONNECTIVITIES_H
#define HELIUM_EXTENDEDCONNECTIVITIES_H

#include "invariants.h"

namespace Helium {

  namespace impl {

    template<typename MoleculeType>
    void extended_connectivities_iterate(MoleculeType *mol, std::vector<unsigned long> &ec)
    {
      typedef typename molecule_traits<MoleculeType>::mol_atom_iter mol_atom_iter;
      typedef typename molecule_traits<MoleculeType>::atom_atom_iter nbr_iter;

      std::vector<unsigned long> next(ec.size());
      mol_atom_iter atom, end_atoms;
      tie(atom, end_atoms) = get_atoms(mol);
      for (; atom != end_atoms; ++atom) {
        nbr_iter nbr, end_nbrs;
        tie(nbr, end_nbrs) = get_nbrs(mol, *atom);
        for (; nbr != end_nbrs; ++nbr)
          next[get_index(mol, *atom)] += ec[get_index(mol, *nbr)];
      }

      ec.swap(next);
    }

    void extended_connectivities_renumber(std::vector<unsigned long> &ec)
    {
      std::set<unsigned long> classes;
      for (std::size_t i = 0; i < ec.size(); ++i)
        classes.insert(ec[i]);
      unsigned long cls = 0;
      for (std::set<unsigned long>::iterator i = classes.begin(); i != classes.end(); ++i) {
        for (std::size_t j = 0; j < ec.size(); ++j)
          if (ec[j] == *i)
            ec[j] = cls;
        ++cls;
      }
    }

  }

  template<typename MoleculeType>
  std::vector<unsigned long> extended_connectivities(MoleculeType *mol)
  {
    typedef typename molecule_traits<MoleculeType>::mol_atom_iter mol_atom_iter;

    // initial atom invariants
    std::vector<unsigned long> ec;
    mol_atom_iter atom, end_atoms;
    tie(atom, end_atoms) = get_atoms(mol);
    for (; atom != end_atoms; ++atom)
      ec.push_back(atom_invariant(mol, *atom));
      

    unsigned int numClasses = unique_elements(ec);
    for (int i = 0; i < 100; ++i) { // should never reach 100...
      impl::extended_connectivities_iterate(mol, ec);
      unsigned int nextNumClasses = unique_elements(ec);
      if (numClasses == nextNumClasses)
        break;
      numClasses = nextNumClasses;      
    }

    impl::extended_connectivities_renumber(ec);
    return ec;
  }

}

#endif
