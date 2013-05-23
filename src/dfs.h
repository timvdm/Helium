#ifndef HELIUM_DFS_H
#define HELIUM_DFS_H

#include <iostream>

#include <Helium/molecule.h>

namespace Helium {

  template<typename MoleculeType>
  struct DFSVisitor
  {
    typedef MoleculeType molecule_type;
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

    void initialize(MoleculeType &mol) {}
    void atom(MoleculeType &mol, atom_type prev, atom_type atom) {}
    void bond(MoleculeType &mol, atom_type prev, bond_type bond) {}
    void backtrack(MoleculeType &mol, atom_type atom) {}
    void back_bond(MoleculeType &mol, bond_type bond) {}
  };

  namespace impl {

    template<typename MoleculeType, typename AtomType, typename DFSVisitorType>
    void dfs_visit(MoleculeType &mol, AtomType atom, DFSVisitorType &visitor, std::vector<bool> &visited,
        AtomType prev = molecule_traits<MoleculeType>::null_atom())
    {
      typedef AtomType atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
      typedef typename molecule_traits<MoleculeType>::incident_iter incident_iter;

      // mark atom as visited
      visited[get_index(mol, atom)] = true;
      // invoke atom visitor
      visitor.atom(mol, prev, atom);

      // call dfs_visit for all unvisited neighbors of v
      FOREACH_INCIDENT (bond, atom, mol, MoleculeType) {
        atom_type nbr = get_other(mol, *bond, atom);

        if (visited[get_index(mol, nbr)]) {
          // if this bond has not been visited before, a back_bond has been found
          if (!visited[num_atoms(mol) + get_index(mol, *bond)])
            visitor.back_bond(mol, *bond);
          // mark bond as visited
          visited[num_atoms(mol) + get_index(mol, *bond)] = true;
          continue;
        }

        // mark bond as visited
        visited[num_atoms(mol) + get_index(mol, *bond)] = true;
        // invoke bond visitor
        visitor.bond(mol, prev, *bond);

        dfs_visit(mol, nbr, visitor, visited, atom);
      }

      // invoke backtrack visitor
      visitor.backtrack(mol, atom);
    }

  }

  template<typename MoleculeType, typename DFSVisitorType>
  void depth_first_search(MoleculeType &mol, DFSVisitorType &visitor)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;

    visitor.initialize(mol);

    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    FOREACH_ATOM (atom, mol, MoleculeType) {
      if (!visited[get_index(mol, *atom)])
        dfs_visit(mol, *atom, visitor, visited);
    }
  }


  template<typename MoleculeType>
  struct DFSAtomOrderVisitor : public DFSVisitor<MoleculeType>
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

    void atom(MoleculeType &mol, atom_type prev, atom_type atom)
    {
      atoms.push_back(get_index(mol, atom));
    }

    std::vector<Index> atoms;
  };

  template<typename MoleculeType>
  struct DFSBondOrderVisitor : public DFSVisitor<MoleculeType>
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

    void bond(MoleculeType &mol, atom_type prev, bond_type bond)
    {
      bonds.push_back(get_index(mol, bond));
    }

    std::vector<Index> bonds;
  };

  template<typename MoleculeType>
  struct DFSClosureRecorderVisitor : public DFSVisitor<MoleculeType>
  {
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

    void back_bond(MoleculeType &mol, bond_type bond)
    {
      back_bonds.push_back(get_index(mol, bond));
    }

    std::vector<Index> back_bonds;
  };

  template<typename MoleculeType>
  struct DFSDebugVisitor : public DFSVisitor<MoleculeType>
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

    void atom(MoleculeType &mol, atom_type prev, atom_type atom)
    {
      std::cout << "atom(" << get_index(mol, atom) << ")" << std::endl;
    }

    void bond(MoleculeType &mol, atom_type prev, bond_type bond)
    {
      std::cout << "bond(" << get_index(mol, bond) << ")" << std::endl;
    }

    void backtrack(MoleculeType &mol, atom_type atom)
    {
      std::cout << "backtrack(" << get_index(mol, atom) << ")" << std::endl;
    }

    void back_bond(MoleculeType &mol, bond_type bond)
    {
      std::cout << "back_bond(" << get_index(mol, bond) << ")" << std::endl;
    }
  };



}

#endif
