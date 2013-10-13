/*
 * Copyright (c) 2013, Tim Vandermeersch
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef HELIUM_DFS_H
#define HELIUM_DFS_H

#include <iostream>

#include <Helium/molecule.h>

namespace Helium {

  /**
   * @brief Base class for DFS functors.
   */
  template<typename MoleculeType>
  struct DFSVisitor
  {
    /**
     * @brief The molecule type.
     */
    typedef MoleculeType molecule_type;
    /**
     * @brief The atom type.
     */
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    /**
     * @brief The bond type.
     */
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

    /**
     * @brief Initialize the functor.
     *
     * This function is called once when the DFS search is started.
     *
     * @param mol The molecule.
     */
    void initialize(const MoleculeType &mol) {}

    /**
     * @brief Visit an atom.
     *
     * This function is called every time a new atom is traversed.
     *
     * @param mol The molecule.
     * @param prev The previous atom on the DFS path.
     * @param atom The traversed atom.
     */
    void atom(const MoleculeType &mol, atom_type prev, atom_type atom) {}

    /**
     * @brief Visit a bond.
     *
     * This function is called every time a new bond is traversed.
     *
     * @param mol The molecule.
     * @param prev The previous atom on the DFS path.
     * @param bond The traversed bond.
     */
    void bond(const MoleculeType &mol, atom_type prev, bond_type bond) {}

    /**
     * @brief Backtrack.
     *
     * This function is called when the current path is exhausted and the DFS
     * search needs to backtrack.
     *
     * @param mol The molecule.
     * @param atom The atom.
     */
    void backtrack(const MoleculeType &mol, atom_type atom) {}

    /**
     * @brief Visit a "back" bond.
     *
     * The "back" bonds are the chords (or ring closure bonds) that would make
     * the DFS spanning graph cyclic.
     *
     * @param mol The molecule.
     * @param bond The "back" bond.
     */
    void back_bond(const MoleculeType &mol, bond_type bond) {}
  };

  namespace impl {

    template<typename MoleculeType, typename AtomType, typename DFSVisitorType>
    void dfs_visit(const MoleculeType &mol, AtomType atom, DFSVisitorType &visitor, std::vector<bool> &visited,
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

    template<typename MoleculeType, typename AtomType, typename DFSVisitorType>
    void exhaustive_dfs_visit(const MoleculeType &mol, AtomType atom, DFSVisitorType &visitor, std::vector<bool> &visited,
        AtomType prev = molecule_traits<MoleculeType>::null_atom())
    {
      typedef AtomType atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
      typedef typename molecule_traits<MoleculeType>::incident_iter incident_iter;

      // mark atom as visited
      visited[get_index(mol, atom)] = true;
      // invoke atom visitor
      visitor.atom(mol, prev, atom);

      // create list of bonds to visit next
      std::vector<bond_type> bonds;
      FOREACH_INCIDENT (bond, atom, mol, MoleculeType) {
        if (!visited[num_atoms(mol) + get_index(mol, *bond)])
          bonds.push_back(*bond);
      }

      // sort the list of bonds
      std::sort(bonds.begin(), bonds.end());

      std::vector<bool> oldVisited = visited;
      do {
        if (bonds.size() > 1)
        visited = oldVisited;

        // call dfs_visit for all unvisited neighbors of v
        for (typename std::vector<bond_type>::iterator bond = bonds.begin(); bond != bonds.end(); ++bond) {
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

          exhaustive_dfs_visit(mol, nbr, visitor, visited, atom);
        }

      } while (std::next_permutation(bonds.begin(), bonds.end()));

      // invoke backtrack visitor
      visitor.backtrack(mol, atom);
    }

  }

  /**
   * @brief Perform a depth-first search (DFS).
   *
   * This depth-first search function considers all fragments and walks a single
   * spanning tree for each fragment. The atom with the lowest index from each
   * fragment is used as root of the resulting spanning tree to start the
   * search.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * @note Complexity: O(n)
   *
   * @param mol The molecule.
   * @param visitor The DFS visitor functor.
   */
  template<typename MoleculeType, typename DFSVisitorType>
  void depth_first_search(MoleculeType &mol, DFSVisitorType &visitor)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;

    visitor.initialize(mol);

    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    FOREACH_ATOM (atom, mol, MoleculeType) {
      if (!visited[get_index(mol, *atom)])
        impl::dfs_visit(mol, *atom, visitor, visited);
    }
  }

  /**
   * @brief Perform an exhaustive depth-first search (DFS).
   *
   * This depth-first search function considers a single fragment and starts
   * walking from the specified atom. If there are multiple neighbors to visit
   * next, all permutations are tried making this an exhaustive search.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * @note Complexity: O(2^n)
   *
   * @param mol The molecule.
   * @param atom The root atom for the DFS spanning tree.
   * @param visitor The DFS visitor functor.
   */
  template<typename MoleculeType, typename AtomType, typename DFSVisitorType>
  void exhaustive_depth_first_search(MoleculeType &mol, AtomType atom, DFSVisitorType &visitor)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;

    visitor.initialize(mol);

    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    impl::exhaustive_dfs_visit(mol, atom, visitor, visited);
  }

  /**
   * @brief A DFS visitor that records the order in which the atoms are visited.
   */
  template<typename MoleculeType>
  struct DFSAtomOrderVisitor : public DFSVisitor<MoleculeType>
  {
    /**
     * @brief The atom type.
     */
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

    /**
     * @brief Visit an atom.
     *
     * This function records the atom index.
     *
     * @param mol The molecule.
     * @param prev The previous atom on the DFS path.
     * @param atom The traversed atom.
     */
    void atom(const MoleculeType &mol, atom_type prev, atom_type atom)
    {
      atoms.push_back(get_index(mol, atom));
    }

    /**
     * @brief The list of atom indices.
     */
    std::vector<Index> atoms;
  };

  /**
   * @brief A DFS visitor that records the order in which the bonds are visited.
   */
  template<typename MoleculeType>
  struct DFSBondOrderVisitor : public DFSVisitor<MoleculeType>
  {
    /**
     * @brief The atom type.
     */
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    /**
     * @brief The bond type.
     */
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

    /**
     * @brief Visit a bond.
     *
     * This function records the bond index.
     *
     * @param mol The molecule.
     * @param prev The previous atom on the DFS path.
     * @param bond The traversed bond.
     */
    void bond(const MoleculeType &mol, atom_type prev, bond_type bond)
    {
      bonds.push_back(get_index(mol, bond));
    }

    /**
     * @brief The list of bond indices.
     */
    std::vector<Index> bonds;
  };

  /**
   * @brief A DFS visitor that records the "back" bonds.
   */
  template<typename MoleculeType>
  struct DFSClosureRecorderVisitor : public DFSVisitor<MoleculeType>
  {
    /**
     * @brief The bond type.
     */
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

    /**
     * @brief Visit a "back" bond.
     *
     * This function records the index of the "back" (or closure) bond.
     *
     * @param mol The molecule.
     * @param bond The "back" bond.
     */
    void back_bond(const MoleculeType &mol, bond_type bond)
    {
      back_bonds.push_back(get_index(mol, bond));
    }

    /**
     * @brief The list of back bonds.
     */
    std::vector<Index> back_bonds;
  };

  /**
   * @brief A DFS visitor that prints out debug information.
   */
  template<typename MoleculeType>
  struct DFSDebugVisitor : public DFSVisitor<MoleculeType>
  {
    /**
     * @brief The atom type.
     */
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    /**
     * @brief The bond type.
     */
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

    /**
     * @brief Constructor.
     *
     * @param os The STL output stream to print to.
     */
    DFSDebugVisitor(std::ostream &os_ = std::cout) : os(os_)
    {
    }

    /**
     * @brief Visit an atom.
     *
     * This function prints the string "atom(i)" to os where i is the
     * atom index.
     *
     * @param mol The molecule.
     * @param prev The previous atom on the DFS path.
     * @param atom The traversed atom.
     */
    void atom(const MoleculeType &mol, atom_type prev, atom_type atom)
    {
      os << "atom(" << get_index(mol, atom) << ")" << std::endl;
    }

    /**
     * @brief Visit a bond.
     *
     * This function prints the string "bond(i)" to os where i is the
     * bond index.
     *
     * @param mol The molecule.
     * @param prev The previous atom on the DFS path.
     * @param bond The traversed bond.
     */
    void bond(const MoleculeType &mol, atom_type prev, bond_type bond)
    {
      os << "bond(" << get_index(mol, bond) << ")" << std::endl;
    }

    /**
     * @brief Backtrack.
     *
     * This function prints the string "backtrack(i)" to os where i is the
     * atom index.
     *
     * @param mol The molecule.
     * @param atom The atom.
     */
    void backtrack(const MoleculeType &mol, atom_type atom)
    {
      os << "backtrack(" << get_index(mol, atom) << ")" << std::endl;
    }

    /**
     * @brief Visit a "back" bond.
     *
     * This function prints the string "back_bond(i)" to os where i is the
     * "back" bond index.
     *
     * @param mol The molecule.
     * @param bond The "back" bond.
     */
    void back_bond(const MoleculeType &mol, bond_type bond)
    {
      os << "back_bond(" << get_index(mol, bond) << ")" << std::endl;
    }

    std::ostream &os; //!< The STL output stream.
  };

}

#endif
