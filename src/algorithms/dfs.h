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
 * DIRECT, INDIRECT, INCIDENT_TAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
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
#include <Helium/util.h>

namespace Helium {

  /**
   * @file algorithms/dfs.h
   * @brief Depth-first search.
   */

  //////////////////////////////////////////////////////////////////////////////
  //
  // DFSVisitor base class
  //
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @struct DFSVisitor algorithms/dfs.h <Helium/algorithms/dfs.h>
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
     * @brief A new component starts. Not called for exhaustive_depth_first_search().
     *
     * @param i The component (i.e. [0,n])
     */
    void component(int i) {}

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

    /**
     * @brief Check if the search should be stopped.
     *
     * By overloading this function it is possible to abort the search before
     * it is complete (e.g. when some goal is achieved).
     *
     * @return True if the search should be stopped.
     */
    bool stop() const
    {
      return false;
    }
  };

  //////////////////////////////////////////////////////////////////////////////
  //
  // basic depth first search
  //
  //////////////////////////////////////////////////////////////////////////////

  namespace impl {

    // basic dfs recursive function
    template<typename MoleculeType, typename AtomType, typename DFSVisitorType>
    void dfs_visit(const MoleculeType &mol, AtomType atom, DFSVisitorType &visitor, std::vector<bool> &visited,
        AtomType prev = molecule_traits<MoleculeType>::null_atom())
    {
      typedef AtomType atom_type;

      if (visitor.stop())
        return;

      // mark atom as visited
      visited[get_index(mol, atom)] = true;
      // invoke atom visitor
      visitor.atom(mol, prev, atom);

      // call dfs_visit for all unvisited neighbors of v
      FOREACH_INCIDENT (bond, atom, mol) {
        if (visitor.stop())
          return;

        if (visited[num_atoms(mol) + get_index(mol, *bond)])
          continue;

        atom_type nbr = get_other(mol, *bond, atom);

        if (visited[get_index(mol, nbr)]) {
          // if this bond has not been visited before, a back_bond has been found
          if (!visited[num_atoms(mol) + get_index(mol, *bond)]) {
            visitor.back_bond(mol, *bond);
            // mark bond as visited
            visited[num_atoms(mol) + get_index(mol, *bond)] = true;
          }
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

  } // namespace impl

  /**
   * @brief Perform a depth-first search (DFS).
   *
   * This depth-first search function considers all components and walks a single
   * spanning tree for each component. The atom with the lowest index from each
   * component is used as root of the resulting spanning tree to start the
   * search.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * @param mol The molecule.
   * @param visitor The DFS visitor functor.
   */
  template<typename MoleculeType, typename DFSVisitorType>
  void depth_first_search(const MoleculeType &mol, DFSVisitorType &visitor)
  {
    visitor.initialize(mol);

    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    int c = 0;
    FOREACH_ATOM (atom, mol) {
      if (visitor.stop())
        return;

      if (!visited[get_index(mol, *atom)]) {
        visitor.component(c++);
        impl::dfs_visit(mol, *atom, visitor, visited);
      }
    }
  }

  /**
   * @brief Perform a depth-first search (DFS).
   *
   * This depth-first search function considers all components and walks a single
   * spanning tree for each component. The atom with the lowest index from each
   * component is used as root of the resulting spanning tree to start the
   * search.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * The @p atomMask specifies the atoms and bonds that will be considered.
   * The element in the @p atomMask that correspond with atoms to be included
   * in the search should be set to true. If an atom is not included in the
   * search, all bonds around the atom are also excluded.
   *
   * @param mol The molecule.
   * @param visitor The DFS visitor functor.
   * @param atomMask The atom mask.
   */
  template<typename MoleculeType, typename DFSVisitorType>
  void depth_first_search_mask(const MoleculeType &mol, DFSVisitorType &visitor,
      const std::vector<bool> &atomMask)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

    // initialize the visitor
    visitor.initialize(mol);

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    // initialize visited to take mask into account
    for (std::size_t i = 0; i < num_atoms(mol); ++i)
      if (!atomMask[i]) {
        // mark masked atoms as visited
        visited[i] = true;
        // mark bonds around atom as visited
        atom_type atom = get_atom(mol, i);
        FOREACH_INCIDENT (bond, atom, mol)
          visited[num_atoms(mol) + get_index(mol, *bond)] = true;
      }

    int c = 0;
    FOREACH_ATOM (atom, mol) {
      if (visitor.stop())
        return;

      if (!visited[get_index(mol, *atom)]) {
        // let the visitor know a new component (i.e. tree) starts
        visitor.component(c++);
        // initiate recursive search
        impl::dfs_visit(mol, *atom, visitor, visited);
      }
    }
  }

  /**
   * @brief Perform a depth-first search (DFS).
   *
   * This depth-first search function considers all components and walks a single
   * spanning tree for each component. The atom with the lowest index from each
   * component is used as root of the resulting spanning tree to start the
   * search.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * The @p atomMask specifies the atoms that will be considered.
   * The element in the @p atomMask that correspond with atoms to be included
   * in the search should be set to true. The @p bondMask works in the same way
   * but specifies the bonds that should be considered.
   *
   * @param mol The molecule.
   * @param visitor The DFS visitor functor.
   * @param atomMask The atom mask.
   * @param bondMask The bond mask.
   */
  template<typename MoleculeType, typename DFSVisitorType>
  void depth_first_search_mask(const MoleculeType &mol, DFSVisitorType &visitor,
      const std::vector<bool> &atomMask, const std::vector<bool> &bondMask)
  {
    // initialize the visitor
    visitor.initialize(mol);

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    // initialize visited to take masks into account
    for (std::size_t i = 0; i < num_atoms(mol); ++i)
      if (!atomMask[i])
        visited[i] = true;
    for (std::size_t i = 0; i < num_bonds(mol); ++i)
      if (!bondMask[i])
        visited[num_atoms(mol) + i] = true;

    int c = 0;
    FOREACH_ATOM (atom, mol) {
      if (visitor.stop())
        return;

      if (!visited[get_index(mol, *atom)]) {
        // let the visitor know a new component (i.e. tree) starts
        visitor.component(c++);
        // initiate recursive search
        impl::dfs_visit(mol, *atom, visitor, visited);
      }
    }
  }

  /**
   * @brief Perform a depth-first search (DFS).
   *
   * This depth-first search function considers a single component and walks a
   * single spanning tree starting at @p atom.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * @param mol The molecule.
   * @param atom The start atom (i.e. root of the spanning tree).
   * @param visitor The DFS visitor functor.
   */
  template<typename MoleculeType, typename DFSVisitorType>
  void depth_first_search(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom, DFSVisitorType &visitor)
  {
    // initialize the visitor
    visitor.initialize(mol);

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    visitor.component(0);
    impl::dfs_visit(mol, atom, visitor, visited);
  }

  /**
   * @brief Perform a depth-first search (DFS).
   *
   * This depth-first search function considers a single component and walks a
   * single spanning tree starting at @p atom.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * The @p atomMask specifies the atoms and bonds that will be considered.
   * The element in the @p atomMask that correspond with atoms to be included
   * in the search should be set to true. If an atom is not included in the
   * search, all bonds around the atom are also excluded.
   *
   * @param mol The molecule.
   * @param atom The start atom (i.e. root of the spanning tree).
   * @param visitor The DFS visitor functor.
   * @param atomMask The atom mask.
   */
  template<typename MoleculeType, typename DFSVisitorType>
  void depth_first_search_mask(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom, DFSVisitorType &visitor,
      const std::vector<bool> &atomMask)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

    // initialize the visitor
    visitor.initialize(mol);

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    // initialize visited to take mask into account
    for (std::size_t i = 0; i < num_atoms(mol); ++i)
      if (!atomMask[i]) {
        // mark masked atoms as visited
        visited[i] = true;
        // mark bonds around atom as visited
        atom_type atom = get_atom(mol, i);
        FOREACH_INCIDENT (bond, atom, mol)
          visited[num_atoms(mol) + get_index(mol, *bond)] = true;
      }

    visitor.component(0);
    impl::dfs_visit(mol, atom, visitor, visited);
  }

  /**
   * @brief Perform a depth-first search (DFS).
   *
   * This depth-first search function considers a single component and walks a
   * single spanning tree starting at @p atom.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * The @p atomMask specifies the atoms that will be considered.
   * The element in the @p atomMask that correspond with atoms to be included
   * in the search should be set to true. The @p bondMask works in the same way
   * but specifies the bonds that should be considered.
   *
   * @param mol The molecule.
   * @param atom The start atom (i.e. root of the spanning tree).
   * @param visitor The DFS visitor functor.
   * @param atomMask The atom mask.
   * @param bondMask The bond mask.
   */
  template<typename MoleculeType, typename DFSVisitorType>
  void depth_first_search_mask(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom, DFSVisitorType &visitor,
      const std::vector<bool> &atomMask, const std::vector<bool> &bondMask)
  {
    // initialize the visitor
    visitor.initialize(mol);

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    // initialize visited to take masks into account
    for (std::size_t i = 0; i < num_atoms(mol); ++i)
      if (!atomMask[i])
        visited[i] = true;
    for (std::size_t i = 0; i < num_bonds(mol); ++i)
      if (!bondMask[i])
        visited[num_atoms(mol) + i] = true;

    visitor.component(0);
    impl::dfs_visit(mol, atom, visitor, visited);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // depth first search with specified atom order
  //
  //////////////////////////////////////////////////////////////////////////////

  namespace impl {

    // dfs recursive function taking specified order into account
    template<typename MoleculeType, typename AtomType, typename DFSVisitorType>
    void ordered_dfs_visit(const MoleculeType &mol, AtomType atom, const std::vector<Index> &order,
        DFSVisitorType &visitor, std::vector<bool> &visited,
        AtomType prev = molecule_traits<MoleculeType>::null_atom())
    {
      typedef AtomType atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      if (visitor.stop())
        return;

      // mark atom as visited
      visited[get_index(mol, atom)] = true;
      // invoke atom visitor
      visitor.atom(mol, prev, atom);

      std::vector<std::pair<std::size_t, bond_type> > bonds;
      FOREACH_INCIDENT (bond, atom, mol) {
        if (visited[num_atoms(mol) + get_index(mol, *bond)])
          continue;
        atom_type nbr = get_other(mol, *bond, atom);
        bonds.push_back(std::make_pair(std::find(order.begin(), order.end(),
                get_index(mol, nbr)) - order.begin(), *bond));
      }

      std::sort(bonds.begin(), bonds.end(), compare_first<std::size_t, bond_type>());

      // call ordered_dfs_visit for all unvisited neighbors of v
      for (std::size_t i = 0; i < bonds.size(); ++i) {
        if (visitor.stop())
          return;

        bond_type bond = bonds[i].second;
        atom_type nbr = get_other(mol, bond, atom);

        if (visited[get_index(mol, nbr)]) {
          // if this bond has not been visited before, a back_bond has been found
          if (!visited[num_atoms(mol) + get_index(mol, bond)])
            visitor.back_bond(mol, bond);
          // mark bond as visited
          visited[num_atoms(mol) + get_index(mol, bond)] = true;
          continue;
        }

        // mark bond as visited
        visited[num_atoms(mol) + get_index(mol, bond)] = true;
        // invoke bond visitor
        visitor.bond(mol, prev, bond);

        ordered_dfs_visit(mol, nbr, order, visitor, visited, atom);
      }

      // invoke backtrack visitor
      visitor.backtrack(mol, atom);
    }

  } // namespace impl

  /**
   * @brief Perform a depth-first search (DFS).
   *
   * This depth-first search function considers all components and walks a single
   * spanning tree for each component. The atom that comes first in @p order from
   * each component is used as root of the resulting spanning tree to start the
   * search. If an atom has multiple neighbors the specified @p order is also used
   * to determine the order in which these are visited.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * @param mol The molecule.
   * @param order The order for visiting atoms.
   * @param visitor The DFS visitor functor.
   */
  template<typename MoleculeType, typename DFSVisitorType, typename T>
  void ordered_depth_first_search(const MoleculeType &mol, const std::vector<T> &order,
      DFSVisitorType &visitor)
  {
    visitor.initialize(mol);

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    int c = 0;
    for (std::size_t i = 0; i < order.size(); ++i) {
      if (visitor.stop())
        return;

      if (!visited[order[i]]) {
        visitor.component(c);
        impl::ordered_dfs_visit(mol, get_atom(mol, order[i]), order, visitor, visited);
      }
    }
  }

  /**
   * @brief Perform a depth-first search (DFS).
   *
   * This depth-first search function considers all components and walks a single
   * spanning tree for each component. The atom that comes first in @p order from
   * each component is used as root of the resulting spanning tree to start the
   * search. If an atom has multiple neighbors the specified @p order is also used
   * to determine the order in which these are visited.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * The @p atomMask specifies the atoms and bonds that will be considered.
   * The element in the @p atomMask that correspond with atoms to be included
   * in the search should be set to true. If an atom is not included in the
   * search, all bonds around the atom are also excluded.
   *
   * @param mol The molecule.
   * @param order The order for visiting atoms.
   * @param visitor The DFS visitor functor.
   * @param atomMask The atom mask.
   */
  template<typename MoleculeType, typename DFSVisitorType, typename T>
  void ordered_depth_first_search_mask(const MoleculeType &mol, const std::vector<T> &order,
      DFSVisitorType &visitor, const std::vector<bool> &atomMask)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

    visitor.initialize(mol);

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    // initialize visited to take mask into account
    for (std::size_t i = 0; i < num_atoms(mol); ++i)
      if (!atomMask[i]) {
        // mark masked atoms as visited
        visited[i] = true;
        // mark bonds around atom as visited
        atom_type atom = get_atom(mol, i);
        FOREACH_INCIDENT (bond, atom, mol)
          visited[num_atoms(mol) + get_index(mol, *bond)] = true;
      }

    int c = 0;
    for (std::size_t i = 0; i < order.size(); ++i) {
      if (visitor.stop())
        return;

      if (!visited[order[i]]) {
        visitor.component(c);
        impl::ordered_dfs_visit(mol, get_atom(mol, order[i]), order, visitor, visited);
      }
    }
  }

  /**
   * @brief Perform a depth-first search (DFS).
   *
   * This depth-first search function considers all components and walks a single
   * spanning tree for each component. The atom that comes first in @p order from
   * each component is used as root of the resulting spanning tree to start the
   * search. If an atom has multiple neighbors the specified @p order is also used
   * to determine the order in which these are visited.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * The @p atomMask specifies the atoms that will be considered.
   * The element in the @p atomMask that correspond with atoms to be included
   * in the search should be set to true. The @p bondMask works in the same way
   * but specifies the bonds that should be considered.
   *
   * @param mol The molecule.
   * @param order The order for visiting atoms.
   * @param visitor The DFS visitor functor.
   * @param atomMask The atom mask.
   * @param bondMask The bond mask.
   */
  template<typename MoleculeType, typename DFSVisitorType, typename T>
  void ordered_depth_first_search_mask(const MoleculeType &mol, const std::vector<T> &order,
      DFSVisitorType &visitor, const std::vector<bool> &atomMask,
      const std::vector<bool> &bondMask)
  {
    visitor.initialize(mol);

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    // initialize visited to take mask into account
    for (std::size_t i = 0; i < num_atoms(mol); ++i)
      if (!atomMask[i])
        visited[i] = true;
    for (std::size_t i = 0; i < num_bonds(mol); ++i)
      if (!bondMask[i])
        visited[num_atoms(mol) + i] = true;

    int c = 0;
    for (std::size_t i = 0; i < order.size(); ++i) {
      if (visitor.stop())
        return;

      if (!visited[order[i]]) {
        visitor.component(c);
        impl::ordered_dfs_visit(mol, get_atom(mol, order[i]), order, visitor, visited);
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // exhaustive depth first search
  //
  //////////////////////////////////////////////////////////////////////////////

  namespace impl {

    template<typename MoleculeType, typename AtomType, typename DFSVisitorType>
    void exhaustive_dfs_visit(const MoleculeType &mol, AtomType atom, DFSVisitorType &visitor,
        std::vector<bool> &visited, AtomType prev = molecule_traits<MoleculeType>::null_atom())
    {
      typedef AtomType atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      if (visitor.stop())
        return;

      // mark atom as visited
      visited[get_index(mol, atom)] = true;
      // invoke atom visitor
      visitor.atom(mol, prev, atom);

      // create list of bonds to visit next
      std::vector<bond_type> bonds;
      FOREACH_INCIDENT (bond, atom, mol) {
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
          if (visitor.stop())
            return;

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

  } // namespace impl

  /**
   * @brief Perform an exhaustive depth-first search (DFS).
   *
   * This depth-first search function considers a single component and starts
   * walking from the specified atom. If there are multiple neighbors to visit
   * next, all permutations are tried making this an exhaustive search.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. DFSAtomOrderVisitor,
   * DFSBondOrderVisitor, DFSClosureRecorderVisitor, DFSDebugVisitor). Custom
   * functor can be implemented by inheriting the DFSVisitor struct and
   * reimplementing the required functions.
   *
   * @param mol The molecule.
   * @param atom The root atom for the DFS spanning tree.
   * @param visitor The DFS visitor functor.
   */
  template<typename MoleculeType, typename AtomType, typename DFSVisitorType>
  void exhaustive_depth_first_search(const MoleculeType &mol, AtomType atom, DFSVisitorType &visitor)
  {
    visitor.initialize(mol);

    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    visitor.component(0);
    impl::exhaustive_dfs_visit(mol, atom, visitor, visited);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Depth first search visitors
  //
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @struct DFSAtomOrderVisitor algorithms/dfs.h <Helium/algorithms/dfs.h>
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
   * @struct DFSBondOrderVisitor algorithms/dfs.h <Helium/algorithms/dfs.h>
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
   * @struct DFSClosureRecorderVisitor algorithms/dfs.h <Helium/algorithms/dfs.h>
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
   * @struct DFSDebugVisitor algorithms/dfs.h <Helium/algorithms/dfs.h>
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
     * @param os_ The STL output stream to print to.
     */
    DFSDebugVisitor(std::ostream &os_ = std::cout) : os(os_)
    {
    }

    /**
     * @brief Print initialize().
     */
    void initialize(const MoleculeType &mol)
    {
      os << "initialize()" << std::endl;
    }

    /**
     * @brief Print component(i).
     */
    void component(int i)
    {
      os << "component(" << i << ")" << std::endl;
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
