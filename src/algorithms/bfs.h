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
#ifndef HELIUM_BFS_H
#define HELIUM_BFS_H

#include <Helium/molecule.h>

#include <iostream>
#include <queue>

namespace Helium {

  /**
   * @file algorithms/bfs.h
   * @brief Breadth-first search.
   */

  /**
   * @struct BFSVisitor algorithms/bfs.h <Helium/algorithms/bfs.h>
   * @brief Base class for BFS functors.
   */
  template<typename MoleculeType>
  struct BFSVisitor
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
     * This function is called once when the BFS search is started.
     *
     * @param mol The molecule.
     */
    void initialize(const MoleculeType &mol) {}

    /**
     * @brief A new component starts.
     *
     * @param i The component (i.e. [0,n])
     */
    void component(int i) {}

    /**
     * @brief A new depth is reached.
     *
     * @param d The new depth.
     */
    void depth(int d) {}

    /**
     * @brief Visit an atom.
     *
     * This function is called every time a new atom is traversed.
     *
     * @param mol The molecule.
     * @param prev The previous atom on the BFS path.
     * @param atom The traversed atom.
     */
    void atom(const MoleculeType &mol, atom_type prev, atom_type atom) {}

    /**
     * @brief Visit a bond.
     *
     * This function is called every time a new bond is traversed.
     *
     * @param mol The molecule.
     * @param prev The previous atom on the BFS path.
     * @param bond The traversed bond.
     */
    void bond(const MoleculeType &mol, atom_type prev, bond_type bond) {}

    /**
     * @brief Visit a "back" bond.
     *
     * The "back" bonds are the chords (or ring closure bonds) that would make
     * the BFS spanning graph cyclic.
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

  namespace impl {

    template<typename MoleculeType>
    struct BFSQueueItem
    {
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      BFSQueueItem(const bond_type &bond_, const atom_type &atom_, int depth_)
        : bond(bond_), atom(atom_), depth(depth_)
      {
      }

      bond_type bond;
      atom_type atom;
      int depth;
    };

    template<typename MoleculeType, typename BFSVisitorType>
    void process_bfs_queue(const MoleculeType &mol, BFSVisitorType &visitor,
        std::queue<BFSQueueItem<MoleculeType> > &queue, std::vector<bool> &visited)
    {
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename impl::BFSQueueItem<MoleculeType> QueueItem;

      int lastDepth = -1;
      while (queue.size()) {
        if (visitor.stop())
          return;

        // pop from stack
        QueueItem next = queue.front();
        queue.pop();

        // is the atom already visited?
        if (visited[get_index(mol, next.atom)])
          continue;

        if (next.depth > lastDepth) {
          visitor.depth(next.depth);
          lastDepth = next.depth;
        }

        // mark atom as visited
        visited[get_index(mol, next.atom)] = true;
        // mark bond as visited
        if (next.bond != molecule_traits<MoleculeType>::null_bond()) {
          visited[num_atoms(mol) + get_index(mol, next.bond)] = true;
          // invoke visitor on bond and atom
          atom_type prev = get_other(mol, next.bond, next.atom);
          visitor.bond(mol, prev, next.bond);
          visitor.atom(mol, prev, next.atom);
        } else {
          visitor.atom(mol, molecule_traits<MoleculeType>::null_atom(), next.atom);
        }

        // add unvisited neighbors to queue
        FOREACH_INCIDENT (bond, next.atom, mol) {
          if (visitor.stop())
            return;

          if (visited[num_atoms(mol) + get_index(mol, *bond)])
            continue;

          atom_type nbr = get_other(mol, *bond, next.atom);

          // is the neighbor atom already visited?
          if (visited[get_index(mol, nbr)]) {
            // if the bond is not visited yet, it is a back bond
            if (!visited[num_atoms(mol) + get_index(mol, *bond)])
              visitor.back_bond(mol, *bond);
            // mark the back bond as visited
            visited[num_atoms(mol) + get_index(mol, *bond)] = true;
            continue;
          }

          // push neighbor to stack
          queue.push(QueueItem(*bond, get_other(mol, *bond, next.atom), next.depth + 1));
        }
      }
    }

  } // namespace impl

  /**
   * @brief Perform a breadth-first search (BFS).
   *
   * This breadth-first search function considers all fragments and walks a single
   * spanning tree for each fragment. The atom with the lowest index from each
   * fragment is used as root of the resulting spanning tree to start the
   * search.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. BFSAtomOrderVisitor,
   * BFSBondOrderVisitor, BFSClosureRecorderVisitor, BFSDebugVisitor). Custom
   * functor can be implemented by inheriting the BFSVisitor struct and
   * reimplementing the required functions.
   *
   * @param mol The molecule.
   * @param visitor The BFS visitor functor.
   */
  template<typename MoleculeType, typename BFSVisitorType>
  void breadth_first_search(const MoleculeType &mol, BFSVisitorType &visitor)
  {
    typedef typename impl::BFSQueueItem<MoleculeType> QueueItem;

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    visitor.initialize(mol);

    std::queue<QueueItem> queue;

    int component = 0;
    FOREACH_ATOM (atom, mol) {
      if (!visited[get_index(mol, *atom)]) {
        // push initial component atom to stack
        queue.push(QueueItem(molecule_traits<MoleculeType>::null_bond(), *atom, 0));

        // start new compoenent
        visitor.component(component++);

        // process the queue
        impl::process_bfs_queue(mol, visitor, queue, visited);
      }
    }
  }

  /**
   * @brief Perform a breadth-first search (BFS).
   *
   * This breadth-first search function considers all fragments and walks a single
   * spanning tree for each fragment. The atom with the lowest index from each
   * fragment is used as root of the resulting spanning tree to start the
   * search.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. BFSAtomOrderVisitor,
   * BFSBondOrderVisitor, BFSClosureRecorderVisitor, BFSDebugVisitor). Custom
   * functor can be implemented by inheriting the BFSVisitor struct and
   * reimplementing the required functions.
   *
   * The @p atomMask specifies the atoms and bonds that will be considered.
   * The element in the @p atomMask that correspond with atoms to be included
   * in the search should be set to true. If an atom is not included in the
   * search, all bonds around the atom are also excluded.
   *
   * @param mol The molecule.
   * @param visitor The BFS visitor functor.
   * @param atomMask The atom mask.
   */
  template<typename MoleculeType, typename BFSVisitorType>
  void breadth_first_search_mask(const MoleculeType &mol, BFSVisitorType &visitor,
      const std::vector<bool> &atomMask)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    typedef typename impl::BFSQueueItem<MoleculeType> QueueItem;

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

    visitor.initialize(mol);

    std::queue<QueueItem> queue;

    int component = 0;
    FOREACH_ATOM (atom, mol) {
      if (!visited[get_index(mol, *atom)]) {
        // push initial component atom to stack
        queue.push(QueueItem(molecule_traits<MoleculeType>::null_bond(), *atom, 0));

        // start new compoenent
        visitor.component(component++);

        // process the queue
        impl::process_bfs_queue(mol, visitor, queue, visited);
      }
    }
  }

  /**
   * @brief Perform a breadth-first search (BFS).
   *
   * This breadth-first search function considers all fragments and walks a single
   * spanning tree for each fragment. The atom with the lowest index from each
   * fragment is used as root of the resulting spanning tree to start the
   * search.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. BFSAtomOrderVisitor,
   * BFSBondOrderVisitor, BFSClosureRecorderVisitor, BFSDebugVisitor). Custom
   * functor can be implemented by inheriting the BFSVisitor struct and
   * reimplementing the required functions.
   *
   * The @p atomMask specifies the atoms that will be considered.
   * The element in the @p atomMask that correspond with atoms to be included
   * in the search should be set to true. The @p bondMask works in the same way
   * but specifies the bonds that should be considered.
   *
   * @param mol The molecule.
   * @param visitor The BFS visitor functor.
   * @param atomMask The atom mask.
   * @param bondMask The bond mask.
   */
  template<typename MoleculeType, typename BFSVisitorType>
  void breadth_first_search_mask(const MoleculeType &mol, BFSVisitorType &visitor,
      const std::vector<bool> &atomMask, const std::vector<bool> &bondMask)
  {
    typedef typename impl::BFSQueueItem<MoleculeType> QueueItem;

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    // initialize visited to take masks into account
    for (std::size_t i = 0; i < num_atoms(mol); ++i)
      if (!atomMask[i])
        visited[i] = true;
    for (std::size_t i = 0; i < num_bonds(mol); ++i)
      if (!bondMask[i])
        visited[num_atoms(mol) + i] = true;

    visitor.initialize(mol);

    std::queue<QueueItem> queue;

    int component = 0;
    FOREACH_ATOM (atom, mol) {
      if (!visited[get_index(mol, *atom)]) {
        // push initial component atom to stack
        queue.push(QueueItem(molecule_traits<MoleculeType>::null_bond(), *atom, 0));

        // start new compoenent
        visitor.component(component++);

        // process the queue
        impl::process_bfs_queue(mol, visitor, queue, visited);
      }
    }
  }

  /**
   * @brief Perform a breadth-first search (BFS).
   *
   * This breadth-first search function considers a single component and walks a
   * single spanning tree starting at @p atom.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. BFSAtomOrderVisitor,
   * BFSBondOrderVisitor, BFSClosureRecorderVisitor, BFSDebugVisitor). Custom
   * functor can be implemented by inheriting the BFSVisitor struct and
   * reimplementing the required functions.
   *
   * @param mol The molecule.
   * @param atom The start atom.
   * @param visitor The BFS visitor functor.
   */
  template<typename MoleculeType, typename BFSVisitorType>
  void breadth_first_search(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom, BFSVisitorType &visitor)
  {
    typedef typename impl::BFSQueueItem<MoleculeType> QueueItem;

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    visitor.initialize(mol);

    std::queue<QueueItem> queue;

    // push initial component atom to stack
    queue.push(QueueItem(molecule_traits<MoleculeType>::null_bond(), atom, 0));

    // start new compoenent
    visitor.component(0);

    // process the queue
    impl::process_bfs_queue(mol, visitor, queue, visited);
  }

  /**
   * @brief Perform a breadth-first search (BFS).
   *
   * This breadth-first search function considers a single component and walks a
   * single spanning tree starting at @p atom.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. BFSAtomOrderVisitor,
   * BFSBondOrderVisitor, BFSClosureRecorderVisitor, BFSDebugVisitor). Custom
   * functor can be implemented by inheriting the BFSVisitor struct and
   * reimplementing the required functions.
   *
   * The @p atomMask specifies the atoms and bonds that will be considered.
   * The element in the @p atomMask that correspond with atoms to be included
   * in the search should be set to true. If an atom is not included in the
   * search, all bonds around the atom are also excluded.
   *
   * @param mol The molecule.
   * @param atom The start atom.
   * @param visitor The BFS visitor functor.
   * @param atomMask The atom mask.
   */
  template<typename MoleculeType, typename BFSVisitorType>
  void breadth_first_search_mask(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom, BFSVisitorType &visitor,
      const std::vector<bool> &atomMask)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    typedef typename impl::BFSQueueItem<MoleculeType> QueueItem;

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

    visitor.initialize(mol);

    std::queue<QueueItem> queue;

    // push initial component atom to stack
    queue.push(QueueItem(molecule_traits<MoleculeType>::null_bond(), atom, 0));

    // start new compoenent
    visitor.component(0);

    // process the queue
    impl::process_bfs_queue(mol, visitor, queue, visited);
  }

  /**
   * @brief Perform a breadth-first search (BFS).
   *
   * This breadth-first search function considers a single component and walks a
   * single spanning tree starting at @p atom.
   *
   * This function makes use of a visitor functor to allow actions to be
   * performed. A number of visitors are available (e.g. BFSAtomOrderVisitor,
   * BFSBondOrderVisitor, BFSClosureRecorderVisitor, BFSDebugVisitor). Custom
   * functor can be implemented by inheriting the BFSVisitor struct and
   * reimplementing the required functions.
   *
   * The @p atomMask specifies the atoms that will be considered.
   * The element in the @p atomMask that correspond with atoms to be included
   * in the search should be set to true. The @p bondMask works in the same way
   * but specifies the bonds that should be considered.
   *
   * @param mol The molecule.
   * @param atom The start atom.
   * @param visitor The BFS visitor functor.
   * @param atomMask The atom mask.
   * @param bondMask The bond mask.
   */
  template<typename MoleculeType, typename BFSVisitorType>
  void breadth_first_search_mask(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom, BFSVisitorType &visitor,
      const std::vector<bool> &atomMask, const std::vector<bool> &bondMask)
  {
    typedef typename impl::BFSQueueItem<MoleculeType> QueueItem;

    // keep track of visited atoms and bonds
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    // initialize visited to take masks into account
    for (std::size_t i = 0; i < num_atoms(mol); ++i)
      if (!atomMask[i])
        visited[i] = true;
    for (std::size_t i = 0; i < num_bonds(mol); ++i)
      if (!bondMask[i])
        visited[num_atoms(mol) + i] = true;

    visitor.initialize(mol);

    std::queue<QueueItem> queue;

    // push initial component atom to stack
    queue.push(QueueItem(molecule_traits<MoleculeType>::null_bond(), atom, 0));

    // start new compoenent
    visitor.component(0);

    // process the queue
    impl::process_bfs_queue(mol, visitor, queue, visited);
  }

  /**
   * @struct BFSAtomOrderVisitor algorithms/bfs.h <Helium/algorithms/bfs.h>
   * @brief A BFS visitor that records the order in which the atoms are visited.
   */
  template<typename MoleculeType>
  struct BFSAtomOrderVisitor : public BFSVisitor<MoleculeType>
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
     * @param prev The previous atom on the BFS path.
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
   * @struct BFSBondOrderVisitor algorithms/bfs.h <Helium/algorithms/bfs.h>
   * @brief A BFS visitor that records the order in which the bonds are visited.
   */
  template<typename MoleculeType>
  struct BFSBondOrderVisitor : public BFSVisitor<MoleculeType>
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
     * @param prev The previous atom on the BFS path.
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
   * @struct BFSClosureRecorderVisitor algorithms/bfs.h <Helium/algorithms/bfs.h>
   * @brief A BFS visitor that records the "back" bonds.
   */
  template<typename MoleculeType>
  struct BFSClosureRecorderVisitor : public BFSVisitor<MoleculeType>
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
   * @struct BFSDebugVisitor algorithms/bfs.h <Helium/algorithms/bfs.h>
   * @brief A BFS visitor that prints out debug information.
   */
  template<typename MoleculeType>
  struct BFSDebugVisitor : public BFSVisitor<MoleculeType>
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
    BFSDebugVisitor(std::ostream &os_ = std::cout) : os(os_)
    {
    }

    /**
     * @brief Initialize the functor.
     *
     * This function is called once when the BFS search is started.
     *
     * @param mol The molecule.
     */
    void initialize(const MoleculeType &mol)
    {
      os << "initialize()" << std::endl;
    }

    /**
     * @brief A new component starts.
     *
     * @param i The component (i.e. [0,n])
     */
    void component(int i)
    {
      os << "component(" << i << ")" << std::endl;
    }

    /**
     * @brief A new depth is reached.
     *
     * @param d The new depth.
     */
    void depth(int d)
    {
      os << "depth(" << d << ")" << std::endl;
    }

    /**
     * @brief Visit an atom.
     *
     * This function prints the string "atom(i)" to os where i is the
     * atom index.
     *
     * @param mol The molecule.
     * @param prev The previous atom on the BFS path.
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
     * @param prev The previous atom on the BFS path.
     * @param bond The traversed bond.
     */
    void bond(const MoleculeType &mol, atom_type prev, bond_type bond)
    {
      os << "bond(" << get_index(mol, bond) << ")" << std::endl;
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
