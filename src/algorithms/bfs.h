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
  };

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
   * @note Complexity: O(n)
   *
   * @param mol The molecule.
   * @param visitor The BFS visitor functor.
   */
  template<typename MoleculeType, typename BFSVisitorType>
  void breadth_first_search(const MoleculeType &mol, BFSVisitorType &visitor)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

    visitor.initialize(mol);

    std::queue<std::pair<bond_type, atom_type> > queue;
    std::vector<bool> visited(num_atoms(mol) + num_bonds(mol));

    FOREACH_ATOM (atom, mol, MoleculeType) {
      if (!visited[get_index(mol, *atom)]) {
        // push initial component atom to stack
        queue.push(std::make_pair(molecule_traits<MoleculeType>::null_bond(), *atom));

        while (queue.size()) {
          // pop from stack
          std::pair<bond_type, atom_type> next = queue.front();
          queue.pop();

          // is the atom already visited?
          if (visited[get_index(mol, next.second)])
            continue;

          // mark atom as visited
          visited[get_index(mol, next.second)] = true;
          // mark bond as visited
          if (next.first != molecule_traits<MoleculeType>::null_bond())
            visited[num_atoms(mol) + get_index(mol, next.first)] = true;
          // invoke visitor on bond and atom
          atom_type prev = get_other(mol, next.first, next.second);
          if (next.first != molecule_traits<MoleculeType>::null_bond())
            visitor.bond(mol, prev, next.first);
          visitor.atom(mol, prev, next.second);

          // add unvisited neighbors to queue
          FOREACH_INCIDENT (bond, next.second, mol, MoleculeType) {
            atom_type nbr = get_other(mol, *bond, next.second);

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
            queue.push(std::make_pair(*bond, get_other(mol, *bond, next.second)));
          }
        }
      }
    }
  }

  /**
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
