/**
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
#ifndef HELIUM_CANONICAL_H
#define HELIUM_CANONICAL_H

#include <Helium/substructure.h>
#include <Helium/invariants.h>
#include <Helium/util.h>

#define DEBUG_CANON 0

namespace Helium {

  namespace impl {

    struct Closure
    {
      struct compare
      {
        bool operator()(const Closure &a, const Closure &b) const
        {
          if (a.source < b.source)
            return true;
          if (a.target < b.target)
            return true;
          return false;
        }
      };

      Closure(Index bond_, Index source_, Index target_) : bond(bond_), source(source_), target(target_)
      {
      }

      Index bond; // bond index
      Index source; // canonical source index
      Index target; // canonical target index
    };

    template<typename MoleculeType>
    class Canonicalize
    {
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
        typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;
        typedef typename molecule_traits<MoleculeType>::incident_iter incident_iter;

      public:
        Canonicalize(MoleculeType &mol, const std::vector<unsigned long> &symmetry)
            : m_mol(mol), m_symmetry(symmetry), m_visited(num_bonds(mol))
        {
        }

        void canonicalize()
        {
          // select atom(s) with lowest symmetry class
          atom_iter atom, end_atoms;
          tie(atom, end_atoms) = get_atoms(m_mol);
          for (; atom != end_atoms; ++atom) {
            if (m_symmetry[get_index(m_mol, *atom)])
              continue;

            assert(m_atoms.empty());
            assert(m_bonds.empty());
            assert(m_from.empty());
            assert(std::find(m_visited.begin(), m_visited.end(), true) == m_visited.end());

            std::vector<std::pair<bond_type, atom_type> > stack(1, std::make_pair(molecule_traits<MoleculeType>::null_bond(), *atom));
            next(stack);
          }
        }

        const std::vector<Index>& labels() const
        {
          return m_labels;
        }

        const std::vector<unsigned long>& code() const
        {
          return m_code;
        }

      private:
        void createCode()
        {
          std::vector<unsigned long> code;
          //
          // encode graph
          //

          // FROM atoms (encodes spanning tree)
          std::copy(m_from.begin(), m_from.end(), std::back_inserter(code));

          // RING-CLOSURES (encodes ring closures)
          unsigned int numClosures = 0;
          for (std::size_t i = 0; i < m_atoms.size(); ++i) {
            atom_type atom = get_atom(m_mol, m_atoms[i]);
            // still need to sort [1 3] and [1 4]
            std::vector<Closure> closures; // [(bond index, other atom index)]

            incident_iter bond, end_bonds;
            tie(bond, end_bonds) = get_bonds(m_mol, atom);
            for (; bond != end_bonds; ++bond)
              // a closure bond is a bond not found while generating the FROM spanning tree.
              if (!m_visited[get_index(m_mol, *bond)]) {
                atom_type other = get_other(m_mol, *bond, atom);
                Index source = std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, atom)) - m_atoms.begin();
                Index target = std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, other)) - m_atoms.begin();
                if (target < source)
                  std::swap(source, target);
                closures.push_back(Closure(get_index(m_mol, *bond), source, target));
                m_visited[get_index(m_mol, *bond)] = true;
              }

            // do the sorting: [1 3] < [1 4]
            std::sort(closures.begin(), closures.end(), Closure::compare());
            numClosures += closures.size();

            for (std::size_t j = 0; j < closures.size(); ++j) {
              // add the closure bond to the code
              code.push_back(closures[j].source);
              code.push_back(closures[j].target);
              //std::cout << "closure: " << code[code.size() - 2] << " " << code.back() << std::endl;
              // add the bond to the list (needed for BOND-TYPES below)
              m_bonds.push_back(closures[j].bond);
            }
          }

          //
          // encode ATOM attributes
          //
          for (std::size_t i = 0; i < m_atoms.size(); ++i)
            //code.push_back(atom_invariant(m_mol, get_atom(m_mol, m_atoms[i])));
            code.push_back(get_element(m_mol, get_atom(m_mol, m_atoms[i])));

          //
          // encode BOND attributes
          //
          for (std::size_t i = 0; i < m_bonds.size(); ++i)
            //code.push_back(bond_invariant(m_mol, get_bond(m_mol, m_bonds[i])));
            code.push_back(1);

          // backtrack closure bonds
          for (unsigned int i = 0; i < numClosures; ++i) {
            m_visited[m_bonds.back()] = false;
            m_bonds.pop_back();
          }

          if (DEBUG_CANON)
            std::cout << "code: " << code << std::endl;

          if (m_code.empty() || code < m_code) {
            m_labels = m_atoms;
            m_code = code;
          }
        }

        void next(std::vector<std::pair<bond_type, atom_type> > &stack)
        {
          if (DEBUG_CANON)
            std::cout << "stack: " << stack << std::endl;

          if (stack.empty())
            return;

          // pop next atom from stack
          atom_type atom;
          bond_type fromBond;
          tie(fromBond, atom) = stack.back();
          stack.pop_back();


          if (DEBUG_CANON)
            std::cout << "next(" << get_index(m_mol, atom) << ")" << std::endl;

          if (std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, atom)) != m_atoms.end())
            return;

          assert(std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, atom)) == m_atoms.end());

          // map the atom
          m_atoms.push_back(get_index(m_mol, atom));
          if (fromBond != molecule_traits<MoleculeType>::null_bond()) {
            m_bonds.push_back(get_index(m_mol, fromBond));
            // add from atom
            m_from.push_back(std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, get_other(m_mol, fromBond, atom))) - m_atoms.begin());
            // mark bond as visited
            m_visited[get_index(m_mol, fromBond)] = true;
          }


          if (m_atoms.size() == num_atoms(m_mol)) {
            // found a mapping
            if (DEBUG_CANON)
              std::cout << "mapping: " << m_atoms << ", from: " << m_from << std::endl;
            createCode();
          } else {

            std::vector<std::pair<bond_type, atom_type> > stackCopy(stack);
            std::vector<std::pair<bond_type, atom_type> > bonds;


            // append unvisited bonds around atom to stack
            incident_iter bond, end_bonds;
            tie(bond, end_bonds) = get_bonds(m_mol, atom);
            for (; bond != end_bonds; ++bond) {
              atom_type other = get_other(m_mol, *bond, atom);
              if (m_visited[get_index(m_mol, *bond)])
                continue;
              if (std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, other)) != m_atoms.end())
                continue;
              bonds.push_back(std::make_pair(*bond, other));
            }


            // sort the new bonds
            std::sort(bonds.begin(), bonds.end(), compare_first<bond_type, atom_type>());

            // recursive call for each permutation of the new bonds
            bool last = false;
            do {
              stack = stackCopy;
              std::copy(bonds.begin(), bonds.end(), std::back_inserter(stack));
              next(stack);
            } while (!last && (last = std::next_permutation(bonds.begin(), bonds.end(), compare_first<bond_type, atom_type>())));

            // process stack
            while (!stack.empty())
              next(stack);

          }

          if (DEBUG_CANON)
            std::cout << "backtrack..." << std::endl;

          // backtrack
          m_atoms.pop_back();
          if (fromBond != molecule_traits<MoleculeType>::null_bond()) {
            m_bonds.pop_back();
            m_from.pop_back();
            m_visited[get_index(m_mol, fromBond)] = false;
          }
        }

        MoleculeType &m_mol;
        const std::vector<unsigned long> &m_symmetry;
        std::vector<Index> m_atoms; // canonical atom order
        std::vector<Index> m_bonds; // canonical bond order
        std::vector<Index> m_from; // from atoms
        std::vector<bool> m_visited; // visited bonds

        std::vector<Index> m_labels; // currently lowest canonical labels
        std::vector<unsigned long> m_code; // currently lowest canonical code
    };

  }


  /**
   * Canonicalize a molecule. The canonicalization algorithm consists of two
   * steps. In the first step, the atoms are partitioned using graph
   * invariants (e.g. atom degree). This initial partitioning reduces the
   * number of states that need to be visited to find the canonical code.
   *
   * In the second step all automorphic permutations of the graph are
   * investigated and a code is generated for each one. The unique code is
   * selected and the associated atom order is the canonical atom order.
   *
   * @return The canonical atom order and canonical code.
   */
  template<typename MoleculeType, typename T>
  std::pair<std::vector<Index>, std::vector<unsigned long> > canonicalize(MoleculeType &mol, const std::vector<T> &symmetry)
  {
    if (DEBUG_CANON) {
      std::cout << "+---------------------------+" << std::endl;
      std::cout << "| START CAONICALIZATION     |" << std::endl;
      std::cout << "+---------------------------+" << std::endl;
      std::cout << "symmetry: " << symmetry << std::endl;
    }
    impl::Canonicalize<MoleculeType> can(mol, symmetry);
    can.canonicalize();
    if (DEBUG_CANON) {
      std::cout << "labels: " << can.labels() << ", code: " << can.code() << std::endl;
      std::cout << "+---------------------------+" << std::endl;
      std::cout << "| START CAONICALIZATION     |" << std::endl;
      std::cout << "+---------------------------+" << std::endl;
    }
    return std::make_pair(can.labels(), can.code());
  }

}

#endif
