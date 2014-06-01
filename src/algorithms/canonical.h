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
#ifndef HELIUM_CANONICAL_H
#define HELIUM_CANONICAL_H

#include <Helium/algorithms/invariants.h>
#include <Helium/util.h>

#include <map>
#include <limits>

#define DEBUG_CANON 0

namespace Helium {

  namespace impl {

    template<typename ItemType, typename ClassType>
    class Permutations
    {
      public:
        Permutations(const std::vector<ItemType> &items, const std::vector<ClassType> &classes)
          : m_items(items), m_classes(classes)
        {
          // find the unique classes
          std::set<ClassType> uniqueClasses;
          for (std::size_t i = 0; i < m_classes.size(); ++i)
            uniqueClasses.insert(classes[i]);

          // create the groups vector (all items with the same class are placed in the same group)
          for (typename std::set<ClassType>::iterator i = uniqueClasses.begin(); i != uniqueClasses.end(); ++i) {
            m_groups.resize(m_groups.size() + 1);
            for (std::size_t j = 0; j < m_classes.size(); ++j)
              if (*i == m_classes[j])
                m_groups.back().push_back(j);
            std::sort(m_groups.back().begin(), m_groups.back().end());
          }
        }

        bool next(std::size_t i = 0)
        {
          if (i >= m_groups.size())
            return false;
          bool more = std::next_permutation(m_groups[i].begin(), m_groups[i].end());
          if (!more)
            return next(i + 1);
          return true;
        }

        std::vector<ItemType> items() const
        {
          std::vector<ItemType> result;
          for (std::size_t i = 0; i < m_groups.size(); ++i)
            for (std::size_t j = 0; j < m_groups[i].size(); ++j)
              result.push_back(m_items[m_groups[i][j]]);
          return result;
        }

      private:
        const std::vector<ItemType> &m_items;
        const std::vector<ClassType> &m_classes;
        std::vector<std::vector<std::size_t> > m_groups; // indices of items with same class
    };

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

    template<typename MoleculeType, typename AtomInvariant, typename BondInvariant>
    class Canonicalize
    {
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
        typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;
        typedef typename molecule_traits<MoleculeType>::incident_iter incident_iter;

      public:
        Canonicalize(const MoleculeType &mol, const std::vector<unsigned long> &symmetry,
            const AtomInvariant &atomInvariant, const BondInvariant &bondInvariant)
          : m_mol(mol), m_symmetry(symmetry), m_visited(num_bonds(mol)),
            m_atomInvariant(atomInvariant), m_bondInvariant(bondInvariant)
        {
        }

        void canonicalize()
        {
          std::map<unsigned long, std::size_t> symClassCounts;
          for (std::size_t i = 0; i < m_symmetry.size(); ++i)
            symClassCounts[m_symmetry[i]]++;

          std::size_t count = std::numeric_limits<std::size_t>::max();
          unsigned long symClass = 0;
          for (std::map<unsigned long, std::size_t>::iterator i = symClassCounts.begin(); i != symClassCounts.end(); ++i)
            if (i->second == count && i->first < symClass) {
              symClass = i->first;
            } else if (i->second < count) {
              count = i->second;
              symClass = i->first;
            }

          // select atom(s) with lowest symmetry class
          atom_iter atom, end_atoms;
          TIE(atom, end_atoms) = get_atoms(m_mol);
          for (; atom != end_atoms; ++atom) {
            /*
            if (m_symmetry[get_index(m_mol, *atom)])
              continue;
            */
            if (m_symmetry[get_index(m_mol, *atom)] != symClass)
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
            TIE(bond, end_bonds) = get_bonds(m_mol, atom);
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
          // encode ATOM invariants
          //
          for (std::size_t i = 0; i < m_atoms.size(); ++i)
            code.push_back(m_atomInvariant(m_mol, get_atom(m_mol, m_atoms[i])));

          //
          // encode BOND invariants
          //
          for (std::size_t i = 0; i < m_bonds.size(); ++i)
            code.push_back(m_bondInvariant(m_mol, get_bond(m_mol, m_bonds[i])));


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
          TIE(fromBond, atom) = stack.back();
          stack.pop_back();


          if (DEBUG_CANON)
            std::cout << "next(" << get_index(m_mol, atom) << ")" << std::endl;

          // check if the atom is already labeled
          if (std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, atom)) != m_atoms.end())
            return;

          assert(std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, atom)) == m_atoms.end());

          /*
          if (fromBond != molecule_traits<MoleculeType>::null_bond()) {
            Index fromAtom = std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, get_other(m_mol, fromBond, atom))) - m_atoms.begin();

            if (m_code.size() > m_from.size()) {
              if (fromAtom > m_code[m_from.size()])
                return;
            }
          }
          */


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
            TIE(bond, end_bonds) = get_bonds(m_mol, atom);
            for (; bond != end_bonds; ++bond) {
              atom_type other = get_other(m_mol, *bond, atom);
              if (m_visited[get_index(m_mol, *bond)])
                continue;
              if (std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, other)) != m_atoms.end())
                continue;
              bonds.push_back(std::make_pair(*bond, other));
            }



            std::vector<unsigned long> classes;
            for (std::size_t i = 0; i < bonds.size(); ++i)
              classes.push_back(m_symmetry[get_index(m_mol, bonds[i].second)]);

            Permutations<std::pair<bond_type, atom_type>, unsigned long> perms(bonds, classes);
            do {
              // restore stack
              stack = stackCopy;
              // append the new bonds
              std::vector<std::pair<bond_type, atom_type> > orderedBonds = perms.items();
              std::copy(orderedBonds.begin(), orderedBonds.end(), std::back_inserter(stack));
              // recursive call...
              next(stack);
            } while (perms.next());








            /*
            // sort the new bonds
            std::sort(bonds.begin(), bonds.end(), compare_first<bond_type, atom_type>());

            // recursive call for each permutation of the new bonds
            bool last = false;
            do {
              // restore stack
              stack = stackCopy;
              // append the new bonds
              std::copy(bonds.begin(), bonds.end(), std::back_inserter(stack));
              // recursive call...
              next(stack);
            } while (!last && (last = std::next_permutation(bonds.begin(), bonds.end(), compare_first<bond_type, atom_type>())));
             */


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

        const MoleculeType &m_mol;
        const std::vector<unsigned long> &m_symmetry;
        std::vector<Index> m_atoms; // canonical atom order
        std::vector<Index> m_bonds; // canonical bond order
        std::vector<Index> m_from; // from atoms
        std::vector<bool> m_visited; // visited bonds

        std::vector<Index> m_labels; // currently lowest canonical labels
        std::vector<unsigned long> m_code; // currently lowest canonical code

        const AtomInvariant &m_atomInvariant;
        const BondInvariant &m_bondInvariant;
    };

  }

  /**
   * @brief Canonicalize a single component.
   *
   * The canonicalization algorithm consists of two
   * steps. In the first step, the atoms are partitioned using graph
   * invariants (e.g. atom degree). This initial partitioning reduces the
   * number of states that need to be visited to find the canonical code.
   *
   * In the second step all automorphic permutations of the graph are
   * investigated and a code is generated for each one. The unique code is
   * selected and the associated atom order is the canonical atom order.
   *
   * @note Complexity: @f$O(2^n)@f$
   * @ingroup Production
   * @note Phase: Production
   *
   * @param mol The molecule.
   * @param symmetry The extended connectivities (see extended_connectivities()).
   * @param atomInvariant The atom invariants to use for the canonical code
   *        (e.g. AtomElementInvariant). The used invariants determine what kind
   *        of information is considered for making the molecule canonical.
   * @param bondInvariant The bond invariants to use for the canonical code
   *        (e.g. BondOrderInvariant). The used invariants determine what kind
   *        of information is considered for making the molecule canonical.
   *
   * @return The canonical atom order and canonical code.
   */
  template<typename MoleculeType, typename T, typename AtomInvariant, typename BondInvariant>
  std::pair<std::vector<Index>, std::vector<unsigned long> >
  canonicalize_component(const MoleculeType &mol, const std::vector<T> &symmetry,
      const AtomInvariant &atomInvariant, const BondInvariant &bondInvariant)
  {
    if (DEBUG_CANON) {
      std::cout << "+---------------------------+" << std::endl;
      std::cout << "| START CAONICALIZATION     |" << std::endl;
      std::cout << "+---------------------------+" << std::endl;
      std::cout << "symmetry: " << symmetry << std::endl;
    }

    impl::Canonicalize<MoleculeType, AtomInvariant, BondInvariant> can(mol, symmetry, atomInvariant, bondInvariant);
    can.canonicalize();

    if (DEBUG_CANON) {
      std::cout << "labels: " << can.labels() << ", code: " << can.code() << std::endl;
      std::cout << "+---------------------------+" << std::endl;
      std::cout << "| START CAONICALIZATION     |" << std::endl;
      std::cout << "+---------------------------+" << std::endl;
    }

    return std::make_pair(can.labels(), can.code());
  }

  namespace impl {

    struct ComponentOrderAndCode
    {
      ComponentOrderAndCode(const std::vector<bool> &atoms_, const std::vector<Index> &order_,
          const std::vector<unsigned long> &code_) : atoms(atoms_), order(order_), code(code_)
      {
      }

      bool operator<(const ComponentOrderAndCode &other) const
      {
        return code < other.code;
      }

      std::vector<bool> atoms;
      std::vector<Index> order;
      std::vector<unsigned long> code;
    };

  } // namespace impl

  template<typename MoleculeType, typename T, typename AtomInvariant, typename BondInvariant>
  std::pair<std::vector<Index>, std::vector<unsigned long> > canonicalize(const MoleculeType &mol,
      const std::vector<T> &symmetry, const AtomInvariant &atomInvariant,
      const BondInvariant &bondInvariant, const std::vector<unsigned int> &atomComponents,
      const std::vector<unsigned int> &bondComponents)
  {
    if (DEBUG_CANON) {
      std::cout << "+---------------------------+" << std::endl;
      std::cout << "| START CAONICALIZATION     |" << std::endl;
      std::cout << "+---------------------------+" << std::endl;
      std::cout << "symmetry: " << symmetry << std::endl;
    }

    std::vector<impl::ComponentOrderAndCode> ordersAndCodes;

    Size numComponents = unique_elements(atomComponents);
    for (Size i = 0; i < numComponents; ++i) {
      std::vector<bool> atoms(num_atoms(mol));
      std::vector<bool> bonds(num_bonds(mol));
      std::vector<T> componentSymmetry;

      for (std::size_t j = 0; j < atomComponents.size(); ++j)
        if (atomComponents[j] == i) {
          atoms[j] = true;
          componentSymmetry.push_back(symmetry[j]);
        }
      for (std::size_t j = 0; j < bondComponents.size(); ++j)
        if (bondComponents[j] == i)
          bonds[j] = true;

      MoleculeType component;
      make_substructure(component, mol, atoms, bonds);

      impl::Canonicalize<MoleculeType, AtomInvariant, BondInvariant> can(component,
          componentSymmetry, atomInvariant, bondInvariant);
      can.canonicalize();

      ordersAndCodes.push_back(impl::ComponentOrderAndCode(atoms, can.labels(), can.code()));
    }

    std::sort(ordersAndCodes.begin(), ordersAndCodes.end());

    // total order & code
    std::vector<Index> order;
    std::vector<unsigned long> code;

    for (std::size_t i = 0; i < ordersAndCodes.size(); ++i) {
      const impl::ComponentOrderAndCode &component = ordersAndCodes[i];

      std::map<Index, Index> indexMap;
      for (std::size_t j = 0, k = 0; j < component.atoms.size(); ++j)
        if (component.atoms[j])
          indexMap[k++] = j;

      for (std::size_t j = 0; j < component.order.size(); ++j)
        order.push_back(indexMap[component.order[j]]);

      std::copy(component.code.begin(), component.code.end(), std::back_inserter(code));
    }

    assert(order.size() == num_atoms(mol));

    if (DEBUG_CANON) {
      std::cout << "order: " << order << std::endl;
      std::cout << "+---------------------------+" << std::endl;
      std::cout << "| START CAONICALIZATION     |" << std::endl;
      std::cout << "+---------------------------+" << std::endl;
    }

    return std::make_pair(order, code);
  }

}

#endif
