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
#include <Helium/algorithms/extendedconnectivities.h>
#include <Helium/algorithms/cycles.h>
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

        unsigned long complexity() const
        {
          unsigned long result = 1;
          for (std::size_t i = 0; i < m_groups.size(); ++i)
            result *= factorial(m_groups[i].size());
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

    /*
    class Automorphism
    {
      public:
        Automorphism(const std::vector<std::vector<unsigned int> > &orbits)
        {
          // number of orbits
          m_data.push_back(orbits.size());
          // orbit sizes
          for (std::size_t i = 0; i < orbits.size(); ++i)
            m_data.push_back(orbits[i].size());
          // orbits
          for (std::size_t i = 0; i < orbits.size(); ++i)
            std::copy(orbits[i].begin(), orbits[i].end(), std::back_inserter(m_data));
          // sort orbits
          std::size_t begin = 1 + orbits.size();
          for (unsigned int i = 0; i < orbits.size(); ++i) {
            std::sort(m_data.begin() + begin, m_data.begin() + begin + orbits[i].size());
            begin += orbits[i].size();
          }
        }

        std::size_t numOrbits() const
        {
          return m_data[0];
        }

        unsigned int orbitSize(unsigned int i) const
        {
          return m_data[1 + i];
        }

        const unsigned int* orbit(unsigned int i) const
        {
          unsigned int offset = 1 + m_data[0];
          for (unsigned int j = 0; j < i; ++j)
            offset += m_data[1 + j];
          return &m_data[0] + offset;
        }

        bool isAutomorphic(unsigned int v, unsigned int w) const
        {
          std::size_t begin = 1 + m_data[0];
          for (unsigned int i = 0; i < m_data[0]; ++i) {
            if (std::binary_search(m_data.begin() + begin, m_data.begin() + begin + m_data[1 + i], v) &&
                std::binary_search(m_data.begin() + begin, m_data.begin() + begin + m_data[1 + i], w))
              return true;
            begin += m_data[1 + i];
          }
          return false;
        }

      private:
        friend std::ostream& operator<<(std::ostream &os, const Automorphism &automorphism);

        // [n] [sizes] [orbits]
        std::vector<unsigned int> m_data;
    };

    inline std::ostream& operator<<(std::ostream &os, const Automorphism &automorphism)
    {
      std::size_t begin = 1 + automorphism.m_data[0];
      for (unsigned int i = 0; i < automorphism.m_data[0]; ++i) {
        os << "( ";
        for (unsigned int j = begin; j < begin + automorphism.m_data[1 + i]; ++j)
          os << automorphism.m_data[j] << " ";
        os << ") ";
        begin += automorphism.m_data[1 + i];
      }
      return os;
    }

    class Automorphisms
    {
      public:
        void add(const std::vector<unsigned int> &p1, const std::vector<unsigned int> &p2)
        {
          std::vector<std::vector<unsigned int> > orbits;
          std::vector<bool> done(p1.size());
          for (std::size_t i = 0; i < p1.size(); ++i) {
            if (done[p1[i]])
              continue;
            if (p1[i] == p2[i])
              continue;

            std::vector<unsigned int> orbit;
            orbit.push_back(p1[i]);
            done[p1[i]] = true;

            std::size_t j = i;
            while (p2[j] != p1[i]) {
              done[p2[j]] = true;
              orbit.push_back(p2[j]);
              j = std::find(p1.begin(), p1.end(), p2[j]) - p1.begin();
            }

            orbits.push_back(orbit);
          }

          int MAX_AUTO = 250;
          if (m_automorphisms.size() >= MAX_AUTO)
            m_automorphisms[rand() % MAX_AUTO] = Automorphism(orbits);
          else
            m_automorphisms.push_back(Automorphism(orbits));
          //removeRedundant();
        }

        const std::vector<Automorphism>& automorphisms() const
        {
          return m_automorphisms;
        }

      private:
        bool isRedundantOrbit(const unsigned int *orbit1, unsigned int size1,
            const unsigned *orbit2, unsigned int size2) const
        {
          if (size1 > size2)
            return false;

          for (unsigned int i = 0; i < size1; ++i)
            if (!std::binary_search(orbit2, orbit2 + size2, orbit1[i]))
              return false;

          return true;
        }

        void removeRedundant()
        {
          std::vector<std::size_t> remove;
          for (std::size_t i = 0; i < m_automorphisms.size(); ++i) {
            const Automorphism &aut1 = m_automorphisms[i];
            for (std::size_t j = 0; j < m_automorphisms.size(); ++j) {
              if (i == j)
                continue;

              const Automorphism &aut2 = m_automorphisms[m_automorphisms.size() - j - 1];

              bool identical = true;
              int redundant = 0;
              for (std::size_t k = 0; k < aut1.numOrbits(); ++k) {
                for (std::size_t l = 0; l < aut2.numOrbits(); ++l) {
                  if (isRedundantOrbit(aut1.orbit(k), aut1.orbitSize(k), aut2.orbit(l), aut2.orbitSize(l))) {
                    redundant++;
                    if (aut1.orbitSize(k) < aut2.orbitSize(l))
                      identical = false;
                    break;
                  }
                }
                if (redundant != k + 1)
                  break;
              }

              if (i < j && identical)
                continue;

              if (redundant == aut1.numOrbits())
                remove.push_back(i);
            }
          }

          while (remove.size()) {
            //std::cout << "removing automorphism " << remove.back() << std::endl;
            m_automorphisms.erase(m_automorphisms.begin() + remove.back());
            remove.pop_back();
          }
        }

        std::vector<Automorphism> m_automorphisms;
    };
    */

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
            m_atomInvariant(atomInvariant), m_bondInvariant(bondInvariant),
            m_backtrackDepth(0)
        {
        }

        void canonicalize()
        {
          std::vector<bool> cyclicBonds;
          cycle_membership(m_mol, m_cyclicAtoms, cyclicBonds);

          // select a symmetry class to start
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

          // initiate the dfs search
          FOREACH_ATOM (atom, m_mol) {
            if (m_symmetry[get_index(m_mol, *atom)] != symClass)
              continue;

            assert(std::find(m_visited.begin(), m_visited.end(), true) == m_visited.end());
            assert(m_atoms.empty());
            assert(m_bonds.empty());
            assert(m_from.empty());
            assert(m_code.empty());

            m_atoms.push_back(get_index(m_mol, *atom));
            m_code.push_back(m_atomInvariant(m_mol, *atom));
            next(*atom);
            m_atoms.pop_back();
            m_code.pop_back();
          }

         /*
          std::cout << "automorphisms: " << std::endl;
          for (std::size_t i = 0; i < m_automorphisms.automorphisms().size(); ++i)
            std::cout << "    " << m_automorphisms.automorphisms()[i] << std::endl;
          */
        }

        const std::vector<Index>& labels() const
        {
          return m_canAtoms;
        }

        const std::vector<unsigned long>& code() const
        {
          return m_canCode;
        }

      private:
        /*
        bool isAutomorphic(const std::vector<unsigned int> v1,
            const std::vector<unsigned int> &v2) const
        {
          assert(v1.size() >= v2.size());

          for (std::size_t i = 0; i < m_automorphisms.automorphisms().size(); ++i) {
            const Automorphism &automorphism = m_automorphisms.automorphisms()[i];

            bool automorphic = true;
            for (std::size_t j = 0; j < v2.size(); ++j)
              if (!automorphism.isAutomorphic(v1[j], v2[j])) {
                automorphic = false;
                break;
              }

            if (automorphic)
              return true;
          }

          return false;
        }
        */

        std::vector<unsigned long> createCode()
        {
          std::vector<unsigned long> code;

          //
          // encode ATOM invariants
          //
          for (std::size_t i = 0; i < m_atoms.size(); ++i)
            code.push_back(m_atomInvariant(m_mol, get_atom(m_mol, m_atoms[i])));

          //
          // FROM atoms (encodes spanning tree)
          //
          std::copy(m_from.begin(), m_from.end(), std::back_inserter(code));

          //
          // RING-CLOSURES (encodes ring closures)
          //
          unsigned int numClosures = 0;
          for (std::size_t i = 0; i < m_atoms.size(); ++i) {
            atom_type atom = get_atom(m_mol, m_atoms[i]);
            // still need to sort [1 3] and [1 4]
            std::vector<Closure> closures; // [(bond index, other atom index)]

            FOREACH_INCIDENT (bond, atom, m_mol) {
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

          return code;
        }

        void next(const atom_type &atom)
        {
          if (DEBUG_CANON)
            std::cout << "next(" << get_index(m_mol, atom) << ")" << std::endl;

          if (m_backtrackDepth) {
            if (m_atoms.size() > m_backtrackDepth)
              return;
            if (m_atoms.size() == m_backtrackDepth) {
              m_backtrackDepth = 0;
              return;
            }
            if (m_atoms.size() < m_backtrackDepth)
              m_backtrackDepth = 0;
          }



          if (m_atoms.size() == num_atoms(m_mol)) {
            // found a mapping
            if (DEBUG_CANON)
              std::cout << "mapping: " << m_atoms << ", from: " << m_from << std::endl;

            std::vector<unsigned long> code = createCode();

            if (m_canCode.empty()) {
              m_firstAtoms = m_atoms;
              m_firstCode = code;
              m_canAtoms = m_atoms;
              m_canCode = code;
            } else {
              if (code < m_canCode) {
                m_canAtoms = m_atoms;
                m_canCode = code;
              } else if (code == m_canCode) {
                //m_automorphisms.add(m_canAtoms, m_atoms);
                m_backtrackDepth = 0;
                for (std::size_t i = 0; i < m_atoms.size(); ++i) {
                  if (m_atoms[i] != m_canAtoms[i])
                    return;
                  m_backtrackDepth++;
                }
              }

              if (code == m_firstCode) {
                //m_automorphisms.add(m_firstAtoms, m_atoms);
                m_backtrackDepth = 0;
                for (std::size_t i = 0; i < m_atoms.size(); ++i) {
                  if (m_atoms[i] != m_firstAtoms[i])
                    return;
                  m_backtrackDepth++;
                }
              }
            }


          } else {
            std::vector<std::pair<bond_type, atom_type> > bonds;

            // find all unvisited bonds around the current atom
            FOREACH_INCIDENT (bond, atom, m_mol) {
              atom_type other = get_other(m_mol, *bond, atom);
              if (m_visited[get_index(m_mol, *bond)])
                continue;
              if (std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, other)) != m_atoms.end())
                continue;
              bonds.push_back(std::make_pair(*bond, other));
            }

            if (bonds.empty()) {
              for (std::size_t i = 0; i < m_atoms.size(); ++i) {
                atom_type nextAtom = get_atom(m_mol, m_atoms[i]);
                FOREACH_INCIDENT (bond, nextAtom, m_mol) {
                  atom_type other = get_other(m_mol, *bond, nextAtom);
                  if (m_visited[get_index(m_mol, *bond)])
                    continue;
                  if (std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, other)) != m_atoms.end())
                    continue;
                  next(nextAtom);
                  return;
                }
              }
              assert(0);
              return;
            }

            // find the symmetry classes for these bonds' other atom
            std::vector<unsigned long> classes;
            for (std::size_t i = 0; i < bonds.size(); ++i)
              classes.push_back(m_symmetry[get_index(m_mol, bonds[i].second)]);

            // generate all permutations of these bonds with the same symmetry class
            Permutations<std::pair<bond_type, atom_type>, unsigned long> perms(bonds, classes);
            //std::cout << "complexity: " << perms.complexity() << std::endl;
            do {
              // append the new bonds
              std::vector<std::pair<bond_type, atom_type> > orderedBonds = perms.items();

              // check if this order is automorphic
              /*
              if (perms.complexity() > 10) {
                m_atoms.push_back(get_index(m_mol, orderedBonds[0].second));

                bool automorphic = false;
                if (m_canAtoms.size() && isAutomorphic(m_canAtoms, m_atoms)) {
                  automorphic = true;
                  //std::cout << "prune automorphism 1..." << std::endl;
                }
                if (m_firstAtoms.size() && isAutomorphic(m_firstAtoms, m_atoms)) {
                  automorphic = true;
                  //std::cout << "prune automorphism 2..." << std::endl;
                }
                m_atoms.pop_back();

                if (automorphic) {
                  std::cout << "PRUNED" << std::endl;
                  continue;
                }
              }
              */

              // update m_code
              for (std::size_t i = 0; i < orderedBonds.size(); ++i)
                m_code.push_back(m_atomInvariant(m_mol, orderedBonds[i].second));

              if (m_canCode.size() && m_code > m_canCode) {
                //std::cout << "prune code..." << std::endl;
                m_code.resize(m_code.size() - orderedBonds.size());
                continue;
              }

              for (std::size_t i = 0; i < orderedBonds.size(); ++i) {
                bond_type fromBond = orderedBonds[i].first;
                atom_type nbr = orderedBonds[i].second;

                // map the atom
                m_atoms.push_back(get_index(m_mol, nbr));
                if (fromBond != molecule_traits<MoleculeType>::null_bond()) {
                  m_bonds.push_back(get_index(m_mol, fromBond));
                  // add from atom
                  m_from.push_back(std::find(m_atoms.begin(), m_atoms.end(), get_index(m_mol, atom)) - m_atoms.begin());
                  // mark bond as visited
                  m_visited[get_index(m_mol, fromBond)] = true;
                }
              }

              // recursive call...
              next(atom);

              // backtrack...
              for (std::size_t i = 0; i < orderedBonds.size(); ++i) {
                m_atoms.pop_back();

                bond_type fromBond = orderedBonds[i].first;
                if (fromBond != molecule_traits<MoleculeType>::null_bond()) {
                  m_bonds.pop_back();
                  m_from.pop_back();
                  m_visited[get_index(m_mol, fromBond)] = false;
                }
              }

              m_code.resize(m_code.size() - orderedBonds.size());

              if (m_atoms.size() <= m_backtrackDepth)
                  m_backtrackDepth = 0;

              if (!m_cyclicAtoms[get_index(m_mol, atom)] || is_aromatic(m_mol, atom))
                break;

            } while (perms.next());

          }

        }

        const MoleculeType &m_mol;
        const std::vector<unsigned long> &m_symmetry;
        std::vector<bool> m_cyclicAtoms;

        // current
        std::vector<bool> m_visited; // visited bonds
        std::vector<Index> m_atoms; // canonical atom order
        std::vector<Index> m_bonds; // canonical bond order
        std::vector<Index> m_from; // from atoms
        std::vector<unsigned long> m_code;

        // first
        std::vector<Index> m_firstAtoms;
        std::vector<unsigned long> m_firstCode;

        // canonical
        std::vector<Index> m_canAtoms; // currently lowest canonical labels
        std::vector<unsigned long> m_canCode; // currently lowest canonical code

        const AtomInvariant &m_atomInvariant;
        const BondInvariant &m_bondInvariant;

        //Automorphisms m_automorphisms;
        std::size_t m_backtrackDepth;
    };

    template<typename MoleculeType>
    void remove_terminal_symmetry(const MoleculeType &mol, std::vector<unsigned long> &symmetry)
    {
      if (symmetry.empty())
        return;

      unsigned int numColors = *std::max_element(symmetry.begin(), symmetry.end());

      //std::cout << "before: " << symmetry << std::endl;

      for (unsigned int color = 0; color < numColors; ++color) {
        std::vector<std::map<unsigned int, std::vector<Index> > > terminals;
        FOREACH_ATOM (atom, mol) {
          if (symmetry[get_index(mol, *atom)] != color)
            continue;
          terminals.resize(terminals.size() + 1);
          FOREACH_NBR (nbr, *atom, mol) {
            if (get_degree(mol, *nbr) == 1)
              terminals.back()[symmetry[get_index(mol, *nbr)]].push_back(get_index(mol, *nbr));
          }
        }

        unsigned int nextColor = *std::max_element(symmetry.begin(), symmetry.end());

        for (std::size_t i = 0; i < terminals.size(); ++i) {
          unsigned int newColor = nextColor;
          std::map<unsigned int, std::vector<Index> >::const_iterator terminal;
          for (terminal = terminals[i].begin(); terminal != terminals[i].end(); ++terminal) {
            //std::cout << terminal->first << " -> " << terminal->second << std::endl;
            for (std::size_t j = 1; j < terminal->second.size(); ++j)
              symmetry[terminal->second[j]] = newColor++;
          }
        }
      }

      //std::cout << "after:  " << symmetry << std::endl;

    }

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

    std::vector<unsigned long> sym = symmetry;
    impl::extended_connectivities_renumber(sym);
    impl::remove_terminal_symmetry(mol, sym);

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

    MoleculeType component;
    std::vector<impl::ComponentOrderAndCode> ordersAndCodes;
    Size numComponents = unique_elements(atomComponents);

    // canonicalize each component separatly
    for (Size i = 0; i < numComponents; ++i) {
      std::vector<bool> atoms(num_atoms(mol));
      std::vector<bool> bonds(num_bonds(mol));
      std::vector<T> componentSymmetry;

      // create atom & bond mask for the component
      for (std::size_t j = 0; j < atomComponents.size(); ++j)
        if (atomComponents[j] == i) {
          atoms[j] = true;
          // also isolate the symmetry for the fragment
          componentSymmetry.push_back(symmetry[j]);
        }
      for (std::size_t j = 0; j < bondComponents.size(); ++j)
        if (bondComponents[j] == i)
          bonds[j] = true;

      // create a molecule for the component
      make_substructure(component, mol, atoms, bonds);

      impl::extended_connectivities_renumber(componentSymmetry);
      impl::remove_terminal_symmetry(component, componentSymmetry);

      impl::Canonicalize<MoleculeType, AtomInvariant, BondInvariant> can(component,
          componentSymmetry, atomInvariant, bondInvariant);
      can.canonicalize();

      assert(can.labels().size() == num_atoms(component));

      ordersAndCodes.push_back(impl::ComponentOrderAndCode(atoms, can.labels(), can.code()));
    }

    // sort the canonical codes for the components
    std::sort(ordersAndCodes.begin(), ordersAndCodes.end());

    // construct total order & code
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
