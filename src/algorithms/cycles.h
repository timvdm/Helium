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
#ifndef HELIUM_CYCLES_H
#define HELIUM_CYCLES_H

#include <Helium/algorithms/components.h>
#include <Helium/algorithms/dfs.h>
#include <Helium/algorithms/isomorphism.h>
#include <Helium/smiles.h>
#include <Helium/ring.h>

namespace Helium {

  /**
   * @brief Get the cyclomatic number.
   *
   * The cyclomatic number is defined by \f$m - n + c\f$ where \f$n\f$ is the
   * number of atoms, \f$m\f$ the number of bonds and \f$c\f$ the number of
   * connected components. This formula is known as Cauchy's formula. The
   * cyclomatic number is the same as the nullity or first Betti's number.
   *
   * @param mol The molecule.
   * @param numComponents The number of connected components.
   *
   * @return The cyclomatic number.
   */
  template<typename MoleculeType>
  Size cyclomatic_number(MoleculeType &mol, Size numComponents)
  {
    return num_bonds(mol) - num_atoms(mol) + numComponents;
  }

  /**
   * @brief Get the cyclomatic number.
   *
   * The cyclomatic number is defined by \f$m - n + c\f$ where \f$n\f$ is the
   * number of atoms, \f$m\f$ the number of bonds and \f$c\f$ the number of
   * connected components. This formula is known as Cauchy's formula. The
   * cyclomatic number is the same as the nullity or first Betti's number.
   *
   * @param mol The molecule.
   *
   * @return The cyclomatic number.
   */
  template<typename MoleculeType>
  Size cyclomatic_number(MoleculeType &mol)
  {
    return cyclomatic_number(mol, num_connected_components(mol));
  }

  namespace impl {

    template<typename MoleculeType>
    struct CycleMembershipVisitor : public DFSVisitor<MoleculeType>
    {
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      CycleMembershipVisitor(std::vector<bool> &cycle_atoms_,
          std::vector<bool> &cycle_bonds_) : cycle_atoms(cycle_atoms_),
          cycle_bonds(cycle_bonds_)
      {
      }

      void initialize(const MoleculeType &mol)
      {
        path.reserve(num_atoms(mol));
      }

      void atom(const MoleculeType &mol, atom_type prev, atom_type atom)
      {
        path.push_back(get_index(mol, atom));
      }

      void backtrack(const MoleculeType &mol, atom_type atom)
      {
        path.pop_back();
      }

      void back_bond(const MoleculeType &mol, bond_type bond)
      {
        // find the atom that ends the cycle
        atom_type last_atom = get_other(mol, bond, get_atom(mol, path.back()));
        cycle_atoms[get_index(mol, last_atom)] = true;
        cycle_bonds[get_index(mol, bond)] = true;
        // traceback the path to last_atom
        for (std::size_t i = 0; i < path.size(); ++i) {
          atom_type atom = get_atom(mol, path[path.size() - i - 1]);
          // mark bond as cyclic
          if (i)
            cycle_bonds[get_index(mol, get_bond(mol, atom, get_atom(mol, path[path.size() - i])))] = true;
          // stop when the atom closing the cycle is found
          if (atom == last_atom)
            break;
          // mark vertex as cyclic
          cycle_atoms[get_index(mol, atom)] = true;
        }
      }

      std::vector<Index> path;
      std::vector<bool> &cycle_atoms;
      std::vector<bool> &cycle_bonds;
    };

  }

  /**
   * @brief Determine atom and bond cycle membership.
   *
   * @note Complexity: O(n)
   *
   * @param mol The molecule.
   * @param cycle_atoms cyclic atoms output parameter.
   * @param cycle_bonds cyclic bonds output parameter.
   */
  template<typename MoleculeType>
  void cycle_membership(const MoleculeType &mol, std::vector<bool> &cycle_atoms, std::vector<bool> &cycle_bonds)
  {
    cycle_atoms.clear();
    cycle_bonds.clear();
    cycle_atoms.resize(num_atoms(mol));
    cycle_bonds.resize(num_bonds(mol));

    impl::CycleMembershipVisitor<MoleculeType> visitor(cycle_atoms, cycle_bonds);
    depth_first_search(mol, visitor);
  }

  namespace impl {

    template<typename MoleculeType, typename QueryType>
    struct CycleAtomMatcher
    {
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<QueryType>::atom_type query_atom_type;

      CycleAtomMatcher(Index source_) : source(source_)
      {
      }

      bool operator()(const QueryType &query, query_atom_type queryAtom, const MoleculeType &mol, atom_type atom) const
      {
        return get_index(mol, atom) <= source;
      }

      Index source;
    };

    template<typename MoleculeType, typename QueryType>
    struct CycleBondMatcher
    {
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
      typedef typename molecule_traits<QueryType>::bond_type query_bond_type;

      bool operator()(const QueryType &query, query_bond_type queryAtom, const MoleculeType &mol, bond_type bond) const
      {
        return true;
      }
    };

    class IsomorphismCycle
    {
      public:
        IsomorphismCycle()
        {
        }

        template<typename MoleculeType>
        IsomorphismCycle(const MoleculeType &mol, const std::vector<Index> &atoms)
        {
          m_edges.insert(get_index(mol, get_bond(mol, get_atom(mol, atoms[0]), get_atom(mol, atoms[atoms.size() - 1]))));
          for (std::size_t i = 1; i < atoms.size(); ++i)
            m_edges.insert(get_index(mol, get_bond(mol, get_atom(mol, atoms[i - 1]), get_atom(mol, atoms[i]))));
        }

        bool operator==(const IsomorphismCycle &other) const
        {
          return m_edges == other.m_edges;
        }

        IsomorphismCycle& operator+=(const IsomorphismCycle &other)
        {
          std::set<Index> result;
          std::set_symmetric_difference(m_edges.begin(), m_edges.end(),
                                        other.m_edges.begin(), other.m_edges.end(),
                                        std::inserter(result, result.begin()));
          m_edges = result;
          return *this;
        }

        const std::set<Index>& edges() const
        {
          return m_edges;
        }

      private:
        std::set<Index> m_edges;
    };

    bool is_relevant(const std::vector<IsomorphismCycle> &cycles, const IsomorphismCycle &cycle)
    {
      std::size_t n = 0;

      // check for identical
      for (std::size_t i = 0; i < cycles.size(); ++i) {
        if (cycle == cycles[i])
          return false;
        if (cycles[i].edges().size() < cycle.edges().size())
          ++n;
      }

      // check if cycle is sum of smaller cycles
      std::vector<std::size_t> indices;
      for (std::size_t i = 0; i < n; ++i)
        indices.push_back(i);

      for (std::size_t size = 2; size <= indices.size(); ++size)
        do {
          IsomorphismCycle sum;
          for (std::size_t i = 0; i < size; ++i)
            sum += cycles[indices[i]];
          if (cycle == sum)
            return false;
        } while (next_combination(indices.begin(), indices.begin() + size, indices.end()));

      return true;
    }

  }


  /**
   * @brief Find the relevant cycles.
   *
   * The set of relevant cycles is formed by taking the union of all the
   * minimum cycles bases. An alternative definition is that a cycle is
   * relevant if it is not the sum of smaller cycles.
   *
   * This function returns the set of relevant cycles listing the atom
   * indices in sequence.
   *
   * @param mol The molecule.
   * @param cyclomatricNumber The cyclomatic number.
   *
   * @return The set of relevant cycles.
   */
  template<typename MoleculeType>
  RingSet<MoleculeType> relevant_cycles(const MoleculeType &mol, Size cyclomaticNumber)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

    RingSet<MoleculeType> rings(mol);

    // current cycle size
    unsigned int size = 3;
    unsigned int lastSize = 3;

    std::vector<impl::IsomorphismCycle> relevant;
    while (true) {
      std::cout << "size: " << size << ", lastSize: " << lastSize << ", nullity: " << cyclomaticNumber << ", count: " << relevant.size() << std::endl;
      if (rings.size() >= cyclomaticNumber && lastSize < size)
        break;
      // create query
      std::string smiles = "*1" + std::string(size - 1, '*') + "1";
      HeMol cycleMol;
      parse_smiles(smiles, cycleMol);

      // find all cycles of size
      MappingList mappings;
      impl::CycleBondMatcher<MoleculeType, HeMol> bondMatcher;
      FOREACH_ATOM (atom, mol, MoleculeType) {
        impl::CycleAtomMatcher<MoleculeType, HeMol> atomMatcher(get_index(mol, *atom));

        if (isomorphism_search(mol, *atom, cycleMol, mappings, atomMatcher, bondMatcher)) {
          for (std::size_t i = 0; i < mappings.maps.size(); ++i) {
            impl::IsomorphismCycle cycle(mol, mappings.maps[i]);
            if (is_relevant(relevant, cycle)) {
              std::vector<atom_type> atoms;
              for (std::size_t j = 0; j < mappings.maps[i].size(); ++j)
                atoms.push_back(get_atom(mol, mappings.maps[i][j]));
              rings.addRing(Ring<MoleculeType>(mol, atoms));
              relevant.push_back(cycle);
              lastSize = size;
            }
          }
        }
      }

      // increment cycle size
      ++size;
    }

    return rings;
  }

  /**
   * @overload
   */
  template<typename MoleculeType>
  RingSet<MoleculeType> relevant_cycles(const MoleculeType &mol)
  {
    return relevant_cycles(mol, cyclomatic_number(mol));
  }

}

#endif
