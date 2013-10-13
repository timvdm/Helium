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
#ifndef HELIUM_ISOMORPHISM_H
#define HELIUM_ISOMORPHISM_H

#include <Helium/molecule.h>
#include <Helium/tie.h>

#include <vector>
#include <cassert>
#include <algorithm>
#include <set>

#define DEBUG_ISOMORPHISM 0
#include <iostream>

namespace Helium {

  /**
   * Type used to represent a single isomorphism mapping.
   */
  typedef std::vector<Index> IsomorphismMapping;
  /**
   * Type used to represent multiple isomorphism mappings.
   */
  typedef std::vector<IsomorphismMapping> IsomorphismMappings;

  /**
   * NoMapping
   */
  struct NoMapping
  {
    enum { single = true };

    /**
     * @brief Constructor.
     */
    NoMapping() : match(false)
    {
    }

    /**
     * @brief True when a match is found.
     */
    bool match;
  };

  /**
   * CountMapping
   */
  struct CountMapping
  {
    enum { single = false };

    /**
     * @brief Constructor.
     */
    CountMapping() : count(0)
    {
    }

    /**
     * @brief The number of found matches.
     */
    int count;
  };

  /**
   * SingleMapping
   */
  struct SingleMapping
  {
    enum { single = true };

    /**
     * @brief The atom mapping.
     */
    IsomorphismMapping map;
  };

  /**
   * MappingList
   */
  struct MappingList
  {
    enum { single = false };

    /**
     * @brief The atom mappings.
     */
    IsomorphismMappings maps;
  };

  namespace impl {

    /**
     * Clear mapping implementations.
     */
    inline void clear_mappig(NoMapping &mapping)
    {
      mapping.match = false;
    }
    inline void clear_mappig(CountMapping &mapping)
    {
      mapping.count = 0;
    }
    inline void clear_mappig(SingleMapping &mapping)
    {
      mapping.map.clear();
    }
    inline void clear_mappig(MappingList &mapping)
    {
      mapping.maps.clear();
    }

    /**
     * Add mapping implementations.
     */
    inline void add_mapping(NoMapping &mapping, const IsomorphismMapping &map)
    {
      mapping.match = true;
    }
    inline void add_mapping(CountMapping &mapping, const IsomorphismMapping &map)
    {
      mapping.count++;
    }
    inline void add_mapping(SingleMapping &mapping, const IsomorphismMapping &map)
    {
      mapping.map = map;
    }
    inline void add_mapping(MappingList &mapping, const IsomorphismMapping &map)
    {
      mapping.maps.push_back(map);
    }

    /**
     * Empty mapping implementations.
     */
    inline bool empty_mappig(NoMapping &mapping)
    {
      return !mapping.match;
    }
    inline bool empty_mappig(CountMapping &mapping)
    {
      return mapping.count == 0;
    }
    inline bool empty_mappig(SingleMapping &mapping)
    {
      return mapping.map.empty();
    }
    inline bool empty_mappig(MappingList &mapping)
    {
      return mapping.maps.empty();
    }

    inline void print_map(const IsomorphismMapping &map)
    {
      std::cout << "m_map: ";
      for (std::size_t i = 0; i < map.size(); ++i)
        std::cout << map[i] << " ";
      std::cout << std::endl;
    }

    template<typename MoleculeType, typename QueryType, typename MappingType, typename AtomMatcher, typename BondMatcher>
    class Isomorphism
    {
      public:
        typedef typename molecule_traits<QueryType>::atom_type query_atom_type;
        typedef typename molecule_traits<QueryType>::bond_type query_bond_type;
        typedef typename molecule_traits<QueryType>::incident_iter query_incident_iter;

        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;
        typedef typename molecule_traits<MoleculeType>::incident_iter incident_iter;

        Isomorphism(MoleculeType &mol, QueryType &query,
            const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher)
          : m_atomMatcher(atomMatcher), m_bondMatcher(bondMatcher), m_mol(mol), m_query(query)
        {
          dfsBonds();
          m_map.resize(num_atoms(query), -1);
        }

        void dfsBonds(query_atom_type atom, std::vector<bool> &visited)
        {
          query_incident_iter bond, end_bonds;
          TIE(bond, end_bonds) = get_bonds(m_query, atom);
          for (; bond != end_bonds; ++bond) {
            if (visited[get_index(m_query, *bond)])
              continue;
            visited[get_index(m_query, *bond)] = true;

            bool swap = get_source(m_query, *bond) != atom;
            m_bondSwap.push_back(swap);

            m_dfsBonds.push_back(get_index(m_query, *bond));
            dfsBonds(swap ? get_source(m_query, *bond) : get_target(m_query, *bond), visited);
          }
        }

        void dfsBonds()
        {
          std::vector<bool> visited(num_bonds(m_query));

          // TODO extend for disconnected...
          query_atom_type atom = get_atom(m_query, 0);
          dfsBonds(atom, visited);

          if (DEBUG_ISOMORPHISM) {
            std::cout << "DFS bonds: ";
            for (std::size_t i = 0; i < m_dfsBonds.size(); ++i) {
              query_bond_type bond = get_bond(m_query, m_dfsBonds[i]);
              if (!m_bondSwap[i])
                std::cout << get_index(m_query, get_source(m_query, bond)) << "-" << get_index(m_query, get_target(m_query, bond)) << " ";
              else
                std::cout << get_index(m_query, get_target(m_query, bond)) << "-" << get_index(m_query, get_source(m_query, bond)) << " ";
            }
            std::cout << std::endl;
          }
        }

        void match(MappingType &mapping, int bondIndex)
        {
          query_bond_type queryBond = get_bond(m_query, m_dfsBonds[bondIndex]);
          query_atom_type querySource = m_bondSwap[bondIndex] ? get_target(m_query, queryBond) : get_source(m_query, queryBond);

          atom_type atom = get_atom(m_mol, m_map[get_index(m_mol, querySource)]);

          query_atom_type queryTarget = m_bondSwap[bondIndex] ? get_source(m_query, queryBond) : get_target(m_query, queryBond);
          bool isRingClosure = m_map[get_index(m_query, queryTarget)] != -1;

          if (DEBUG_ISOMORPHISM)
            std::cout << "mapping bond: " << get_index(m_query, querySource) << "-" << get_index(m_query, queryTarget) << std::endl;

          incident_iter bond, end_bonds;
          TIE(bond, end_bonds) = get_bonds(m_mol, atom);
          for (; bond != end_bonds; ++bond) {
            if (!m_bondMatcher(m_query, queryBond, m_mol, *bond))
              continue;

            atom_type nbrAtom = get_other(m_mol, *bond, atom);

            if (isRingClosure) {
              if (DEBUG_ISOMORPHISM)
                std::cout << "ring closure..." << std::endl;
              if (m_map[get_index(m_query, queryTarget)] != get_index(m_mol, nbrAtom))
                continue;
            } else {
              if (std::find(m_map.begin(), m_map.end(), get_index(m_mol, nbrAtom)) != m_map.end())
                continue;

              if (!m_atomMatcher(m_query, queryTarget, m_mol, nbrAtom))
                continue;

              if (DEBUG_ISOMORPHISM)
                std::cout << get_index(m_query, queryTarget) << " -> " << get_index(m_mol, nbrAtom) << std::endl;

              // map the target atom
              m_map[get_index(m_query, queryTarget)] = get_index(m_mol, nbrAtom);
              //print_map(m_map);
            }

            if (bondIndex + 1 == m_dfsBonds.size()) {
              if (DEBUG_ISOMORPHISM) {
                std::cout << "found mapping..." << std::endl;
                print_map(m_map);
              }

              // create bit mask of atoms (to ensure uniqueness of mapping)
              std::vector<bool> atoms(num_atoms(m_mol));
              for (std::size_t i = 0; i < m_map.size(); ++i)
                atoms[m_map[i]] = true;

              // add the mapping to the result if it is unique
              if (m_mappings.find(atoms) == m_mappings.end()) {
                m_mappings.insert(atoms);
                add_mapping(mapping, m_map);
              }
            } else
              match(mapping, bondIndex + 1);

            // bracktrack target atom
            if (!isRingClosure) {
              if (DEBUG_ISOMORPHISM)
                std::cout << "backtrack: " << get_index(m_query, queryTarget) << std::endl;

              m_map[get_index(m_query, queryTarget)] = -1;
            }

            // exit as soon as possible if only one match is required
            if (MappingType::single && !empty_mappig(mapping))
              return;
          }


        }

        void match(MappingType &mapping, atom_type atom)
        {
          if (!num_atoms(m_query))
            return;

          if (!num_bonds(m_query)) {
            query_atom_type queryAtom = get_atom(m_query, 0);

            if (!m_atomMatcher(m_query, queryAtom, m_mol, atom))
              return;

            if (DEBUG_ISOMORPHISM)
              std::cout << get_index(m_query, queryAtom) << " -> " << get_index(m_mol, atom) << std::endl;

            // map the source atom
            m_map[get_index(m_query, queryAtom)] = get_index(m_mol, atom);

            add_mapping(mapping, m_map);

            if (DEBUG_ISOMORPHISM)
              std::cout << "backtrack: " << get_index(m_query, queryAtom) << std::endl;
            m_map[get_index(m_query, queryAtom)] = -1;

            if (MappingType::single && !empty_mappig(mapping))
              return;
          } else {
            query_atom_type queryAtom = m_bondSwap[0] ? get_target(m_query, get_bond(m_query, m_dfsBonds[0])) : get_source(m_query, get_bond(m_query, m_dfsBonds[0]));

            if (!m_atomMatcher(m_query, queryAtom, m_mol, atom))
              return;

            if (DEBUG_ISOMORPHISM)
              std::cout << get_index(m_query, queryAtom) << " -> " << get_index(m_mol, atom) << std::endl;

            // map the source atom
            m_map[get_index(m_query, queryAtom)] = get_index(m_mol, atom);

            match(mapping, 0);

            if (DEBUG_ISOMORPHISM)
              std::cout << "backtrack: " << get_index(m_query, queryAtom) << std::endl;
            m_map[get_index(m_query, queryAtom)] = -1;

            if (MappingType::single && !empty_mappig(mapping))
              return;
          }
        }

        void match(MappingType &mapping)
        {
          if (!num_atoms(m_query))
            return;

          if (!num_bonds(m_query)) {
            atom_iter atom, end_atoms;
            TIE(atom, end_atoms) = get_atoms(m_mol);
            for (; atom != end_atoms; ++atom)
              match(mapping, *atom);
          } else {
            // try to match each atom in the molecule against the first atom
            // epxression in the SMARTS
            atom_iter atom, end_atoms;
            TIE(atom, end_atoms) = get_atoms(m_mol);
            for (; atom != end_atoms; ++atom)
              match(mapping, *atom);
          }
        }

      private:
        const AtomMatcher &m_atomMatcher;
        const BondMatcher &m_bondMatcher;
        MoleculeType &m_mol; // the queried molecule
        QueryType &m_query; // the query
        IsomorphismMapping  m_map; // current mapping: query atom index -> queried atom index
        std::set<std::vector<bool> > m_mappings; // keep track of unique mappins
        std::vector<Index> m_dfsBonds; // dfs bond order
        std::vector<bool> m_bondSwap; // should source/target atoms be swapped
    };

  }

  /**
   * @brief The default atom matcher for isomorphism searches.
   */
  template<typename MoleculeType, typename QueryType>
  struct DefaultAtomMatcher
  {
    /**
     * @brief The molecule atom type.
     */
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    /**
     * @brief The query atom type.
     */
    typedef typename molecule_traits<QueryType>::atom_type query_atom_type;

    /**
     * @brief Atom match function.
     *
     * @brief query The query.
     * @brief queryAtom The query atom.
     * @brief mol The queried molecule.
     * @brief atom The queried atom.
     *
     * @return True if the query atom matches the molecule atom.
     */
    bool operator()(QueryType &query, query_atom_type queryAtom, MoleculeType &mol, atom_type atom) const
    {
      return get_element(query, queryAtom) == get_element(mol, atom);
    }
  };

  /**
   * @brief The default bond matcher for isomorphism searches.
   */
  template<typename MoleculeType, typename QueryType>
  struct DefaultBondMatcher
  {
    /**
     * @brief The molecule bond type.
     */
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
    /**
     * @brief The query bond type.
     */
    typedef typename molecule_traits<QueryType>::bond_type query_bond_type;

    /**
     * @brief Bond match function.
     *
     * @brief query The query.
     * @brief queryBond The query bond.
     * @brief mol The queried molecule.
     * @brief bond The queried bond.
     *
     * @return True if the query atom matches the molecule atom.
     */
    bool operator()(QueryType &query, query_bond_type queryBond, MoleculeType &mol, bond_type bond) const
    {
      return get_order(query, queryBond) == get_order(mol, bond);
    }
  };


  /**
   * Perform a subgraph isomorphism search for the specified query in the
   * molecule.
   *
   * @param mol The molecule (queried).
   * @param atom The first atom to be matched against query atom 0.
   * @param query The query.
   * @param mapping The desired mapping (e.g. NoMapping, SingleMapping, ...).
   *
   * @return True if the query is a substructure of @p mol.
   */
  template<typename MoleculeType, typename AtomType, typename QueryType, typename MappingType, typename AtomMatcher, typename BondMatcher>
  bool isomorphism_search(MoleculeType &mol, AtomType atom, QueryType &query, MappingType &mapping,
      const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher)
  {
    impl::clear_mappig(mapping);

    if (!num_atoms(query))
      return false;

    impl::Isomorphism<MoleculeType, QueryType, MappingType, AtomMatcher, BondMatcher> iso(mol, query, atomMatcher, bondMatcher);
    iso.match(mapping, atom);

    return !impl::empty_mappig(mapping);
  }

  /**
   * Perform a subgraph isomorphism search for the specified query in the
   * molecule.
   *
   * @param mol The molecule (queried).
   * @param query The query.
   * @param mapping The desired mapping (e.g. NoMapping, SingleMapping, ...).
   *
   * @return True if the query is a substructure of @p mol.
   */
  template<typename MoleculeType, typename QueryType, typename MappingType, typename AtomMatcher, typename BondMatcher>
  bool isomorphism_search(MoleculeType &mol, QueryType &query, MappingType &mapping,
      const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher)
  {
    impl::clear_mappig(mapping);

    if (!num_atoms(query))
      return false;

    impl::Isomorphism<MoleculeType, QueryType, MappingType, AtomMatcher, BondMatcher> iso(mol, query, atomMatcher, bondMatcher);
    iso.match(mapping);

    return !impl::empty_mappig(mapping);
  }

  /**
   * @overload
   */
  template<typename MoleculeType, typename QueryType, typename AtomMatcher, typename BondMatcher>
  bool isomorphism_search(MoleculeType &mol, QueryType &query,
      const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher)
  {
    NoMapping mapping;
    return isomorphism_search(mol, query, mapping, atomMatcher, bondMatcher);
  }

}

#endif
