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
#ifndef HELIUM_ISOMORPHISM_H
#define HELIUM_ISOMORPHISM_H

#include <Helium/molecule.h>
#include <Helium/stereo.h>

#include <vector>
#include <cassert>
#include <algorithm>
#include <set>

//#define DEBUG_ISOMORPHISM 0
#include <iostream>

namespace Helium {

  /**
   * @file algorithms/isomorphism.h
   * @brief Subgraph isomorphism algorithm.
   */

  /**
   * Type used to represent a single isomorphism mapping.
   */
  typedef std::vector<Index> IsomorphismMapping;

  /**
   * Type used to represent multiple isomorphism mappings.
   */
  typedef std::vector<IsomorphismMapping> IsomorphismMappings;

  /**
   * @struct NoMapping algorithms/isomorphism.h <Helium/algorithms/isomorphism.h>
   * @brief Do not map an isomorphism match.
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
   * @struct CountMapping algorithms/isomorphism.h <Helium/algorithms/isomorphism.h>
   * @brief Count the number of isomorphism mappings.
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
   * @struct SingleMapping algorithms/isomorphism.h <Helium/algorithms/isomorphism.h>
   * @brief Record only a single isomorphism mapping.
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
   * @struct MappingList algorithms/isomorphism.h <Helium/algorithms/isomorphism.h>
   * @brief Record all isomorphism mappings.
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

    // defined in stereo.cpp
    bool stereo_compare_undirected_cycles(const Stereo::Ref *refs, Stereo::Ref ref1, Stereo::Ref ref2,
        Stereo::Ref ref3, Stereo::Ref ref4);

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

        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;

        Isomorphism(const MoleculeType &mol, const QueryType &query,
            const Stereochemistry &molStereo, const Stereochemistry &queryStereo,
            const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher)
          : m_atomMatcher(atomMatcher), m_bondMatcher(bondMatcher), m_mol(mol), m_query(query),
            m_molStereo(molStereo), m_queryStereo(queryStereo)
        {
          dfsBonds();
          m_map.resize(num_atoms(query), -1);
        }

        void dfsBonds(query_atom_type atom, std::vector<bool> &visited)
        {
          for (auto &bond : get_bonds(m_query, atom)) {
            if (visited[get_index(m_query, bond)])
              continue;
            visited[get_index(m_query, bond)] = true;

            bool swap = get_source(m_query, bond) != atom;
            m_bondSwap.push_back(swap);

            m_dfsBonds.push_back(get_index(m_query, bond));
            dfsBonds(swap ? get_source(m_query, bond) : get_target(m_query, bond), visited);
          }
        }

        void dfsBonds()
        {
          std::vector<bool> visited(num_bonds(m_query));

          // TODO extend for disconnected...
          query_atom_type atom = get_atom(m_query, 0);
          dfsBonds(atom, visited);

          /*
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
          */
        }

        bool stereoMatches()
        {
          //std::cout << "stereoMatches()" << std::endl;
          //std::cout << "  m_map: " << m_map << std::endl;
          for (auto &queryStereo : m_queryStereo.allStereo()) {
            if (queryStereo.type() == Stereo::None)
              continue;

            //std::cout << "  queryStereo: " << queryStereo << std::endl;

            Stereo::Ref molCenter;
            if (queryStereo.type() == Stereo::CisTrans) {
              auto queryBond = get_bond(m_query, queryStereo.center());
              auto querySource = get_source(m_query, queryBond);
              auto queryTarget = get_target(m_query, queryBond);

              auto molSource = get_atom(m_mol, m_map[get_index(m_query, querySource)]);
              auto molTarget = get_atom(m_mol, m_map[get_index(m_query, queryTarget)]);
              auto molBond = get_bond(m_mol, molSource, molTarget);
              assert(molBond != molecule_traits<MoleculeType>::null_bond());

              molCenter = get_index(m_mol, molBond);
            } else
              molCenter = m_map[queryStereo.center()];

            //std::cout << "  molCenter: " << molCenter << std::endl;

            // find matching stereo in mol
            bool found = false;
            for (auto &molStereo : m_molStereo.allStereo()) {
              if (queryStereo.type() != molStereo.type())
                continue;
              if (molCenter != molStereo.center())
                continue;

              //std::cout << "  molStereo: " << molStereo << std::endl;

              int numRefs = queryStereo.numRefs();

              // map the query refs using the found mapping
              Stereo::Ref queryRefs[6] = {Stereo::implRef(), Stereo::implRef(), Stereo::implRef(),
                Stereo::implRef(), Stereo::implRef(), Stereo::implRef()};
              for (int i = 0; i < numRefs; ++i)
                if (queryStereo.ref(i) != Stereo::implRef())
                  if (!is_hydrogen(m_query, get_atom(m_query, queryStereo.ref(i))))
                    queryRefs[i] = m_map[queryStereo.ref(i)];

              // make hydrogens implicit in molStereo
              Stereo::Ref molRefs[6] = {Stereo::implRef(), Stereo::implRef(), Stereo::implRef(),
                Stereo::implRef(), Stereo::implRef(), Stereo::implRef()};
              for (int i = 0; i < numRefs; ++i)
                if (molStereo.ref(i) != Stereo::implRef())
                  if (!is_hydrogen(m_mol, get_atom(m_mol, molStereo.ref(i))))
                    molRefs[i] = molStereo.ref(i);

              StereoStorage storage(molStereo.type(), molStereo.center(), molRefs, molRefs + numRefs);
              //std::cout << "  storage: " << storage << std::endl;

              switch (queryStereo.type()) {
                case Stereo::Tetrahedral:
                case Stereo::Allene:
                  found = (Stereo::TH1 == tetrahedral_class(storage, queryRefs[0], queryRefs[1],
                        queryRefs[2], queryRefs[3]));
                  break;
                case Stereo::SquarePlanar:
                  found = (Stereo::SP1 == squareplanar_class(storage, queryRefs[0], queryRefs[1],
                        queryRefs[2], queryRefs[3]));
                  break;
                case Stereo::TrigonalBipyramidal:
                  found = (Stereo::TB1 == trigonalbipyramidal_class(storage, queryRefs[0], queryRefs[1],
                        queryRefs[2], queryRefs[3], queryRefs[4]));
                  break;
                case Stereo::Octahedral:
                  found = (Stereo::OH1 == octahedral_class(storage, queryRefs[0], queryRefs[1], queryRefs[2],
                        queryRefs[3], queryRefs[4], queryRefs[5]));
                  break;
                case Stereo::CisTrans:
                  found = stereo_compare_undirected_cycles(queryRefs, molRefs[0], molRefs[1], molRefs[2], molRefs[3]);
                  break;
                default:
                  break;
              }

              if (found)
                break;
            }

            // if there is no stereo in mol, return false
            if (!found)
              return false;
          }

          return true;
        }

        void match(MappingType &mapping, int bondIndex)
        {
          query_bond_type queryBond = get_bond(m_query, m_dfsBonds[bondIndex]);
          query_atom_type querySource = m_bondSwap[bondIndex] ? get_target(m_query, queryBond) : get_source(m_query, queryBond);

          atom_type atom = get_atom(m_mol, m_map[get_index(m_query, querySource)]);

          query_atom_type queryTarget = m_bondSwap[bondIndex] ? get_source(m_query, queryBond) : get_target(m_query, queryBond);
          bool isRingClosure = m_map[get_index(m_query, queryTarget)] != -1;

          /*
          if (DEBUG_ISOMORPHISM)
            std::cout << "mapping bond: " << get_index(m_query, querySource) << "-" << get_index(m_query, queryTarget) << std::endl;
          */

          for (auto &bond : get_bonds(m_mol, atom)) {
            if (!m_bondMatcher(m_query, queryBond, m_mol, bond))
              continue;

            atom_type nbrAtom = get_other(m_mol, bond, atom);

            if (isRingClosure) {
              /*
              if (DEBUG_ISOMORPHISM)
                std::cout << "ring closure..." << std::endl;
              */
              if (m_map[get_index(m_query, queryTarget)] != get_index(m_mol, nbrAtom))
                continue;
            } else {
              if (std::find(m_map.begin(), m_map.end(), get_index(m_mol, nbrAtom)) != m_map.end())
                continue;

              if (!m_atomMatcher(m_query, queryTarget, m_mol, nbrAtom))
                continue;

              /*
              if (DEBUG_ISOMORPHISM)
                std::cout << get_index(m_query, queryTarget) << " -> " << get_index(m_mol, nbrAtom) << std::endl;
              */

              // map the target atom
              m_map[get_index(m_query, queryTarget)] = get_index(m_mol, nbrAtom);
              //print_map(m_map);
            }

            if (bondIndex + 1 == m_dfsBonds.size()) {
              /*
              if (DEBUG_ISOMORPHISM) {
                std::cout << "found mapping..." << std::endl;
                print_map(m_map);
              }
              */

              // create bit mask of atoms (to ensure uniqueness of mapping)
              std::vector<bool> atoms(num_atoms(m_mol));
              for (std::size_t i = 0; i < m_map.size(); ++i)
                atoms[m_map[i]] = true;

              // add the mapping to the result if it is unique
              if (m_mappings.find(atoms) == m_mappings.end()) {
                m_mappings.insert(atoms);
                if (stereoMatches())
                  add_mapping(mapping, m_map);
              }
            } else
              match(mapping, bondIndex + 1);

            // bracktrack target atom
            if (!isRingClosure) {
              /*
              if (DEBUG_ISOMORPHISM)
                std::cout << "backtrack: " << get_index(m_query, queryTarget) << std::endl;
              */

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

            /*
            if (DEBUG_ISOMORPHISM)
              std::cout << get_index(m_query, queryAtom) << " -> " << get_index(m_mol, atom) << std::endl;
            */

            // map the source atom
            m_map[get_index(m_query, queryAtom)] = get_index(m_mol, atom);

            add_mapping(mapping, m_map);

            /*
            if (DEBUG_ISOMORPHISM)
              std::cout << "backtrack: " << get_index(m_query, queryAtom) << std::endl;
            */
            m_map[get_index(m_query, queryAtom)] = -1;

            if (MappingType::single && !empty_mappig(mapping))
              return;
          } else {
            query_atom_type queryAtom = m_bondSwap[0] ? get_target(m_query, get_bond(m_query, m_dfsBonds[0])) : get_source(m_query, get_bond(m_query, m_dfsBonds[0]));

            if (!m_atomMatcher(m_query, queryAtom, m_mol, atom))
              return;

            /*
            if (DEBUG_ISOMORPHISM)
              std::cout << get_index(m_query, queryAtom) << " -> " << get_index(m_mol, atom) << std::endl;
            */

            // map the source atom
            m_map[get_index(m_query, queryAtom)] = get_index(m_mol, atom);

            match(mapping, 0);

            /*
            if (DEBUG_ISOMORPHISM)
              std::cout << "backtrack: " << get_index(m_query, queryAtom) << std::endl;
            */
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
            for (auto &atom : get_atoms(m_mol))
              match(mapping, atom);
          } else {
            // try to match each atom in the molecule against the first atom
            // epxression in the SMARTS
            for (auto &atom : get_atoms(m_mol))
              match(mapping, atom);
          }
        }

      private:
        const AtomMatcher &m_atomMatcher; // the atom matcher
        const BondMatcher &m_bondMatcher; // the bond matcher
        const MoleculeType &m_mol; // the queried molecule
        const QueryType &m_query; // the query
        const Stereochemistry &m_molStereo;
        const Stereochemistry &m_queryStereo;
        IsomorphismMapping  m_map; // current mapping: query atom index -> queried atom index
        std::set<std::vector<bool> > m_mappings; // keep track of unique mappins
        std::vector<Index> m_dfsBonds; // dfs bond order
        std::vector<bool> m_bondSwap; // should source/target atoms be swapped
    };

  }

  /**
   * @struct DefaultAtomMatcher algorithms/isomorphism.h <Helium/algorithms/isomorphism.h>
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
    bool operator()(const QueryType &query, query_atom_type queryAtom,
        const MoleculeType &mol, atom_type atom) const
    {
      return get_element(query, queryAtom) == get_element(mol, atom);
    }
  };

  /**
   * @struct DefaultBondMatcher algorithms/isomorphism.h <Helium/algorithms/isomorphism.h>
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
    bool operator()(const QueryType &query, query_bond_type queryBond,
        const MoleculeType &mol, bond_type bond) const
    {
      return get_order(query, queryBond) == get_order(mol, bond);
    }
  };


  /**
   * Perform a subgraph isomorphism search for the specified query in the
   * molecule.
   *
   * @note This function only works for single component queries.
   *
   * @param mol The molecule (queried).
   * @param atom The first atom to be matched against query atom 0.
   * @param query The query.
   * @param molStereo The molecule stereochemistry.
   * @param queryStereo The query stereochemistry.
   * @param mapping The desired mapping (e.g. NoMapping, SingleMapping, ...).
   * @param atomMatcher The AtomMatcher functor.
   * @param bondMatcher The BondMatcher functor.
   *
   * @return True if the query is a substructure of @p mol.
   */
  template<typename MoleculeType, typename QueryType, typename MappingType, typename AtomMatcher, typename BondMatcher>
  bool isomorphism_search(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom, QueryType &query,
      const Stereochemistry &molStereo, const Stereochemistry &queryStereo, MappingType &mapping,
      const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher)
  {
    impl::clear_mappig(mapping);

    if (!num_atoms(query))
      return false;

    impl::Isomorphism<MoleculeType, QueryType, MappingType, AtomMatcher, BondMatcher>
      iso(mol, query, molStereo, queryStereo, atomMatcher, bondMatcher);
    iso.match(mapping, atom);

    return !impl::empty_mappig(mapping);
  }

  /**
   * @overload
   */
  template<typename MoleculeType, typename QueryType, typename MappingType, typename AtomMatcher, typename BondMatcher>
  bool isomorphism_search(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom, QueryType &query, MappingType &mapping,
      const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher)
  {
    Stereochemistry molStereo, queryStereo;
    return isomorphism_search(mol, atom, query, molStereo, queryStereo, mapping, atomMatcher, bondMatcher);
  }

  /**
   * Perform a subgraph isomorphism search for the specified query in the
   * molecule.
   *
   * @note This function only works for single component queries.
   *
   * @param mol The molecule (queried).
   * @param query The query.
   * @param molStereo The molecule stereochemistry.
   * @param queryStereo The query stereochemistry.
   * @param mapping The desired mapping (e.g. NoMapping, SingleMapping, ...).
   * @param atomMatcher The AtomMatcher functor.
   * @param bondMatcher The BondMatcher functor.
   *
   * @return True if the query is a substructure of @p mol.
   */
  template<typename MoleculeType, typename QueryType, typename MappingType, typename AtomMatcher, typename BondMatcher>
  bool isomorphism_search(const MoleculeType &mol, const QueryType &query,
      const Stereochemistry &molStereo, const Stereochemistry &queryStereo, MappingType &mapping,
      const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher)
  {
    impl::clear_mappig(mapping);

    if (!num_atoms(query))
      return false;

    impl::Isomorphism<MoleculeType, QueryType, MappingType, AtomMatcher, BondMatcher>
      iso(mol, query, molStereo, queryStereo, atomMatcher, bondMatcher);
    iso.match(mapping);

    return !impl::empty_mappig(mapping);
  }

  /**
   * @overload
   */
  template<typename MoleculeType, typename QueryType, typename MappingType, typename AtomMatcher, typename BondMatcher>
  bool isomorphism_search(const MoleculeType &mol, const QueryType &query, MappingType &mapping,
      const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher)
  {
    Stereochemistry molStereo, queryStereo;
    return isomorphism_search(mol, query, molStereo, queryStereo, mapping, atomMatcher, bondMatcher);
  }

  /**
   * @overload
   */
  template<typename MoleculeType, typename QueryType, typename AtomMatcher, typename BondMatcher>
  bool isomorphism_search(const MoleculeType &mol, const QueryType &query,
      const Stereochemistry &molStereo, const Stereochemistry &queryStereo,
      const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher)
  {
    NoMapping mapping;
    return isomorphism_search(mol, query, molStereo, queryStereo, mapping, atomMatcher, bondMatcher);
  }

  /**
   * @overload
   */
  template<typename MoleculeType, typename QueryType, typename AtomMatcher, typename BondMatcher>
  bool isomorphism_search(const MoleculeType &mol, const QueryType &query,
      const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher)
  {
    NoMapping mapping;
    return isomorphism_search(mol, query, mapping, atomMatcher, bondMatcher);
  }

}

#endif
