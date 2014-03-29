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
#ifndef HELIUM_DIJKSTRA_H
#define HELIUM_DIJKSTRA_H

#include <vector>
#include <limits>
#include <algorithm>

#include <Helium/molecule.h>

namespace Helium {

  namespace impl {

    template<typename MoleculeType>
    struct DijkstraHeapCompare
    {
      DijkstraHeapCompare(const MoleculeType &mol, const std::vector<Size> &dist)
        : m_mol(mol), m_dist(dist)
      {
      }

      bool operator()(typename molecule_traits<MoleculeType>::atom_type first,
                      typename molecule_traits<MoleculeType>::atom_type second) const
      {
        return m_dist[get_index(m_mol, first)] >= m_dist[get_index(m_mol, second)];
      }

      const MoleculeType &m_mol;
      const std::vector<Size> &m_dist;
    };

  }

  /**
   * @brief Class for running Dijkstra's shortest path algorithm.
   *
   * The Dijkstra shortest path algorithm finds a shortest path from a source
   * atom to all other atoms in a molecule. The algorithm is executed when the
   * constructor is executed. Later, the distances and paths can be retrieved
   * using the distance() and path() member functions.
   *
   * @note Complexity: @f$O(n^2)@f$
   * @ingroup Beta
   * @note Phase: Beta
   */
  template<typename MoleculeType>
  class Dijkstra
  {
    public:
      /**
       * @brief Constructor.
       *
       * Using this constructor, all atoms will be considered.
       *
       * @param mol The molecule.
       * @param source The source atom.
       * @param preferSmallerIndexPaths When true, same length paths are
       *        preferred that include only atoms with indices smaller than the
       *        source's index.
       */
      template<typename AtomType>
      Dijkstra(const MoleculeType &mol, AtomType source, bool preferSmallerIndexPaths = false)
          : m_mol(mol), m_source(source)
      {
        // add all atoms in mol to Q
        std::vector<AtomType> Q;
        FOREACH_ATOM_T (atom, mol, MoleculeType)
          Q.push_back(*atom);

        dijkstra(mol, source, Q, preferSmallerIndexPaths);
      }

      /**
       * @brief Get the infinity value that is returned by distance().
       *
       * @return The infinity value.
       */
      Size infinity() const
      {
        return std::numeric_limits<Size>::max();
      }

      /**
       * @brief Get the distance between source and target atoms.
       *
       * The distance is the number of bonds between the source and target
       * atoms. If there is no path, the infinity() value is returned.
       *
       * @param target the target atom.
       */
      template<typename AtomType>
      Size distance(AtomType target) const
      {
        return m_dist[get_index(m_mol, target)];
      }

      /**
       * @brief Reconstruct the path between source and target atoms.
       *
       * The path includes the source and target atoms (i.e. [source, ..., target]).
       *
       * @param target the target atom.
       */
      template<typename AtomType>
      std::vector<AtomType> path(AtomType target) const
      {
        std::vector<AtomType> S;
        AtomType u = target;

        while (m_prev[get_index(m_mol, u)] != molecule_traits<MoleculeType>::null_atom()) {
          S.insert(S.begin(), u);
          u = m_prev[get_index(m_mol, u)];
        }

        S.insert(S.begin(), m_source);

        return S;
      }

    private:
      template<typename AtomType>
      void dijkstra(const MoleculeType &mol, AtomType source, std::vector<AtomType> &Q, bool preferSmallerIndexPaths)
      {
        // distance from source to other atoms
        m_dist.resize(num_atoms(mol), std::numeric_limits<Size>::max());
        // previous atom in shortest path from source
        m_prev.resize(num_atoms(mol), molecule_traits<MoleculeType>::null_atom());

        // distance from source to source
        m_dist[get_index(mol, source)] = 0;

        // make Q a heap [heap implementation only]
        typedef typename impl::DijkstraHeapCompare<MoleculeType> HeapCompare;
        std::make_heap(Q.begin(), Q.end(), HeapCompare(mol, m_dist));

        while (Q.size()) {
          /* [array implementation]
          // find atom in Q with smallest distance in m_dist
          std::size_t uIndex = 0;
          for (std::size_t i = 1; i < Q.size(); ++i)
            if (m_dist[get_index(mol, Q[i])] < m_dist[get_index(mol, Q[uIndex])])
              uIndex = i;
          AtomType u = Q[uIndex];
          // remove u from Q
          Q.erase(Q.begin() + uIndex);
          */

          // find atom in Q with smallest distance in m_dist
          AtomType u = Q.front();
          // remove u from Q
          std::pop_heap(Q.begin(), Q.end(), HeapCompare(mol, m_dist));
          Q.pop_back();

          if (m_dist[get_index(mol, u)] == std::numeric_limits<Size>::max())
            // all remaining atoms are inaccessible from source
            break;

          FOREACH_NBR_T (v, u, mol, MoleculeType) {
            Size alt = m_dist[get_index(mol, u)] + 1;

            if (alt < m_dist[get_index(mol, *v)]) {
              m_dist[get_index(mol, *v)] = alt;
              m_prev[get_index(mol, *v)] = u;
              // [heap implementation only]
              std::make_heap(Q.begin(), Q.end(), HeapCompare(mol, m_dist));
            }

            if (preferSmallerIndexPaths && alt == m_dist[get_index(mol, *v)] &&
                contains_larger_index(mol, path(*v), get_index(mol, m_source)))
              m_prev[get_index(mol, *v)] = u;
          }
        }
      }

      const MoleculeType &m_mol;
      typename molecule_traits<MoleculeType>::atom_type m_source;
      std::vector<Size> m_dist;
      std::vector<typename molecule_traits<MoleculeType>::atom_type> m_prev;
  };

}

#endif
