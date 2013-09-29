#ifndef HELIUM_DIJKSTRA_H
#define HELIUM_DIJKSTRA_H

#include <vector>
#include <limits>

#include <Helium/molecule.h>

namespace Helium {

  struct DijkstraSmallerIndexAtoms
  {
    DijkstraSmallerIndexAtoms(Index index) : m_index(index)
    {
    }

    template<typename MoleculeType, typename AtomType>
    bool operator()(const MoleculeType &mol, AtomType atom) const
    {
      return get_index(mol, atom) <= m_index;
    }

    Index m_index;
  };

  /**
   * @brief Class for running Dijkstra's shortest path algorithm.
   *
   * The Dijkstra shortest path algorithm finds a shortest path from a source
   * atom to all other atoms in a molecule. The algorithm is executed when the
   * constructor is executed. Later, the distances and paths can be retrieved
   * using the distance() and path() member functions.
   *
   * @note Complexity: O(n^2)
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
       */
      template<typename AtomType>
      Dijkstra(const MoleculeType &mol, AtomType source)
          : m_mol(mol), m_source(source)
      {
        // add all atoms in mol to Q
        std::vector<AtomType> Q;
        FOREACH_ATOM (atom, mol, MoleculeType)
          Q.push_back(*atom);

        dijkstra(mol, source, Q);
      }

      /**
       * @brief Constructor.
       *
       * Using this constructor, some atoms can be excluded using the
       * appropriate predicate functor. The functor should implement the
       * function call operator as illustrated below:
       *
       * @code
       * struct MyDijkstraAtomPredicate
       * {
       *   template<typename MoleculeType, typename AtomType>
       *   bool operator()(const MoleculeType &mol, AtomType atom) const
       *   {
       *     return true; // include all atoms..
       *   }
       * };
       * @endcode
       *
       * @param mol The molecule.
       * @param source The source atom.
       * @param predicate The atoms predicate.
       */
      template<typename AtomType, typename AtomPredicate>
      Dijkstra(const MoleculeType &mol, AtomType source, const AtomPredicate &predicate)
          : m_mol(mol), m_source(source)
      {
        // add all atoms in mol to Q
        std::vector<AtomType> Q;
        FOREACH_ATOM (atom, mol, MoleculeType)
          if (predicate(mol, *atom))
            Q.push_back(*atom);

        dijkstra(mol, source, Q);
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
      void dijkstra(const MoleculeType &mol, AtomType source, std::vector<AtomType> &Q)
      {
        // distance from source to other atoms
        m_dist.resize(num_atoms(mol), std::numeric_limits<Size>::max());
        // previous atom in shortest path from source
        m_prev.resize(num_atoms(mol), molecule_traits<MoleculeType>::null_atom());

        // distance from source to source
        m_dist[get_index(mol, source)] = 0;

        while (Q.size()) {
          // find atom in Q with smallest distance in m_dist
          std::size_t uIndex = 0;
          for (std::size_t i = 1; i < Q.size(); ++i)
            if (m_dist[get_index(mol, Q[i])] < m_dist[get_index(mol, Q[uIndex])])
              uIndex = i;
          AtomType u = Q[uIndex];
          // remove u from Q
          Q.erase(Q.begin() + uIndex);

          if (m_dist[get_index(mol, u)] == std::numeric_limits<Size>::max())
            // all remaining atoms are inaccessible from source
            break;

          FOREACH_NBR (v, u, mol, MoleculeType) {
            Size alt = m_dist[get_index(mol, u)] + 1;

            if (alt < m_dist[get_index(mol, *v)]) {
              m_dist[get_index(mol, *v)] = alt;
              m_prev[get_index(mol, *v)] = u;
            }
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
