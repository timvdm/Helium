#ifndef HELIUM_FLOYDWARSHALL_H
#define HELIUM_FLOYDWARSHALL_H

#include <vector>
#include <limits>

#include <Helium/molecule.h>
#include <Helium/distancematrix.h>

namespace Helium {

  /**
   * @brief Floyd-Warshall shortest-path algorithm.
   *
   * @param mol The molecule.
   */
  template<typename MoleculeType>
  DistanceMatrix floyd_warshall(const MoleculeType &mol)
  {
    DistanceMatrix dist(num_atoms(mol), 0, DistanceMatrix::infinity());

    FOREACH_BOND (bond, mol)
      dist(get_index(mol, get_source(mol, *bond)), get_index(mol, get_target(mol, *bond))) = 1;

    for (Size k = 0; k < num_atoms(mol); ++k)
      for (Size i = 0; i < num_atoms(mol); ++i)
        for (Size j = 0; j < num_atoms(mol); ++j) {
          if (dist(i, k) == DistanceMatrix::infinity() || dist(k, j) == DistanceMatrix::infinity())
            continue;
          if (dist(i, k) + dist(k, j) < dist(i, j))
            dist(i, j) = dist(i, k) + dist(k, j);
        }
    
    return dist;
  }

  /*
  template<typename MoleculeType>
  class FloydWarshall
  {
    public:
      template<typename AtomType>
      FloydWarshall(const MoleculeType &mol) : m_mol(mol)
      {
        // distance from source to other atoms
        m_dist.resize(num_atoms(mol), std::numeric_limits<Size>::max());
        // previous atom in shortest path from source
        m_prev.resize(num_atoms(mol), molecule_traits<MoleculeType>::null_atom());

        // distance from source to source
        m_dist[get_index(mol, source)] = 0;
        // add all atoms in mol to Q
        std::vector<AtomType> Q;
        FOREACH_ATOM (atom, mol)
          Q.push_back(*atom);

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

          FOREACH_NBR (v, u, mol) {
            Size alt = m_dist[get_index(mol, u)] + 1;

            if (alt < m_dist[get_index(mol, *v)]) {
              m_dist[get_index(mol, *v)] = alt;
              m_prev[get_index(mol, *v)] = u;
            }
          }
        }
      }

      template<typename AtomType>
      Size distance(AtomType target) const
      {
        return m_dist[get_index(m_mol, target)];
      }

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
      const MoleculeType &m_mol;
      std::vector<Size> m_dist;
      std::vector<typename molecule_traits<MoleculeType>::atom_type> m_prev;
  };
*/

}

#endif
