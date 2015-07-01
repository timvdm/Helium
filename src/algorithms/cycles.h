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
#ifndef HELIUM_CYCLES_H
#define HELIUM_CYCLES_H

#include <Helium/algorithms/components.h>
#include <Helium/algorithms/dfs.h>
#include <Helium/algorithms/isomorphism.h>
#include <Helium/algorithms/dijkstra.h>
#include <Helium/smiles.h>
#include <Helium/ring.h>
#include <Helium/hemol.h>

namespace Helium {

  /**
   * @file cycles.h
   * @brief Cycle detection and perception.
   */

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
  Size cyclomatic_number(const MoleculeType &mol, Size numComponents)
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
  Size cyclomatic_number(const MoleculeType &mol)
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

        int size() const
        {
          return m_edges.size();
        }

        const std::set<Index>& edges() const
        {
          return m_edges;
        }

      private:
        std::set<Index> m_edges;
    };

    class CycleBitMatrix
    {
      public:
        CycleBitMatrix(int cols) : m_rows(0), m_cols(cols)
        {
        }

        CycleBitMatrix(int rows, int cols) : m_rows(rows), m_cols(cols)
        {
          m_data.resize(rows * cols);
        }

        void addRow()
        {
          ++m_rows;
          m_data.resize(m_data.size() + m_cols);
        }

        void popRow()
        {
          --m_rows;
          m_data.resize(m_data.size() - m_cols);
        }

        template<typename MoleculeType>
        void addRow(const MoleculeType &mol, const std::vector<Index> &atoms)
        {
          int offset = m_data.size();
          addRow();
          m_data[offset + get_index(mol, get_bond(mol, get_atom(mol, atoms[0]), get_atom(mol, atoms[atoms.size() - 1])))] = true;
          for (std::size_t i = 1; i < atoms.size(); ++i)
            m_data[offset + get_index(mol, get_bond(mol, get_atom(mol, atoms[i - 1]), get_atom(mol, atoms[i])))] = true;
        }

        int rows() const
        {
          return m_rows;
        }

        int cols() const
        {
          return m_cols;
        }

        bool get(int row, int col) const
        {
          return m_data[index(row, col)];
        }

        void set(int row, int col, bool value)
        {
          m_data[index(row, col)] = value;
        }

        int eliminate(int x = 0, int y = 0)
        {
          while (x < m_cols && y < m_rows) {
            int i = indexOf(x, y);

            // this column is done, continue with next column
            if (i < 0)
              return eliminate(x + 1, y);

            // swap rows if needed
            if (i != y)
              swapRows(i, y);

            for (int j = y + 1; j < m_rows; ++j)
              if (m_data[index(j, x)])
                xorRows(y, j);

            ++y;
          }
          return y;
        }

      private:
        std::size_t index(int row, int col) const
        {
          return m_cols * row + col;
        }

        int indexOf(int x, int y) const
        {
          for (int i = y; i < m_rows; ++i)
            if (m_data[index(i, x)])
              return i;
          return -1;
        }

        void swapRows(int i, int j)
        {
          std::vector<bool> tmp(m_cols);
          // i -> tmp
          std::copy(m_data.begin() + i * m_cols, m_data.begin() + (i + 1) * m_cols, tmp.begin());
          // j -> i
          std::copy(m_data.begin() + j * m_cols, m_data.begin() + (j + 1) * m_cols, m_data.begin() + i * m_cols);
          // tmp -> j
          std::copy(tmp.begin(), tmp.end(), m_data.begin() + j * m_cols);
        }

        void xorRows(int i, int j)
        {
          for (int col = 0; col < m_cols; ++col)
            m_data[index(j, col)] = m_data[index(i, col)] ^ m_data[index(j, col)];
        }

        std::vector<bool> m_data;
        int m_rows;
        int m_cols;
    };

    inline std::ostream& operator<<(std::ostream &os, const CycleBitMatrix &m)
    {
      for (int i = 0; i < m.rows(); ++i) {
        os << "[ ";
        for (int j = 0; j < m.cols(); ++j)
          os << m.get(i, j) << " ";
        os << "]" << std::endl;
      }
      return os;
    }

    inline bool is_relevant(const CycleBitMatrix &matrix, const IsomorphismCycle &cycle)
    {
      CycleBitMatrix m(matrix);

      std::set<Index>::const_iterator e;

      int row = m.rows();
      m.addRow();
      for (e = cycle.edges().begin(); e != cycle.edges().end(); ++e)
        m.set(row, *e, true);

      int rank = m.rows();
      bool relevant = (rank == m.eliminate());
      return relevant;
    }

    template<typename MoleculeType>
    bool is_ringset_complete(const MoleculeType &mol, const RingSet<MoleculeType> &rings,
        const std::vector<bool> &cyclicAtoms, const std::vector<bool> &cyclicBonds)
    {
      // check if all cyclic atoms/bonds are covered
      for (auto &atom : get_atoms(mol))
        if (cyclicAtoms[get_index(mol, atom)] && !rings.isAtomInRing(atom))
          return false;
      for (auto &bond : get_bonds(mol))
        if (cyclicBonds[get_index(mol, bond)] && !rings.isBondInRing(bond))
          return false;
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
   * @param mol The molecule.
   * @param cyclomaticNumber The cyclomatic number.
   * @param cyclicAtoms Atom cycle membership.
   * @param cyclicBonds Bond cycle membership.
   *
   * @return The set of relevant cycles.
   */
  template<typename MoleculeType>
  RingSet<MoleculeType> relevant_cycles_isomorphism(const MoleculeType &mol, Size cyclomaticNumber,
      const std::vector<bool> &cyclicAtoms, const std::vector<bool> &cyclicBonds)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

    RingSet<MoleculeType> rings(mol);

    // current cycle size
    unsigned int size = 3;
    unsigned int lastSize = 3;

    impl::CycleBitMatrix matrix(num_bonds(mol));

    while (true) {
      //std::cout << "size: " << size << ", lastSize: " << lastSize << ", nullity: " << cyclomaticNumber << ", count: " << relevant.size() << std::endl;
      if (rings.size() >= cyclomaticNumber && lastSize < size &&
          impl::is_ringset_complete(mol, rings, cyclicAtoms, cyclicBonds))
        break;
      // sanity check
      if (size > num_atoms(mol))
        break;
      // create query
      std::string smiles = "*1" + std::string(size - 1, '*') + "1";

      //Smiles SMILES;
      HeMol cycleMol;
      //SMILES.read(smiles, cycleMol);

      for (unsigned int i = 0; i < size; ++i)
        add_atom(cycleMol);
      for (unsigned int i = 1; i < size; ++i)
        add_bond(cycleMol, get_atom(cycleMol, i - 1), get_atom(cycleMol, i));
      add_bond(cycleMol, get_atom(cycleMol, 0), get_atom(cycleMol, size - 1));

      // find all cycles of size
      MappingList mappings;
      impl::CycleBondMatcher<MoleculeType, HeMol> bondMatcher;
      for (auto &atom : get_atoms(mol)) {
        impl::CycleAtomMatcher<MoleculeType, HeMol> atomMatcher(get_index(mol, atom));

        if (isomorphism_search(mol, atom, cycleMol, mappings, atomMatcher, bondMatcher)) {
          for (std::size_t i = 0; i < mappings.maps.size(); ++i) {

            impl::IsomorphismCycle cycle(mol, mappings.maps[i]);

            if (is_relevant(matrix, cycle)) {
              std::vector<atom_type> atoms;
              for (std::size_t j = 0; j < mappings.maps[i].size(); ++j)
                atoms.push_back(get_atom(mol, mappings.maps[i][j]));
              rings.addRing(Ring<MoleculeType>(mol, atoms));
              lastSize = size;
            }
          }
        }
      }

      // add cycles of size to matrix
      for (std::size_t i = 0; i < rings.size(); ++i) {
        if (rings.ring(i).size() < size)
          continue;
        int row = matrix.rows();
        matrix.addRow();
        for (std::size_t j = 0; j < rings.ring(i).size(); ++j)
          matrix.set(row, get_index(mol, rings.ring(i).bond(j)), true);
      }

      int rank = matrix.eliminate();
      while (matrix.rows() > rank)
        matrix.popRow();

      // increment cycle size
      ++size;
    }

    return rings;
  }


  /**
   * @overload
   */
  template<typename MoleculeType>
  RingSet<MoleculeType> relevant_cycles_isomorphism(const MoleculeType &mol)
  {
    std::vector<bool> cyclicAtoms, cyclicBonds;
    cycle_membership(mol, cyclicAtoms, cyclicBonds);
    return relevant_cycles_isomorphism(mol, cyclomatic_number(mol), cyclicAtoms, cyclicBonds);
  }


  namespace impl {

    template<typename MoleculeType, typename AtomType>
    bool paths_intersection_is_source(const Dijkstra<MoleculeType> &dijkstra,
        const MoleculeType &mol, const AtomType &y, const AtomType &z)
    {
      AtomType u = y;
      while (dijkstra.prev()[get_index(mol, u)] != molecule_traits<MoleculeType>::null_atom()) {
        AtomType v = z;
        while (dijkstra.prev()[get_index(mol, v)] != molecule_traits<MoleculeType>::null_atom()) {
          if (u == v)
            return false;
          v = dijkstra.prev()[get_index(mol, v)];
        }
        u = dijkstra.prev()[get_index(mol, u)];
      }

      return true;
      /*
      std::vector<AtomType> Py = dijkstra.path(y);
      std::vector<AtomType> Pz = dijkstra.path(z);
      std::sort(Py.begin(), Py.end());
      std::sort(Pz.begin(), Pz.end());
      std::vector<AtomType> intersection;
      std::set_intersection(Py.begin(), Py.end(), Pz.begin(), Pz.end(),
          std::back_inserter(intersection));
      return intersection.size() == 1 && intersection[0] == dijkstra.source();
      */
    }

    class CycleFamily
    {
      public:
        CycleFamily(Index r, Index p, Index q, const std::vector<Index> &prototype)
          : m_r(r), m_p(p), m_q(q), m_x(std::numeric_limits<Index>::max()),
            m_prototype(prototype)
        {
        }

        CycleFamily(Index r, Index p, Index q, Index x, const std::vector<Index> &prototype)
          : m_r(r), m_p(p), m_q(q), m_x(x), m_prototype(prototype)
        {
        }

        Index r() const
        {
          return m_r;
        }

        Index p() const
        {
          return m_p;
        }

        Index q() const
        {
          return m_q;
        }

        Index x() const
        {
          return m_x;
        }

        const std::vector<Index>& prototype() const
        {
          return m_prototype;
        }

        bool isOdd() const
        {
          return m_x == std::numeric_limits<Index>::max();
        }

        bool isEven() const
        {
          return m_x != std::numeric_limits<Index>::max();
        }

        bool operator<(const CycleFamily &other) const
        {
          return m_prototype.size() < other.m_prototype.size();
        }

      private:
        Index m_r;
        Index m_p;
        Index m_q;
        Index m_x;
        std::vector<Index> m_prototype;
    };

    inline void list_paths(const CycleBitMatrix &Dr, Index r, Index x, std::vector<Index> &current,
        std::vector<std::vector<Index> > &result)
    {
      // add x at the head-end of current path
      current.insert(current.begin(), x);

      //std::cout << current << std::endl;

      if (x == r) {
        result.push_back(current);
      } else {
        // for any z such that (x, z) is in Ur
        for (std::size_t z = 0; z <= r; ++z)
          if (Dr.get(x, z)) {
            std::vector<Index> currentCopy(current);
            list_paths(Dr, r, z, currentCopy, result);
          }
      }

      //current.pop_back();
    }

    inline std::vector<std::vector<Index> > list_paths(const CycleBitMatrix &Dr, Index r, Index x)
    {
      std::vector<std::vector<Index> > result;
      std::vector<Index> current;
      list_paths(Dr, r, x, current, result);
      return result;
    }

    template<typename MoleculeType>
    void add_atom_cycle_to_matrix(const MoleculeType &mol, CycleBitMatrix &B,
        const std::vector<Index> &cycle)
    {
      int row = B.rows();
      B.addRow();
      for (std::size_t i = 1; i < cycle.size(); ++i)
        B.set(row, get_index(mol, get_bond(mol, get_atom(mol, cycle[i - 1]), get_atom(mol, cycle[i]))), true);
      B.set(row, get_index(mol, get_bond(mol, get_atom(mol, cycle.front()),
              get_atom(mol, cycle.back()))), true);
    }

    template<typename MoleculeType>
    bool is_relevant(const MoleculeType &mol, const CycleBitMatrix &matrix,
        const std::vector<Index> &cycle)
    {
      CycleBitMatrix m(matrix);

      add_atom_cycle_to_matrix(mol, m, cycle);

      int rank = m.rows();
      bool relevant = (rank == m.eliminate());
      return relevant;
    }

  }


  /**
   * @brief Find the relevant cycles.
   *
   * The set of relevant cycles is formed by taking the union of all the
   * minimum cycles bases. An alternative definition is that a cycle is
   * relevant if it is not the sum of smaller cycles.
   *
   * @param mol The molecule.
   * @param cyclomaticNumber The cyclomatic number.
   * @param cyclicAtoms Atom cycle membership.
   * @param cyclicBonds Bond cycle membership.
   *
   * @return The set of relevant cycles.
   */
  template<typename MoleculeType>
  RingSet<MoleculeType> relevant_cycles_vismara(const MoleculeType &mol, Size cyclomaticNumber,
      const std::vector<bool> &cyclicAtoms, const std::vector<bool> &cyclicBonds)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

    // V = {1, ..., n}  |V| = n
    //
    // E = { (x, y) | there is an edge between x and y }  |E| = m
    //
    // w((x, y)): weight of edge (x, y) (always 1)
    //
    // G = (V, E)
    //
    // pi: a random total order on V (i.e. atom indexes)
    //
    // Vr = { z in V | there is a shortest path in G from r to z that passes only
    //                 through vertices which precede r in the ordering pi }

    // digraph Dr = (Vr, Ur)
    //
    // Ur = { directed edge (y, z) | (z, y) belongs to a shortest path in G from
    //                               r to y that passes only to through vertices
    //                               which precede r in the ordering pi }
    std::vector<impl::CycleBitMatrix> D;

    std::vector<impl::CycleFamily> families;

    // for all r in V do
    for (auto &r : get_atoms(mol)) {
      // compute Vr and for all t in Vr find a shortest path P(r, t) from r to t
      /*
      std::vector<bool> atomMask(get_index(mol, r) + 1, true);
      atomMask.resize(num_atoms(mol));
      */
      std::vector<bool> atomMask(num_atoms(mol));
      for (std::size_t i = 0; i < get_index(mol, r) + 1; ++i)
        atomMask[i] = cyclicAtoms[i];
      Dijkstra<MoleculeType> dijkstra(mol, r, atomMask);

      //std::cout << atomMask << std::endl;

      std::vector<Index> Vr;
      for (std::size_t i = 0; i <= get_index(mol, r); ++i)
        if (dijkstra.distance(get_atom(mol, i)) < dijkstra.infinity())
          Vr.push_back(i);

      //std::cout << "Vr: " << Vr << std::endl;

      D.push_back(impl::CycleBitMatrix(num_atoms(mol), num_atoms(mol)));
      impl::CycleBitMatrix &Dr = D.back();

      // for all y in Vr do
      for (std::size_t i = 0; i < Vr.size(); ++i) {
        atom_type y = get_atom(mol, Vr[i]);

        // S <- {}
        std::vector<atom_type> S;

        // for all z in Vr such that z is adjacent to y
        for (auto &z : get_nbrs(mol, y)) {
          if (!std::binary_search(Vr.begin(), Vr.end(), get_index(mol, z)))
            continue;

          // if d(r, z) + w((z, y)) = d(r, y) then
          if (dijkstra.distance(z) + 1 == dijkstra.distance(y)) {
            // S <- S U {z}
            S.push_back(z);
            //std::cout << "edge: (" << get_index(mol, y) << ", " << get_index(mol, z) << ")" << std::endl;
            Dr.set(get_index(mol, y), get_index(mol, z), true);
          // else if d(r, z) != d(r, y) + w((z, y))
          //     and pi(z) < pi(y)
          //     and P(r, y) ^ P(r, z) = {r}
          // then
          } else if (dijkstra.distance(z) != dijkstra.distance(y)  + 1 &&
                     get_index(mol, z) < get_index(mol, y) &&
                     impl::paths_intersection_is_source(dijkstra, mol, y, z)) {
            // add to CI' the odd cycle C = P(r, y) + P(r, z) + (z, y)
            //std::cout << "odd cycle family: " << get_index(mol, r) << ", "
            //          << get_index(mol, y) << ", " << get_index(mol, z) << std::endl;
            //std::cout << "    |P(r, y)| = " << dijkstra.distance(y) << std::endl;
            //std::cout << "    |P(r, z)| = " << dijkstra.distance(z) << std::endl;

            std::vector<Index> prototype;
            std::vector<atom_type> Py = dijkstra.path(y);
            std::vector<atom_type> Pz = dijkstra.path(z);
            for (std::size_t j = 0; j < Py.size(); ++j)
              prototype.push_back(get_index(mol, Py[j]));
            for (std::size_t j = Pz.size() - 1; j > 0; --j)
              prototype.push_back(get_index(mol, Pz[j]));

            assert(prototype.size() == Py.size() + Pz.size() - 1);

            families.push_back(impl::CycleFamily(get_index(mol, r), get_index(mol, y),
                get_index(mol, z), prototype));
          }
        }

        // for any pair of vertices p, q in S such that P(r, p) ^ P(r, q) = {r} do
        for (std::size_t j = 0; j < S.size(); ++j)
          for (std::size_t k = j + 1; k < S.size(); ++k) {
            const atom_type &p = S[j];
            const atom_type &q = S[k];
            if (!impl::paths_intersection_is_source(dijkstra, mol, p, q))
              continue;

            // add to CI' the even cycle C = P(r, p) + P(r, q) + (p, y, q)
            //std::cout << "even cycle family: " << get_index(mol, r) << ", "
            //          << get_index(mol, p) << ", " << get_index(mol, y) << ", "
            //          << get_index(mol, q) << std::endl;

            std::vector<Index> prototype;
            std::vector<atom_type> Pp = dijkstra.path(p);
            std::vector<atom_type> Pq = dijkstra.path(q);
            for (std::size_t l = 0; l < Pp.size(); ++l)
              prototype.push_back(get_index(mol, Pp[l]));
            prototype.push_back(get_index(mol, y));
            for (std::size_t l = Pq.size() - 1; l > 0; --l)
              prototype.push_back(get_index(mol, Pq[l]));

            assert(prototype.size() == Pp.size() + Pq.size());

            families.push_back(impl::CycleFamily(get_index(mol, r), get_index(mol, p),
                get_index(mol, q), get_index(mol, y), prototype));
          }
      }
    }

    if (families.empty())
      return RingSet<MoleculeType>(mol);

    //
    // select relevant families
    //

    std::sort(families.begin(), families.end());

    unsigned int lastSize = families.front().prototype().size();

    std::vector<std::size_t> remove;
    impl::CycleBitMatrix B_less(num_bonds(mol));
    for (std::size_t i = 0; i < families.size(); ++i) {
      const impl::CycleFamily &family = families[i];

      // update B_less
      if (lastSize < family.prototype().size()) {
        for (std::size_t j = 0; j < families.size(); ++j) {
          if (families[j].prototype().size() < lastSize)
            continue;
          if (families[j].prototype().size() >= family.prototype().size())
            break;
          add_atom_cycle_to_matrix(mol, B_less, families[j].prototype());
        }

        int rank = B_less.eliminate();
        while (B_less.rows() > rank)
          B_less.popRow();

        lastSize = family.prototype().size();
      }

      if (!is_relevant(mol, B_less, family.prototype()))
        remove.push_back(i);
    }

    //std::cout << "removing " << remove.size() << " families..." << std::endl;

    for (std::size_t i = 0; i < remove.size(); ++i)
      families.erase(families.begin() + remove[remove.size() - i - 1]);

    //
    // enumerate relevant cycles
    //

    lastSize = 0;
    RingSet<MoleculeType> rings(mol);
    for (std::size_t i = 0; i < families.size(); ++i) {
      const impl::CycleFamily &family = families[i];

      if (lastSize < family.prototype().size()) {
        if (rings.size() >= cyclomaticNumber &&
            impl::is_ringset_complete(mol, rings, cyclicAtoms, cyclicBonds))
          break;
        lastSize = family.prototype().size();
      }

      const impl::CycleBitMatrix &Dr = D[family.r()];

      //std::cout << "Dr:" << std::endl << Dr << std::endl;

      std::vector<std::vector<Index> > pPaths = impl::list_paths(Dr, family.r(), family.p());
      std::vector<std::vector<Index> > qPaths = impl::list_paths(Dr, family.r(), family.q());


      //std::cout << "    (r, ..., p): " << pPaths << std::endl;
      //std::cout << "    (r, ..., q): " << qPaths << std::endl;

      if (family.isOdd()) {
        //std::cout << "odd cycle: " << family.r() << " " << family.p() << " "
        //          << family.q() << "  " << family.prototype() << std::endl;
        for (std::size_t j = 0; j < pPaths.size(); ++j)
          for (std::size_t k = 0; k < qPaths.size(); ++k) {
            std::vector<atom_type> atoms;

            for (std::size_t l = 0; l < pPaths[j].size(); ++l)
              atoms.push_back(get_atom(mol, pPaths[j][l]));
            for (std::size_t l = qPaths[k].size() - 1; l > 0; --l)
              atoms.push_back(get_atom(mol, qPaths[k][l]));
            //std::cout << atoms << std::endl;

            rings.addRing(Ring<MoleculeType>(mol, atoms));
          }
      } else {
        //std::cout << "even cycle: " << family.r() << " " << family.p() << " "
        //          << family.q() << " " << family.x() << "  "
        //          << family.prototype() << std::endl;
        for (std::size_t j = 0; j < pPaths.size(); ++j)
          for (std::size_t k = 0; k < qPaths.size(); ++k) {
            std::vector<atom_type> atoms;

            for (std::size_t l = 0; l < pPaths[j].size(); ++l)
              atoms.push_back(get_atom(mol, pPaths[j][l]));
            atoms.push_back(get_atom(mol, family.x()));
            for (std::size_t l = qPaths[k].size() - 1; l > 0; --l)
              atoms.push_back(get_atom(mol, qPaths[k][l]));
            //std::cout << atoms << std::endl;

            rings.addRing(Ring<MoleculeType>(mol, atoms));
          }
      }

    }

    return rings;
  }

  /**
   * @overload
   */
  template<typename MoleculeType>
  RingSet<MoleculeType> relevant_cycles_vismara(const MoleculeType &mol)
  {
    std::vector<bool> cyclicAtoms, cyclicBonds;
    cycle_membership(mol, cyclicAtoms, cyclicBonds);
    return relevant_cycles_vismara(mol, cyclomatic_number(mol),
        cyclicAtoms, cyclicBonds);
  }

  /**
   * @brief Find the relevant cycles.
   *
   * The set of relevant cycles is formed by taking the union of all the
   * minimum cycles bases. An alternative definition is that a cycle is
   * relevant if it is not the sum of smaller cycles.
   *
   * @param mol The molecule.
   *
   * @return The set of relevant cycles.
   */
  template<typename MoleculeType>
  RingSet<MoleculeType> relevant_cycles(const MoleculeType &mol)
  {
    return relevant_cycles_vismara(mol);
  }

}

#endif
