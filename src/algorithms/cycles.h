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
#include <Helium/smiles.h>
#include <Helium/ring.h>
#include <Helium/hemol.h>

namespace Helium {

  /**
   * @brief Get the cyclomatic number.
   *
   * The cyclomatic number is defined by \f$m - n + c\f$ where \f$n\f$ is the
   * number of atoms, \f$m\f$ the number of bonds and \f$c\f$ the number of
   * connected components. This formula is known as Cauchy's formula. The
   * cyclomatic number is the same as the nullity or first Betti's number.
   *
   * @note Complexity: @f$O(1)@f$
   * @ingroup Production
   * @note Phase: Production
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
   * @note Complexity: @f$O(1)@f$
   * @ingroup Production
   * @note Phase: Production
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
   * @note Complexity: @f$O(n)@f$
   * @ingroup Production
   * @note Phase: Production
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
   * @note Complexity: @f$O(2^n)@f$
   * @ingroup Production
   * @note Phase: Production
   *
   * @param mol The molecule.
   * @param cyclomaticNumber The cyclomatic number.
   * @param cyclicAtoms Atom cycle membership.
   * @param cyclicBonds Bond cycle membership.
   *
   * @return The set of relevant cycles.
   */
  template<typename MoleculeType>
  RingSet<MoleculeType> relevant_cycles(const MoleculeType &mol, Size cyclomaticNumber,
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
      if (rings.size() >= cyclomaticNumber && lastSize < size) {
        // check if all cyclic atoms/bonds are convered
        bool done = true;
        FOREACH_ATOM_T (atom, mol, MoleculeType)
          if (cyclicAtoms[get_index(mol, *atom)] && !rings.isAtomInRing(*atom)) {
            done = false;
            break;
          }
        FOREACH_BOND_T (bond, mol, MoleculeType)
          if (cyclicBonds[get_index(mol, *bond)] && !rings.isBondInRing(*bond)) {
            done = false;
            break;
          }

        if (done)
          break;
      }
      // sanity check
      if (size > num_atoms(mol))
        break;
      // create query
      std::string smiles = "*1" + std::string(size - 1, '*') + "1";

      Smiles SMILES;
      HeMol cycleMol;
      SMILES.read(smiles, cycleMol);

      // find all cycles of size
      MappingList mappings;
      impl::CycleBondMatcher<MoleculeType, HeMol> bondMatcher;
      FOREACH_ATOM_T (atom, mol, MoleculeType) {
        impl::CycleAtomMatcher<MoleculeType, HeMol> atomMatcher(get_index(mol, *atom));

        if (isomorphism_search(mol, *atom, cycleMol, mappings, atomMatcher, bondMatcher)) {
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
  RingSet<MoleculeType> relevant_cycles(const MoleculeType &mol)
  {
    std::vector<bool> cyclicAtoms, cyclicBonds;
    cycle_membership(mol, cyclicAtoms, cyclicBonds);
    return relevant_cycles(mol, cyclomatic_number(mol), cyclicAtoms, cyclicBonds);
  }

  /* SLOWER than regular relevant_cycles()
  template<typename MoleculeType>
  RingSet<MoleculeType> relevant_cycles_substructure(const MoleculeType &mol, Size cyclomaticNumber)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

    std::vector<bool> cycle_atoms;
    std::vector<bool> cycle_bonds;
    cycle_membership(mol, cycle_atoms, cycle_bonds);

    Substructure<MoleculeType> substructure(mol, cycle_atoms, cycle_bonds);

    std::vector<unsigned int> components = connected_bond_components(substructure);
    unsigned int numComponents = unique_elements(components);

    RingSet<MoleculeType> cycles(mol);

    for (unsigned int i = 0; i < numComponents; ++i) {
      std::vector<bool> atoms(num_atoms(mol));
      std::vector<bool> bonds(num_bonds(mol));

      FOREACH_BOND_T (bond, substructure, Substructure<MoleculeType>) {
        if (i == components[get_index(substructure, *bond)]) {
          bonds[substructure.oldBondIndex(*bond)] = true;
          atoms[substructure.oldAtomIndex(get_source(substructure, *bond))] = true;
          atoms[substructure.oldAtomIndex(get_target(substructure, *bond))] = true;
        }
      }

      Substructure<MoleculeType> component(mol, atoms, bonds);

      RingSet<Substructure<MoleculeType> > componentCycles = relevant_cycles(component);

      for (std::size_t j = 0; j < componentCycles.size(); ++j) {
        const Ring<Substructure<MoleculeType> > &ring = componentCycles.ring(j);
        std::vector<atom_type> atoms;
        for (std::size_t k = 0; k < ring.size(); ++k)
          atoms.push_back(get_atom(mol, component.oldAtomIndex(ring.atom(k))));
        cycles.addRing(Ring<MoleculeType>(mol, atoms));
      }
    }

    return cycles;
  }

  template<typename MoleculeType>
  RingSet<MoleculeType> relevant_cycles_substructure(const MoleculeType &mol)
  {
    return relevant_cycles2(mol, cyclomatic_number(mol));
  }
  */

}

#endif
