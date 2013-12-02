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
#ifndef HELIUM_SUBSTRUCTURE_H
#define HELIUM_SUBSTRUCTURE_H

#include <vector>
#include <algorithm>
#include <cassert>

#include <Helium/tie.h>
#include <Helium/hemol.h>

namespace Helium {

  //@cond dev

  namespace impl {

    template<typename SubstructureType>
    class substructure_atom_iterator
    {
      public:
        typedef typename SubstructureType::molecule_type molecule_type;
        typedef typename molecule_traits<molecule_type>::atom_type atom_type;
        typedef typename molecule_traits<molecule_type>::atom_iter atom_iter;

        substructure_atom_iterator()
        {
        }

        substructure_atom_iterator(const SubstructureType *substructure, const atom_iter &iter, const atom_iter &end)
          : m_substructure(substructure), m_iter(iter), m_end(end)
        {
          while (m_iter != m_end && m_substructure->isHidden(*m_iter))
            ++m_iter;
        }

        atom_type operator*() const
        {
          return *m_iter;
        }

        substructure_atom_iterator<SubstructureType>& operator++()
        {
          ++m_iter;
          while (m_iter != m_end && m_substructure->isHidden(*m_iter))
            ++m_iter;
          return *this;
        }

        substructure_atom_iterator<SubstructureType> operator++(int)
        {
          substructure_atom_iterator<SubstructureType> tmp = *this;
          ++m_iter;
          while (m_iter != m_end && m_substructure->isHidden(*m_iter))
            ++m_iter;
          return tmp;
        }

        bool operator!=(const substructure_atom_iterator<SubstructureType> &other)
        {
          return m_iter != other.m_iter;
        }

      private:
        const SubstructureType *m_substructure;
        atom_iter m_iter;
        atom_iter m_end;
    };

    template<typename SubstructureType>
    class substructure_bond_iterator
    {
        public:
          typedef typename SubstructureType::molecule_type molecule_type;
          typedef typename molecule_traits<molecule_type>::bond_type bond_type;
          typedef typename molecule_traits<molecule_type>::bond_iter bond_iter;

          substructure_bond_iterator()
          {
          }

          substructure_bond_iterator(const SubstructureType *substructure, const bond_iter &iter, const bond_iter &end)
            : m_substructure(substructure), m_iter(iter), m_end(end)
          {
            while (m_iter != m_end && m_substructure->isHidden(*m_iter))
              ++m_iter;
          }

          bond_type operator*() const
          {
            return *m_iter;
          }

          substructure_bond_iterator<SubstructureType>& operator++()
          {
            ++m_iter;
            while (m_iter != m_end && m_substructure->isHidden(*m_iter))
              ++m_iter;
            return *this;
          }

          substructure_bond_iterator<SubstructureType> operator++(int)
          {
            substructure_bond_iterator<SubstructureType> tmp = *this;
            ++m_iter;
            while (m_iter != m_end && m_substructure->isHidden(*m_iter))
              ++m_iter;
            return tmp;
          }

          bool operator!=(const substructure_bond_iterator<SubstructureType> &other)
          {
            return m_iter != other.m_iter;
          }

        private:
          const SubstructureType *m_substructure;
          bond_iter m_iter;
          bond_iter m_end;
    };

    template<typename SubstructureType>
    class substructure_incident_iterator
    {
      public:
        typedef typename SubstructureType::molecule_type molecule_type;
        typedef typename molecule_traits<molecule_type>::bond_type bond_type;
        typedef typename molecule_traits<molecule_type>::incident_iter incident_iter;

        substructure_incident_iterator()
        {
        }

        substructure_incident_iterator(const SubstructureType *substructure, const incident_iter &iter, const incident_iter &end)
          : m_substructure(substructure), m_iter(iter), m_end(end)
        {
          while (m_iter != m_end && m_substructure->isHidden(*m_iter))
            ++m_iter;
        }

        bond_type operator*() const
        {
          return *m_iter;
        }

        substructure_incident_iterator<SubstructureType>& operator++()
        {
          ++m_iter;
          while (m_iter != m_end && m_substructure->isHidden(*m_iter))
            ++m_iter;
          return *this;
        }

        substructure_incident_iterator<SubstructureType> operator++(int)
        {
          substructure_incident_iterator<SubstructureType> tmp = *this;
          ++m_iter;
          while (m_iter != m_end && m_substructure->isHidden(*m_iter))
            ++m_iter;
          return tmp;
        }

        bool operator!=(const substructure_incident_iterator<SubstructureType> &other)
        {
          return m_iter != other.m_iter;
        }

      private:
        const SubstructureType *m_substructure;
        incident_iter m_iter;
        incident_iter m_end;
    };

    template<typename SubstructureType>
    class substructure_nbr_iterator
    {
        typedef typename molecule_traits<SubstructureType>::atom_type atom_type;
        typedef typename molecule_traits<SubstructureType>::bond_type bond_type;
        typedef typename molecule_traits<SubstructureType>::incident_iter incident_iter;
      public:
        substructure_nbr_iterator()
        {
        }

        substructure_nbr_iterator(atom_type atom, incident_iter iter) : m_atom(atom), m_iter(iter)
        {
        }

        atom_type operator*() const
        {
          // FIXME
          return (*m_iter).other(m_atom);
        }

        substructure_nbr_iterator& operator++()
        {
          ++m_iter;
          return *this;
        }

        substructure_nbr_iterator operator++(int)
        {
          substructure_nbr_iterator tmp = *this;
          ++m_iter;
          return tmp;
        }

        bool operator!=(const substructure_nbr_iterator &other)
        {
          return m_iter != other.m_iter;
        }

      private:
        atom_type m_atom;
        incident_iter m_iter;
    };

  }

  //@endcond

  template<typename MoleculeType>
  class Substructure
  {
    public:
      // @cond dev

      typedef MoleculeType molecule_type;

      typedef typename molecule_traits<molecule_type>::atom_type atom_type;
      typedef typename molecule_traits<molecule_type>::bond_type bond_type;

      typedef impl::substructure_atom_iterator<Substructure> atom_iter;
      typedef impl::substructure_bond_iterator<Substructure> bond_iter;
      typedef impl::substructure_incident_iterator<Substructure> incident_iter;
      typedef impl::substructure_nbr_iterator<Substructure> nbr_iter;

      typedef impl::substructure_atom_iterator<Substructure> const_atom_iter;
      typedef impl::substructure_bond_iterator<Substructure> const_bond_iter;

      Substructure(const molecule_type &mol, const std::vector<bool> &atoms,
          const std::vector<bool> &bonds) : m_mol(mol),
          m_atoms(atoms), m_bonds(bonds)
      {
        assert(atoms.size() == num_atoms(mol));
        assert(bonds.size() == num_bonds(mol));
        m_numAtoms = std::count(atoms.begin(), atoms.end(), true);
        m_numBonds = std::count(bonds.begin(), bonds.end(), true);

        Index index = 0;
        FOREACH_ATOM (atom, mol, MoleculeType)
          if (m_atoms[get_index(mol, *atom)])
            m_atomIndices.push_back(index++);
          else
            m_atomIndices.push_back(-1);

        index = 0;
        FOREACH_BOND (bond, mol, MoleculeType)
          if (m_bonds[get_index(mol, *bond)])
            m_bondIndices.push_back(index++);
          else
            m_bondIndices.push_back(-1);

        assert(m_atomIndices.size() == num_atoms(mol));
        assert(m_bondIndices.size() == num_bonds(mol));
      }

      Size numAtoms() const
      {
        return m_numAtoms;
      }

      Size numBonds() const
      {
        return m_numBonds;
      }

      std::pair<atom_iter, atom_iter> atoms() const
      {
        typename molecule_traits<molecule_type>::atom_iter begin, end;
        TIE(begin, end) = get_atoms(m_mol);
        return std::make_pair(atom_iter(this, begin, end), atom_iter(this, end, end));
      }

      std::pair<bond_iter, bond_iter> bonds() const
      {
        typename molecule_traits<molecule_type>::bond_iter begin, end;
        TIE(begin, end) = get_bonds(m_mol);
        return std::make_pair(bond_iter(this, begin, end), bond_iter(this, end, end));
      }

      atom_type atom(Index index) const
      {
        ++index;
        for (std::size_t i = 0; i < m_atoms.size(); ++i) {
          if (m_atoms[i])
            --index;
          if (!index)
            return get_atom(m_mol, i);
        }
        assert(0);
        return atom_type();
      }

      bond_type bond(Index index) const
      {
        ++index;
        for (std::size_t i = 0; i < m_bonds.size(); ++i) {
          if (m_bonds[i])
            --index;
          if (!index)
            return get_bond(m_mol, i);
        }
        assert(0);
        return bond_type();
      }

      static Index null_index()
      {
        return -1;
      }

      static atom_type null_atom()
      {
        return molecule_traits<molecule_type>::null_atom();
      }

      static bond_type null_bond()
      {
        return molecule_traits<molecule_type>::null_bond();
      }

      const molecule_type& mol() const
      {
        return m_mol;
      }

      bool isHidden(atom_type atom) const
      {
        return !m_atoms[get_index(m_mol, atom)];
      }

      bool isHidden(bond_type bond) const
      {
        return !m_bonds[get_index(m_mol, bond)];
      }

      Index newAtomIndex(atom_type atom) const
      {
        return m_atomIndices[get_index(m_mol, atom)];
      }

      Index oldAtomIndex(atom_type atom) const
      {
        return index_of(m_atomIndices, newAtomIndex(atom));
      }

      Index newBondIndex(bond_type bond) const
      {
        return m_bondIndices[get_index(m_mol, bond)];
      }

      Index oldBondIndex(bond_type bond) const
      {
        return index_of(m_bondIndices, newBondIndex(bond));
      }

    private:
      const molecule_type &m_mol;
      Size m_numAtoms;
      Size m_numBonds;
      std::vector<bool> m_atoms;
      std::vector<bool> m_bonds;
      std::vector<Index> m_atomIndices;
      std::vector<Index> m_bondIndices;

      //@endcond
  };

  //@cond dev

  //////////////////////////////////////////////////////////////////////////////
  //
  // Molecule
  //
  //////////////////////////////////////////////////////////////////////////////

  /////////////////////////////
  //
  // Atoms
  //
  /////////////////////////////

  template<typename MoleculeType>
  Size num_atoms(const Substructure<MoleculeType> &mol)
  {
    return mol.numAtoms();
  }

  template<typename MoleculeType>
  std::pair<typename molecule_traits<Substructure<MoleculeType> >::atom_iter,
            typename molecule_traits<Substructure<MoleculeType> >::atom_iter>
  get_atoms(const Substructure<MoleculeType> &mol)
  {
    return mol.atoms();
  }

  template<typename MoleculeType>
  typename molecule_traits<Substructure<MoleculeType> >::atom_type
  get_atom(const Substructure<MoleculeType> &mol, Index index)
  {
    return mol.atom(index);
  }


  /////////////////////////////
  //
  // Bonds
  //
  /////////////////////////////

  template<typename MoleculeType>
  Size num_bonds(const Substructure<MoleculeType> &mol)
  {
    return mol.numBonds();
  }

  template<typename MoleculeType>
  std::pair<typename molecule_traits<Substructure<MoleculeType> >::bond_iter,
            typename molecule_traits<Substructure<MoleculeType> >::bond_iter>
  get_bonds(const Substructure<MoleculeType> &mol)
  {
    return mol.bonds();
  }

  template<typename MoleculeType>
  typename molecule_traits<Substructure<MoleculeType> >::bond_type
  get_bond(const Substructure<MoleculeType> &mol, Index index)
  {
    return mol.bond(index);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Atom
  //
  //////////////////////////////////////////////////////////////////////////////

  template<typename MoleculeType>
  Index get_index(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type atom)
  {
    return mol.newAtomIndex(atom);
  }

  template<typename MoleculeType>
  std::pair<typename molecule_traits<Substructure<MoleculeType> >::incident_iter,
            typename molecule_traits<Substructure<MoleculeType> >::incident_iter>
  get_bonds(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type atom)
  {
    typename molecule_traits<HeMol>::incident_iter begin, end; // <- use of HeMol is intended!!
    TIE(begin, end) = get_bonds(mol.mol(), atom);
    typedef typename molecule_traits<Substructure<MoleculeType> >::incident_iter incident_iter;
    return std::make_pair(incident_iter(&mol, begin, end), incident_iter(&mol, end, end));
  }

  template<typename MoleculeType>
  std::pair<typename molecule_traits<Substructure<MoleculeType> >::nbr_iter,
            typename molecule_traits<Substructure<MoleculeType> >::nbr_iter>
  get_nbrs(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type atom)
  {
    typename molecule_traits<Substructure<MoleculeType> >::incident_iter begin, end;
    TIE(begin, end) = get_bonds(mol, atom);
    typedef typename molecule_traits<Substructure<MoleculeType> >::nbr_iter nbr_iter;
    return std::make_pair(nbr_iter(atom, begin), nbr_iter(atom, end));
  }

  template<typename MoleculeType>
  bool is_aromatic(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type atom)
  {
    return is_aromatic(mol.mol(), atom);
  }

  template<typename MoleculeType>
  int get_element(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type atom)
  {
    return get_element(mol.mol(), atom);
  }

  template<typename MoleculeType>
  int get_mass(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type atom)
  {
    return get_mass(mol.mol(), atom);
  }

  template<typename MoleculeType>
  int get_degree(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type atom)
  {
    typename molecule_traits<Substructure<MoleculeType> >::incident_iter bond, end_bonds;
    TIE(bond, end_bonds) = get_bonds(const_cast<Substructure<MoleculeType>&>(mol), atom);
    int d = 0;
    for (; bond != end_bonds; ++bond)
      ++d;
    return d;
  }

  template<typename MoleculeType>
  int num_hydrogens(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type atom)
  {
    return num_hydrogens(mol.mol(), atom);
  }

  template<typename MoleculeType>
  int get_charge(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type atom)
  {
    return get_charge(mol.mol(), atom);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Bond
  //
  //////////////////////////////////////////////////////////////////////////////

  template<typename MoleculeType>
  Index get_index(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::bond_type bond)
  {
    return mol.newBondIndex(bond);
  }

  template<typename MoleculeType>
  typename molecule_traits<Substructure<MoleculeType> >::atom_type
  get_source(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::bond_type bond)
  {
    return get_source(mol.mol(), bond);
  }

  template<typename MoleculeType>
  typename molecule_traits<Substructure<MoleculeType> >::atom_type
  get_target(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::bond_type bond)
  {
    return get_target(mol.mol(), bond);
  }

  template<typename MoleculeType>
  typename molecule_traits<Substructure<MoleculeType> >::atom_type
  get_other(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::bond_type bond,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type atom)
  {
    return get_other(mol.mol(), bond, atom);
  }

  template<typename MoleculeType>
  bool is_aromatic(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::bond_type bond)
  {
    return is_aromatic(mol.mol(), bond);
  }

  template<typename MoleculeType>
  int get_order(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::bond_type bond)
  {
    return get_order(mol.mol(), bond);
  }

  template<typename MoleculeType>
  typename molecule_traits<Substructure<MoleculeType> >::bond_type
  get_bond(const Substructure<MoleculeType> &mol,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type source,
      typename molecule_traits<Substructure<MoleculeType> >::atom_type target)
  {
    typename molecule_traits<Substructure<MoleculeType> >::incident_iter bond, end_bonds;
    TIE(bond, end_bonds) = get_bonds(mol, source);
    for (; bond != end_bonds; ++bond)
      if (get_other(mol.mol(), *bond, source) == target)
        return *bond;
    return 0;
  }

  //@endcond

}

#endif
