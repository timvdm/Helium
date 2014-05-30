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
#ifndef HELIUM_HEMOL_H
#define HELIUM_HEMOL_H

#include <Helium/config.h>
#include <Helium/molecule.h>
#include <Helium/tie.h>
#include <Helium/util/vector.h>

#include <vector>
#include <istream>
#include <algorithm>
#include <cassert>

namespace Helium {

  //@cond dev

  namespace impl {

    template<typename MoleculeType>
    class nbr_iterator : public std::iterator<std::forward_iterator_tag,
        typename molecule_traits<MoleculeType>::atom_type>
    {
      public:
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

        nbr_iterator()
        {
        }

        nbr_iterator(const atom_type &atom, typename std::vector<bond_type>::iterator iter)
          : m_atom(atom), m_iter(iter)
        {
        }

        atom_type& operator*() const
        {
          m_nbr = (*m_iter).other(m_atom);
          return m_nbr;
        }

        nbr_iterator<MoleculeType>& operator++()
        {
          ++m_iter;
          return *this;
        }

        nbr_iterator<MoleculeType> operator++(int)
        {
          nbr_iterator<MoleculeType> tmp = *this;
          ++m_iter;
          return tmp;
        }

        bool operator==(const nbr_iterator<MoleculeType> &other) const
        {
          return m_iter == other.m_iter;
        }

        bool operator!=(const nbr_iterator<MoleculeType> &other) const
        {
          return m_iter != other.m_iter;
        }

      private:
        atom_type m_atom;
        typename std::vector<bond_type>::iterator m_iter;
        mutable atom_type m_nbr;
    };


    /**
     * @brief Class representing an atom in an HeMol.
     */
    template<typename MoleculeType>
    class AtomImpl
    {
      public:
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

        AtomImpl(MoleculeType *mol = 0, Index index = -1) : m_mol(mol), m_index(index)
        {
        }

        iterator_pair<typename std::vector<bond_type>::iterator> bonds()
        {
          return make_iterator_pair(m_mol->m_adjList[m_index].begin(), m_mol->m_adjList[m_index].end());
        }

        iterator_pair<impl::nbr_iterator<MoleculeType> > nbrs()
        {
          return make_iterator_pair(impl::nbr_iterator<MoleculeType>(*this, m_mol->m_adjList[m_index].begin()),
                                    impl::nbr_iterator<MoleculeType>(*this, m_mol->m_adjList[m_index].end()));
        }

        Index index() const
        {
          return m_index;
        }

        bool isAromatic() const
        {
          return m_mol->m_atomAromatic[m_index];
        }

        void setAromatic(bool value)
        {
          m_mol->m_atomAromatic[m_index] = value;
        }

        int element() const
        {
          return m_mol->m_element[m_index];
        }

        void setElement(int value)
        {
          m_mol->m_element[m_index] = value;
        }

        int mass() const
        {
          return m_mol->m_mass[m_index];
        }

        void setMass(int value)
        {
          m_mol->m_mass[m_index] = value;
        }

        int degree() const
        {
          return m_mol->m_adjList[m_index].size();
        }

        int hydrogens() const
        {
          return m_mol->m_hydrogens[m_index];
        }

        void setHydrogens(int value)
        {
          m_mol->m_hydrogens[m_index] = value;
        }

        int charge() const
        {
          return m_mol->m_charge[m_index];
        }

        void setCharge(int value)
        {
          m_mol->m_charge[m_index] = value;
        }

        bool operator==(const atom_type &other) const
        {
          return m_index == other.m_index;
        }

        bool operator!=(const atom_type &other) const
        {
          return m_index != other.m_index;
        }

        bool operator<(const atom_type &other) const
        {
          return m_index < other.m_index;
        }

        bool operator>(const atom_type &other) const
        {
          return m_index > other.m_index;
        }

      private:
        MoleculeType *m_mol;
        Index m_index;
    };

    template<typename MoleculeType>
    class BondImpl
    {
      public:
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

        BondImpl(MoleculeType *mol = 0, Index index = -1) : m_mol(mol), m_index(index)
        {
        }

        Index index() const
        {
          return m_index;
        }

        atom_type source() const
        {
          return atom_type(m_mol, m_mol->m_source[m_index]);
        }

        atom_type target() const
        {
          return atom_type(m_mol, m_mol->m_target[m_index]);
        }

        atom_type other(const atom_type &atom) const
        {
          Index source = m_mol->m_source[m_index];
          Index target = m_mol->m_target[m_index];
          return source == atom.index() ? atom_type(m_mol, target) : atom_type(m_mol, source);
        }

        bool isAromatic() const
        {
          return m_mol->m_bondAromatic[m_index];
        }

        void setAromatic(bool value)
        {
          m_mol->m_bondAromatic[m_index] = value;
        }

        int order() const
        {
          return m_mol->m_order[m_index];
        }

        void setOrder(int value)
        {
          m_mol->m_order[m_index] = value;
        }

        bool operator==(const bond_type &other) const
        {
          return m_index == other.m_index;
        }

        bool operator!=(const bond_type &other) const
        {
          return m_index != other.m_index;
        }

        bool operator<(const bond_type &other) const
        {
          return m_index < other.m_index;
        }

        bool operator>(const bond_type &other) const
        {
          return m_index > other.m_index;
        }

      private:
        MoleculeType *m_mol;
        Index m_index;
    };

    template<typename MoleculeType>
    class atom_iterator : public std::iterator<std::forward_iterator_tag,
        typename molecule_traits<MoleculeType>::atom_type>
    {
      public:
        atom_iterator()
        {
        }

        atom_iterator(const MoleculeType *mol, Index index = 0) : m_mol(mol), m_index(index)
        {
        }

        typename molecule_traits<MoleculeType>::atom_type& operator*() const
        {
          typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
          m_atom = atom_type(const_cast<MoleculeType*>(m_mol), m_index);
          return m_atom;
        }

        atom_iterator<MoleculeType>& operator++()
        {
          ++m_index;
          return *this;
        }

        atom_iterator<MoleculeType> operator++(int)
        {
          atom_iterator<MoleculeType> tmp = *this;
          ++m_index;
          return tmp;
        }

        bool operator==(const atom_iterator<MoleculeType> &other) const
        {
          return m_index == other.m_index;
        }

        bool operator!=(const atom_iterator<MoleculeType> &other) const
        {
          return m_index != other.m_index;
        }

      private:
        const MoleculeType *m_mol;
        Index m_index;
        mutable typename molecule_traits<MoleculeType>::atom_type m_atom;
    };

    template<typename MoleculeType>
    class bond_iterator : public std::iterator<std::forward_iterator_tag,
        typename molecule_traits<MoleculeType>::bond_type>
    {
      public:
        bond_iterator()
        {
        }

        bond_iterator(const MoleculeType *mol, Index index = 0) : m_mol(mol), m_index(index)
        {
        }

        typename molecule_traits<MoleculeType>::bond_type& operator*() const
        {
          typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
          m_bond = bond_type(const_cast<MoleculeType*>(m_mol), m_index);
          return m_bond;
        }

        bond_iterator<MoleculeType>& operator++()
        {
          ++m_index;
          return *this;
        }

        bond_iterator<MoleculeType> operator++(int)
        {
          bond_iterator<MoleculeType> tmp = *this;
          ++m_index;
          return tmp;
        }

        bool operator==(const bond_iterator<MoleculeType> &other) const
        {
          return m_index == other.m_index;
        }

        bool operator!=(const bond_iterator<MoleculeType> &other) const
        {
          return m_index != other.m_index;
        }

      private:
        const MoleculeType *m_mol;
        Index m_index;
        mutable typename molecule_traits<MoleculeType>::bond_type m_bond;
    };

    template<typename T>
    void apply_permutation(std::vector<T> &elements, const std::vector<Index> &permutation)
    {
        std::vector<T> copy;
        for (std::size_t i = 0; i < permutation.size(); ++i)
          copy.push_back(elements[permutation[i]]);
        elements = copy;
    }

    /**
     * @brief Class representing a molecule.
     */
    template<template<typename> class AtomType, template<typename> class BondType>
    class MolImpl
    {
      public:
        typedef AtomType<MolImpl<AtomType, BondType> > atom_type;
        typedef BondType<MolImpl<AtomType, BondType> > bond_type;

        // iterators
        typedef impl::atom_iterator<MolImpl<AtomType, BondType> > atom_iter;
        typedef impl::bond_iterator<MolImpl<AtomType, BondType> > bond_iter;
        typedef typename std::vector<bond_type>::iterator incident_iter;
        typedef impl::nbr_iterator<MolImpl<AtomType, BondType> > nbr_iter;

        MolImpl()
        {
        }

        MolImpl(const MolImpl &other);

        MolImpl& operator=(const MolImpl &other);

        ~MolImpl()
        {
        }

        Size numAtoms() const
        {
          return m_element.size();
        }

        Size numBonds() const
        {
          return m_order.size();
        }

        iterator_pair<atom_iter> atoms() const
        {
          return make_iterator_pair(atom_iter(this, 0), atom_iter(this, m_element.size()));
        }

        iterator_pair<bond_iter> bonds() const
        {
          return make_iterator_pair(bond_iter(this, 0), bond_iter(this, m_order.size()));
        }

        atom_type atom(Index index) const
        {
          return atom_type(const_cast<MolImpl*>(this), index);
        }

        bond_type bond(Index index) const
        {
          return bond_type(const_cast<MolImpl*>(this), index);
        }

        static Index null_index()
        {
          return -1;
        }

        static atom_type null_atom()
        {
          return atom_type(0, -1);
        }

        static bond_type null_bond()
        {
          return bond_type(0, -1);
        }

        //atom_type addAtom();
        typename molecule_traits<MolImpl<AtomType, BondType> >::atom_type addAtom();
        void removeAtom(const atom_type &atom);

        //bond_type addBond(const atom_type &source, const atom_type &target);
        typename molecule_traits<MolImpl<AtomType, BondType> >::bond_type
          addBond(const atom_type &source, const atom_type &target);

        void removeBond(const bond_type &bond);

        void clear();

        void renumberAtoms(const std::vector<Index> &permutation);

      private:
        template<typename> friend class impl::AtomImpl;
        template<typename> friend class impl::BondImpl;

        // atoms
        std::vector<std::vector<bond_type> > m_adjList;
        std::vector<bool> m_atomAromatic;
        std::vector<unsigned char> m_element;
        std::vector<unsigned char> m_mass;
        std::vector<unsigned char> m_hydrogens;
        std::vector<signed char> m_charge;

        // bonds
        std::vector<Index> m_source;
        std::vector<Index> m_target;
        std::vector<bool> m_bondAromatic;
        std::vector<unsigned char> m_order;
    };

    template<template<typename> class AtomType, template<typename> class BondType>
    MolImpl<AtomType, BondType>::MolImpl(const MolImpl &other)
    {
      // copy m_adjList
      m_adjList.resize(other.m_adjList.size());
      for (std::size_t i = 0; i < other.m_adjList.size(); ++i)
        for (std::size_t j = 0; j < other.m_adjList[i].size(); ++j)
          m_adjList[i].push_back(bond_type(this, other.m_adjList[i][j].index()));

      // copy atom properties
      m_atomAromatic = other.m_atomAromatic;
      m_element = other.m_element;
      m_mass = other.m_mass;
      m_hydrogens = other.m_hydrogens;
      m_charge = other.m_charge;

      // copy bond properties
      m_source = other.m_source;
      m_target = other.m_target;
      m_bondAromatic = other.m_bondAromatic;
      m_order = other.m_order;
    }

    template<template<typename> class AtomType, template<typename> class BondType>
    MolImpl<AtomType, BondType>& MolImpl<AtomType, BondType>::operator=(const MolImpl &other)
    {
      // copy m_adjList
      m_adjList.clear();
      m_adjList.resize(other.m_adjList.size());
      for (std::size_t i = 0; i < other.m_adjList.size(); ++i)
        for (std::size_t j = 0; j < other.m_adjList[i].size(); ++j)
          m_adjList[i].push_back(bond_type(this, other.m_adjList[i][j].index()));

      // copy atom properties
      m_atomAromatic = other.m_atomAromatic;
      m_element = other.m_element;
      m_mass = other.m_mass;
      m_hydrogens = other.m_hydrogens;
      m_charge = other.m_charge;

      // copy bond properties
      m_source = other.m_source;
      m_target = other.m_target;
      m_bondAromatic = other.m_bondAromatic;
      m_order = other.m_order;

      return *this;
    }

    template<template<typename> class AtomType, template<typename> class BondType>
    typename molecule_traits<MolImpl<AtomType, BondType> >::atom_type MolImpl<AtomType, BondType>::addAtom()
    {
      Index index = m_element.size();
      m_adjList.resize(m_adjList.size() + 1);
      m_atomAromatic.resize(m_atomAromatic.size() + 1);
      m_element.resize(m_element.size() + 1);
      m_mass.resize(m_mass.size() + 1);
      m_hydrogens.resize(m_hydrogens.size() + 1);
      m_charge.resize(m_charge.size() + 1);

      return atom_type(this, index);
    }

    struct SortBondsByDecreasingIndex
    {
      template<typename BondType>
      bool operator()(const BondType &bond1, const BondType &bond2) const
      {
        return bond1.index() > bond2.index();
      }
    };

    template<template<typename> class AtomType, template<typename> class BondType>
    void MolImpl<AtomType, BondType>::removeAtom(const atom_type &atom)
    {
      Index index = atom.index();
      std::vector<bond_type> bonds = m_adjList[index];

      // remove properties
      m_adjList.erase(m_adjList.begin() + index);
      m_atomAromatic.erase(m_atomAromatic.begin() + index);
      m_element.erase(m_element.begin() + index);
      m_mass.erase(m_mass.begin() + index);
      m_hydrogens.erase(m_hydrogens.begin() + index);
      m_charge.erase(m_charge.begin() + index);

      // sort bonds by decreasing bond index so they can be correctly removed
      std::sort(bonds.begin(), bonds.end(), impl::SortBondsByDecreasingIndex());

      // update m_source & m_target
      for (std::size_t i = 0; i < m_source.size(); ++i)
        if (m_source[i] > index)
          --m_source[i];
      for (std::size_t i = 0; i < m_target.size(); ++i)
        if (m_target[i] > index)
          --m_target[i];

      // remove bonds to removed atom
      for (std::size_t i = 0; i < bonds.size(); ++i)
        removeBond(bonds[i]);
    }

    template<template<typename> class AtomType, template<typename> class BondType>
    typename molecule_traits<MolImpl<AtomType, BondType> >::bond_type
    MolImpl<AtomType, BondType>::addBond(const atom_type &source,
                                         const atom_type &target)
    {
      Index index = m_order.size();

      m_source.push_back(source.index());
      m_target.push_back(target.index());
      m_bondAromatic.resize(m_bondAromatic.size() + 1);
      m_order.resize(m_order.size() + 1, 1);

      bond_type bond(this, index);

      m_adjList[source.index()].push_back(bond);
      m_adjList[target.index()].push_back(bond);

      return bond;
    }

    template<template<typename> class AtomType, template<typename> class BondType>
    void MolImpl<AtomType, BondType>::removeBond(const bond_type &bond)
    {
      Index index = bond.index();

      // remove properties
      m_source.erase(m_source.begin() + index);
      m_target.erase(m_target.begin() + index);
      m_bondAromatic.erase(m_bondAromatic.begin() + index);
      m_order.erase(m_order.begin() + index);

      // update m_adjList
      for (std::size_t i = 0; i < m_adjList.size(); ++i) {
        int remove = -1;
        for (std::size_t j = 0; j < m_adjList[i].size(); ++j) {
          if (m_adjList[i][j].index() == index)
            remove = j;
          else if (m_adjList[i][j].index() > index)
            m_adjList[i][j] = bond_type(this, m_adjList[i][j].index() - 1);
        }
        if (remove != -1)
          m_adjList[i].erase(m_adjList[i].begin() + remove);
      }
    }

    template<template<typename> class AtomType, template<typename> class BondType>
    void MolImpl<AtomType, BondType>::clear()
    {
      m_adjList.clear();
      m_atomAromatic.clear();
      m_element.clear();
      m_mass.clear();
      m_hydrogens.clear();
      m_charge.clear();

      m_source.clear();
      m_target.clear();
      m_bondAromatic.clear();
      m_order.clear();
    }

    template<template<typename> class AtomType, template<typename> class BondType>
    void MolImpl<AtomType, BondType>::renumberAtoms(const std::vector<Index> &permutation)
    {
      assert(permutation.size() == m_adjList.size());
      impl::apply_permutation(m_adjList, permutation);
      impl::apply_permutation(m_atomAromatic, permutation);
      impl::apply_permutation(m_element, permutation);
      impl::apply_permutation(m_mass, permutation);
      impl::apply_permutation(m_hydrogens, permutation);
      impl::apply_permutation(m_charge, permutation);

      std::vector<Index> copy;
      for (std::size_t i = 0; i < m_source.size(); ++i)
        copy.push_back(index_of(permutation, m_source[i]));
      m_source = copy;

      copy.clear();
      for (std::size_t i = 0; i < m_target.size(); ++i)
        copy.push_back(index_of(permutation, m_target[i]));
      m_target = copy;
    }

  }

  typedef impl::MolImpl<impl::AtomImpl, impl::BondImpl> HeMol;
  typedef impl::AtomImpl<HeMol> HeAtom;
  typedef impl::BondImpl<HeMol> HeBond;

  //////////////////////////////////////////////////////////////////////////////
  //
  // Molecule
  //
  //////////////////////////////////////////////////////////////////////////////

  template<>
  inline void clear_molecule<HeMol>(HeMol &mol)
  {
    mol.clear();
  }

  /////////////////////////////
  //
  // HeAtoms
  //
  /////////////////////////////

  template<>
  inline Size num_atoms<HeMol>(const HeMol &mol)
  {
    return mol.numAtoms();
  }

  template<>
  inline iterator_pair<molecule_traits<HeMol>::atom_iter>
  get_atoms<HeMol>(const HeMol &mol)
  {
    return mol.atoms();
  }

  template<>
  inline molecule_traits<HeMol>::atom_type
  get_atom<HeMol>(const HeMol &mol, Index index)
  {
    return mol.atom(index);
  }

  template<>
  inline molecule_traits<HeMol>::atom_type add_atom<HeMol>(HeMol &mol)
  {
    return mol.addAtom();
  }

  template<>
  inline void remove_atom<HeMol>(HeMol &mol, molecule_traits<HeMol>::atom_type atom)
  {
    mol.removeAtom(atom);
  }

  /////////////////////////////
  //
  // HeBonds
  //
  /////////////////////////////

  template<>
  inline Size num_bonds<HeMol>(const HeMol &mol)
  {
    return mol.numBonds();
  }

  template<>
  inline iterator_pair<molecule_traits<HeMol>::bond_iter>
  get_bonds<HeMol>(const HeMol &mol)
  {
    return mol.bonds();
  }

  template<>
  inline molecule_traits<HeMol>::bond_type
  get_bond<HeMol>(const HeMol &mol, Index index)
  {
    return mol.bond(index);
  }

  template<>
  inline molecule_traits<HeMol>::bond_type add_bond<HeMol>(HeMol &mol,
      molecule_traits<HeMol>::atom_type source,
      molecule_traits<HeMol>::atom_type target)
  {
    return mol.addBond(source, target);
  }

  template<>
  inline void remove_bond<HeMol>(HeMol &mol, molecule_traits<HeMol>::bond_type bond)
  {
    mol.removeBond(bond);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // HeAtom
  //
  //////////////////////////////////////////////////////////////////////////////

  template<>
  inline Index get_index<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.index();
  }

  template<>
  inline iterator_pair<molecule_traits<HeMol>::incident_iter>
  get_bonds<HeMol>(const HeMol &mol, molecule_traits<HeMol>::atom_type atom)
  {
    return atom.bonds();
  }

  template<>
  inline iterator_pair<molecule_traits<HeMol>::nbr_iter>
  get_nbrs<HeMol>(const HeMol &mol, molecule_traits<HeMol>::atom_type atom)
  {
    return atom.nbrs();
  }

  template<>
  inline bool is_aromatic<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.isAromatic();
  }

  template<>
  inline void set_aromatic<HeMol>(HeMol &mol,
      molecule_traits<HeMol>::atom_type atom, bool value)
  {
    atom.setAromatic(value);
  }

  template<>
  inline int get_element<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.element();
  }

  template<>
  inline void set_element<HeMol>(HeMol &mol,
      molecule_traits<HeMol>::atom_type atom, int value)
  {
    atom.setElement(value);
  }

  template<>
  inline int get_mass<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.mass();
  }

  template<>
  inline void set_mass<HeMol>(HeMol &mol,
      molecule_traits<HeMol>::atom_type atom, int value)
  {
    atom.setMass(value);
  }

  template<>
  inline int get_degree<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.degree();
  }

  template<>
  inline int get_hydrogens<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.hydrogens();
  }

  template<>
  inline void set_hydrogens<HeMol>(HeMol &mol,
      molecule_traits<HeMol>::atom_type atom, int value)
  {
    atom.setHydrogens(value);
  }

  template<>
  inline int get_charge<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.charge();
  }

  template<>
  inline void set_charge<HeMol>(HeMol &mol,
      molecule_traits<HeMol>::atom_type atom, int value)
  {
    atom.setCharge(value);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // HeBond
  //
  //////////////////////////////////////////////////////////////////////////////

  template<>
  inline Index get_index<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::bond_type bond)
  {
    return bond.index();
  }

  template<>
  inline HeAtom get_source<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::bond_type bond)
  {
    return bond.source();
  }

  template<>
  inline HeAtom get_target<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::bond_type bond)
  {
    return bond.target();
  }

  template<>
  inline HeAtom get_other<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::bond_type bond,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return bond.other(atom);
  }

  template<>
  inline bool is_aromatic<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::bond_type bond)
  {
    return bond.isAromatic();
  }

  template<>
  inline void set_aromatic<HeMol>(HeMol &mol,
      molecule_traits<HeMol>::bond_type bond, bool value)
  {
    bond.setAromatic(value);
  }

  template<>
  inline int get_order<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::bond_type bond)
  {
    return bond.order();
  }

  template<>
  inline void set_order<HeMol>(HeMol &mol,
      molecule_traits<HeMol>::bond_type bond, int value)
  {
    bond.setOrder(value);
  }

  template<>
  inline molecule_traits<HeMol>::bond_type
  get_bond<HeMol>(const HeMol &mol,
      molecule_traits<HeMol>::atom_type source,
      molecule_traits<HeMol>::atom_type target)
  {
    molecule_traits<HeMol>::incident_iter bond, end_bonds;
    TIE(bond, end_bonds) = get_bonds(mol, source);
    for (; bond != end_bonds; ++bond)
      if (get_other(mol, *bond, source) == target)
        return *bond;
    return mol.null_bond();
  }

  /**
   * @brief STL output stream operator for HeMol.
   */
  std::ostream& operator<<(std::ostream &os, HeMol &mol);

  namespace impl {

    inline std::ostream& operator<<(std::ostream &os, const HeAtom &atom)
    {
      os << "HeAtom(" << atom.index() << ")";
      return os;
    }

    inline std::ostream& operator<<(std::ostream &os, const HeBond &bond)
    {
      os << "HeBond(" << bond.index() << ")";
      return os;
    }

  }

  //@endcond


}

#endif
