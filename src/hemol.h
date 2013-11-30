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

    template<typename HeAtomType, typename HeBondType>
    class nbr_iterator
    {
      public:
        nbr_iterator()
        {
        }

        nbr_iterator(HeAtomType atom, typename std::vector<HeBondType>::iterator iter) : m_atom(atom), m_iter(iter)
        {
        }

        HeAtomType operator*() const
        {
          return (*m_iter).other(m_atom);
        }

        nbr_iterator<HeAtomType, HeBondType>& operator++()
        {
          ++m_iter;
          return *this;
        }

        nbr_iterator<HeAtomType, HeBondType> operator++(int)
        {
          nbr_iterator<HeAtomType, HeBondType> tmp = *this;
          ++m_iter;
          return tmp;
        }

        bool operator!=(const nbr_iterator<HeAtomType, HeBondType> &other)
        {
          return m_iter != other.m_iter;
        }

      private:
        HeAtomType m_atom;
        typename std::vector<HeBondType>::iterator m_iter;
    };


    /**
     * @brief Class representing an atom in an HeMol.
     */
    template<typename HeMolType>
    class HeAtom
    {
      public:
        typedef typename molecule_traits<HeMolType>::atom_type atom_type;
        typedef typename molecule_traits<HeMolType>::bond_type bond_type;

        HeAtom(HeMolType *mol = 0, Index index = -1) : m_mol(mol), m_index(index)
        {
        }

        std::pair<typename std::vector<bond_type>::iterator, typename std::vector<bond_type>::iterator> bonds()
        {
          return std::make_pair(m_mol->m_adjList[m_index].begin(), m_mol->m_adjList[m_index].end());
        }

        std::pair<impl::nbr_iterator<atom_type, bond_type>, impl::nbr_iterator<atom_type, bond_type> > nbrs()
        {
          return std::make_pair(impl::nbr_iterator<atom_type, bond_type>(*this, m_mol->m_adjList[m_index].begin()),
                                impl::nbr_iterator<atom_type, bond_type>(*this, m_mol->m_adjList[m_index].end()));
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
        HeMolType *m_mol;
        Index m_index;
    };

    template<typename HeMolType>
    class HeBond
    {
      public:
        typedef typename molecule_traits<HeMolType>::atom_type atom_type;
        typedef typename molecule_traits<HeMolType>::bond_type bond_type;

        HeBond(HeMolType *mol = 0, Index index = -1) : m_mol(mol), m_index(index)
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
        HeMolType *m_mol;
        Index m_index;
    };

    template<typename HeMolType>
    class atom_iterator
    {
      public:
        atom_iterator()
        {
        }

        atom_iterator(const HeMolType *mol, Index index = 0) : m_mol(mol), m_index(index)
        {
        }

        HeAtom<HeMolType> operator*() const
        {
          return HeAtom<HeMolType>(const_cast<HeMolType*>(m_mol), m_index);
        }

        atom_iterator<HeMolType>& operator++()
        {
          ++m_index;
          return *this;
        }

        atom_iterator<HeMolType> operator++(int)
        {
          atom_iterator<HeMolType> tmp = *this;
          ++m_index;
          return tmp;
        }

        bool operator!=(const atom_iterator<HeMolType> &other)
        {
          return m_index != other.m_index;
        }

      private:
        const HeMolType *m_mol;
        Index m_index;
    };

    template<typename HeMolType>
    class bond_iterator
    {
      public:
        bond_iterator()
        {
        }

        bond_iterator(const HeMolType *mol, Index index = 0) : m_mol(mol), m_index(index)
        {
        }

        HeBond<HeMolType> operator*() const
        {
          return HeBond<HeMolType>(const_cast<HeMolType*>(m_mol), m_index);
        }

        bond_iterator<HeMolType>& operator++()
        {
          ++m_index;
          return *this;
        }

        bond_iterator<HeMolType> operator++(int)
        {
          bond_iterator<HeMolType> tmp = *this;
          ++m_index;
          return tmp;
        }

        bool operator!=(const bond_iterator<HeMolType> &other)
        {
          return m_index != other.m_index;
        }

      private:
        const HeMolType *m_mol;
        Index m_index;
    };

    template<typename T>
    void apply_permutation(std::vector<T> &elements, const std::vector<Index> &permutation)
    {
        std::vector<T> copy;
        for (std::size_t i = 0; i < permutation.size(); ++i)
          copy.push_back(elements[permutation[i]]);
        elements = copy;
    }

  }

  /**
   * @brief Class representing a molecule.
   */
  class HeMol
  {
    public:
      typedef impl::HeAtom<HeMol> atom_type;
      typedef impl::HeBond<HeMol> bond_type;

      // iterators
      typedef impl::atom_iterator<HeMol> atom_iter;
      typedef impl::bond_iterator<HeMol> bond_iter;
      typedef std::vector<bond_type>::iterator incident_iter;
      typedef impl::nbr_iterator<atom_type, bond_type> nbr_iter;

      HeMol()
      {
      }

      HeMol(const HeMol &other);

      HeMol& operator=(const HeMol &other);

      ~HeMol()
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

      std::pair<atom_iter, atom_iter> atoms() const
      {
        return std::make_pair(atom_iter(this, 0), atom_iter(this, m_element.size()));
      }

      std::pair<bond_iter, bond_iter> bonds() const
      {
        return std::make_pair(bond_iter(this, 0), bond_iter(this, m_order.size()));
      }

      atom_type atom(Index index) const
      {
        return atom_type(const_cast<HeMol*>(this), index);
      }

      bond_type bond(Index index) const
      {
        return bond_type(const_cast<HeMol*>(this), index);
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

      atom_type addAtom();

      bond_type addBond(const atom_type &source, const atom_type &target);

      void clear();

      void renumberAtoms(const std::vector<Index> &permutation);

     private:
      template<typename> friend class impl::HeAtom;
      template<typename> friend class impl::HeBond;

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

  typedef impl::HeAtom<HeMol> HeAtom;
  typedef impl::HeBond<HeMol> HeBond;

  //////////////////////////////////////////////////////////////////////////////
  //
  // Molecule
  //
  //////////////////////////////////////////////////////////////////////////////

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
  inline std::pair<molecule_traits<HeMol>::atom_iter,
                   molecule_traits<HeMol>::atom_iter>
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
  inline std::pair<molecule_traits<HeMol>::bond_iter,
                   molecule_traits<HeMol>::bond_iter>
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
  inline std::pair<molecule_traits<HeMol>::incident_iter,
                   molecule_traits<HeMol>::incident_iter>
  get_bonds<HeMol>(const HeMol &mol, molecule_traits<HeMol>::atom_type atom)
  {
    return atom.bonds();
  }

  template<>
  inline std::pair<molecule_traits<HeMol>::nbr_iter,
                   molecule_traits<HeMol>::nbr_iter>
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
  inline int get_element<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.element();
  }

  template<>
  inline int get_mass<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.mass();
  }

  template<>
  inline int get_degree<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.degree();
  }

  template<>
  inline int num_hydrogens<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.hydrogens();
  }

  template<>
  inline int get_charge<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::atom_type atom)
  {
    return atom.charge();
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
  inline int get_order<HeMol>(const HeMol &mol,
      const molecule_traits<HeMol>::bond_type bond)
  {
    return bond.order();
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

  //@endcond

  /**
   * @brief Read a SMILES string and construct a HeMol object from it.
   *
   * @param smiles The SMILES string.
   * @param mol The molecule output paramater.
   */
  void hemol_from_smiles(const std::string &smiles, HeMol &mol);

  /**
   * @brief Read a SMILES string and construct a HeMol object from it.
   *
   * @param smiles The SMILES string.
   *
   * @return The HeMol object for the SMILES.
   */
  HeMol hemol_from_smiles(const std::string &smiles);

}

#endif
