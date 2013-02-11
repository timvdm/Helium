#ifndef HELIUM_HEMOL_H
#define HELIUM_HEMOL_H

#include "molecule.h"
#include "tie.h"
#include "util/vector.h"

#include <vector>
#include <istream>
#include <algorithm>
#include <cassert>

namespace Helium {

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

        bool isCyclic() const
        {
          return m_mol->m_atomCyclic[m_index];
        }

        void setCyclic(bool value)
        {
          m_mol->m_atomCyclic[m_index] = value;
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

        bool isCyclic() const
        {
          return m_mol->m_bondCyclic[m_index];
        }

        void setCyclic(bool value)
        {
          m_mol->m_bondCyclic[m_index] = value;
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

        atom_iterator(HeMolType *mol, Index index = 0) : m_mol(mol), m_index(index)
        {
        }

        HeAtom<HeMolType> operator*() const
        {
          return HeAtom<HeMolType>(m_mol, m_index);
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
        HeMolType *m_mol;
        Index m_index;
    };

    template<typename HeMolType>
    class bond_iterator
    {
      public:
        bond_iterator()
        {
        }

        bond_iterator(HeMolType *mol, Index index = 0) : m_mol(mol), m_index(index)
        {
        }

        HeBond<HeMolType> operator*() const
        {
          return HeBond<HeMolType>(m_mol, m_index);
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
        HeMolType *m_mol;
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
      typedef impl::atom_iterator<HeMol> mol_atom_iter;
      typedef impl::bond_iterator<HeMol> mol_bond_iter;
      typedef std::vector<bond_type>::iterator atom_bond_iter;
      typedef impl::nbr_iterator<atom_type, bond_type> atom_atom_iter;

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

      std::pair<mol_atom_iter, mol_atom_iter> atoms()
      {
        return std::make_pair(mol_atom_iter(this, 0), mol_atom_iter(this, m_element.size()));
      }

      std::pair<mol_bond_iter, mol_bond_iter> bonds()
      {
        return std::make_pair(mol_bond_iter(this, 0), mol_bond_iter(this, m_order.size()));
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
      std::vector<bool> m_atomCyclic;
      std::vector<unsigned char> m_element;
      std::vector<unsigned char> m_mass;
      std::vector<unsigned char> m_hydrogens;
      std::vector<signed char> m_charge;

      // bonds
      std::vector<Index> m_source;
      std::vector<Index> m_target;
      std::vector<bool> m_bondAromatic;
      std::vector<bool> m_bondCyclic;
      std::vector<unsigned char> m_order;
  };

  typedef impl::HeAtom<HeMol> HeAtom;
  typedef impl::HeBond<HeMol> HeBond;

  //@cond dev

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

  inline Size num_atoms(const HeMol &mol)
  {
    return mol.numAtoms();
  }

  inline std::pair< molecule_traits<HeMol>::mol_atom_iter,  molecule_traits<HeMol>::mol_atom_iter>
  get_atoms(HeMol &mol)
  {
    return mol.atoms();
  }

  inline  molecule_traits<HeMol>::atom_type get_atom(const HeMol &mol, Index index)
  {
    return mol.atom(index);
  }


  /////////////////////////////
  //
  // HeBonds
  //
  /////////////////////////////

  inline Size num_bonds(const HeMol &mol)
  {
    return mol.numBonds();
  }

  inline std::pair< molecule_traits<HeMol>::mol_bond_iter,  molecule_traits<HeMol>::mol_bond_iter>
  get_bonds(HeMol &mol)
  {
    return mol.bonds();
  }

  inline  molecule_traits<HeMol>::bond_type get_bond(const HeMol &mol, Index index)
  {
    return mol.bond(index);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // HeAtom
  //
  //////////////////////////////////////////////////////////////////////////////

  inline Index get_index(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom.index();
  }

  inline std::pair< molecule_traits<HeMol>::atom_bond_iter,  molecule_traits<HeMol>::atom_bond_iter>
  get_bonds(const HeMol &mol,  molecule_traits<HeMol>::atom_type atom)
  {
    return atom.bonds();
  }

  inline std::pair< molecule_traits<HeMol>::atom_atom_iter,  molecule_traits<HeMol>::atom_atom_iter>
  get_nbrs(const HeMol &mol,  molecule_traits<HeMol>::atom_type atom)
  {
    return atom.nbrs();
  }

  inline bool is_aromatic(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom.isAromatic();
  }

  inline bool is_cyclic(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom.isCyclic();
  }

  inline int get_element(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom.element();
  }

  inline int get_mass(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom.mass();
  }

  inline int get_degree(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom.degree();
  }

  inline int num_hydrogens(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom.hydrogens();
  }

  inline int get_charge(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom.charge();
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // HeBond
  //
  //////////////////////////////////////////////////////////////////////////////

  inline Index get_index(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond.index();
  }

  inline HeAtom get_source(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond.source();
  }

  inline HeAtom get_target(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond.target();
  }

  inline HeAtom get_other(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond, const  molecule_traits<HeMol>::atom_type atom)
  {
    return bond.other(atom);
  }

  inline bool is_aromatic(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond.isAromatic();
  }

  inline bool is_cyclic(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond.isCyclic();
  }

  inline bool get_order(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond.order();
  }

  inline  molecule_traits<HeMol>::bond_type get_bond(const HeMol &mol,  molecule_traits<HeMol>::atom_type source,
                                                                                molecule_traits<HeMol>::atom_type target)
  {
    molecule_traits<HeMol>::atom_bond_iter bond, end_bonds;
    tie(bond, end_bonds) = get_bonds(mol, source);
    for (; bond != end_bonds; ++bond)
      if (get_other(mol, *bond, source) == target)
        return *bond;
    return mol.null_bond();
  }

  //@endcond

  /**
   * @brief STL output stream operator for HeMol.
   */
  inline std::ostream& operator<<(std::ostream &os, HeMol &mol)
  {
    os << "Molecule:" << std::endl;
    os << "    Atoms:\tindex\telement" << std::endl;
    molecule_traits<HeMol>::mol_atom_iter atom, end_atoms;
    tie(atom, end_atoms) = mol.atoms();
    for (; atom != end_atoms; ++atom)
      os << "          \t" << (*atom).index() << "\t" << (*atom).element() << std::endl;

    os << "    Bonds:\tsource\ttarget\torder" << std::endl;
    molecule_traits<HeMol>::mol_bond_iter bond, end_bonds;
    tie(bond, end_bonds) = mol.bonds();
    for (; bond != end_bonds; ++bond)
      os << "          \t" << (*bond).source().index() << "\t" << (*bond).target().index() << "\t" << (*bond).order() << std::endl;
    return os;
  }

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

#endif
