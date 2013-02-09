#ifndef HELIUM_SUBSTRUCTURE_H
#define HELIUM_SUBSTRUCTURE_H

#include <vector>
#include <algorithm>
#include <cassert>

#include "tie.h"
#include "hemol.h"

namespace Helium {

  //@cond dev

  namespace impl {

    template<typename SubstructureType>
    class substructure_atom_iterator
    {
      public:
        typedef typename SubstructureType::molecule_type molecule_type;
        typedef typename molecule_traits<molecule_type>::atom_type atom_type;
        typedef typename molecule_traits<molecule_type>::mol_atom_iter mol_atom_iter;

        substructure_atom_iterator()
        {
        }

        substructure_atom_iterator(SubstructureType *substructure, const mol_atom_iter &iter, const mol_atom_iter &end)
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
        SubstructureType *m_substructure;
        mol_atom_iter m_iter;
        mol_atom_iter m_end;
    };

    template<typename SubstructureType>
    class substructure_bond_iterator
    {
        public:
          typedef typename SubstructureType::molecule_type molecule_type;
          typedef typename molecule_traits<molecule_type>::bond_type bond_type;
          typedef typename molecule_traits<molecule_type>::mol_bond_iter mol_bond_iter;

          substructure_bond_iterator()
          {
          }

          substructure_bond_iterator(SubstructureType *substructure, const mol_bond_iter &iter, const mol_bond_iter &end)
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
          SubstructureType *m_substructure;
          mol_bond_iter m_iter;
          mol_bond_iter m_end;
    };

    template<typename SubstructureType>
    class substructure_atom_bond_iterator
    {
      public:
        typedef typename SubstructureType::molecule_type molecule_type;
        typedef typename molecule_traits<molecule_type>::bond_type bond_type;
        typedef typename molecule_traits<molecule_type>::atom_bond_iter atom_bond_iter;

        substructure_atom_bond_iterator()
        {
        }

        substructure_atom_bond_iterator(const SubstructureType *substructure, const atom_bond_iter &iter, const atom_bond_iter &end)
          : m_substructure(substructure), m_iter(iter), m_end(end)
        {
          while (m_iter != m_end && m_substructure->isHidden(*m_iter))
            ++m_iter;
        }

        bond_type operator*() const
        {
          return *m_iter;
        }

        substructure_atom_bond_iterator<SubstructureType>& operator++()
        {
          ++m_iter;
          while (m_iter != m_end && m_substructure->isHidden(*m_iter))
            ++m_iter;
          return *this;
        }

        substructure_atom_bond_iterator<SubstructureType> operator++(int)
        {
          substructure_atom_bond_iterator<SubstructureType> tmp = *this;
          ++m_iter;
          while (m_iter != m_end && m_substructure->isHidden(*m_iter))
            ++m_iter;
          return tmp;
        }

        bool operator!=(const substructure_atom_bond_iterator<SubstructureType> &other)
        {
          return m_iter != other.m_iter;
        }

      private:
        const SubstructureType *m_substructure;
        atom_bond_iter m_iter;
        atom_bond_iter m_end;
    };

    template<typename SubstructureType>
    class substructure_nbr_iterator
    {
        typedef typename molecule_traits<SubstructureType>::atom_type atom_type;
        typedef typename molecule_traits<SubstructureType>::bond_type bond_type;
        typedef typename molecule_traits<SubstructureType>::atom_bond_iter atom_bond_iter;
      public:
        substructure_nbr_iterator()
        {
        }

        substructure_nbr_iterator(atom_type atom, atom_bond_iter iter) : m_atom(atom), m_iter(iter)
        {
        }

        atom_type operator*() const
        {
          return (*m_iter)->other(m_atom);
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
        atom_bond_iter m_iter;
    };

  }

  //@endcond

  template<typename MoleculeType>
  class Substructure
  {
    public:
      typedef MoleculeType molecule_type;

      typedef typename molecule_traits<molecule_type>::atom_type atom_type;
      typedef typename molecule_traits<molecule_type>::bond_type bond_type;

      typedef impl::substructure_atom_iterator<Substructure> mol_atom_iter;
      typedef impl::substructure_bond_iterator<Substructure> mol_bond_iter;
      typedef impl::substructure_atom_bond_iterator<Substructure> atom_bond_iter;
      typedef impl::substructure_nbr_iterator<Substructure> atom_atom_iter;
      
      typedef impl::substructure_atom_iterator<Substructure> const_atom_iter;
      typedef impl::substructure_bond_iterator<Substructure> const_bond_iter;

      Substructure(molecule_type &mol, const std::vector<bool> &atoms,                      
          const std::vector<bool> &bonds) : m_mol(mol),
          m_atoms(atoms), m_bonds(bonds)
      {
        assert(atoms.size() == num_atoms(mol));
        assert(bonds.size() == num_bonds(mol));
        m_numAtoms = std::count(atoms.begin(), atoms.end(), true);
        m_numBonds = std::count(bonds.begin(), bonds.end(), true);
        
        typename molecule_traits<molecule_type>::mol_atom_iter atom, end_atoms;
        tie(atom, end_atoms) = get_atoms(mol);
        Index index = 0;
        for (; atom != end_atoms; ++atom)
          if (m_atoms[get_index(mol, *atom)])
            m_atomIndices.push_back(index++);
          else
            m_atomIndices.push_back(-1);
        
        typename molecule_traits<molecule_type>::mol_bond_iter bond, end_bonds;
        tie(bond, end_bonds) = get_bonds(mol);
        index = 0;
        for (; bond != end_bonds; ++bond)
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

      std::pair<mol_atom_iter, mol_atom_iter> atoms()
      {
        typename molecule_traits<molecule_type>::mol_atom_iter begin, end;
        tie(begin, end) = get_atoms(m_mol);
        return std::make_pair(mol_atom_iter(this, begin, end), mol_atom_iter(this, end, end));
      }

      std::pair<mol_bond_iter, mol_bond_iter> bonds()
      {
        typename molecule_traits<molecule_type>::mol_bond_iter begin, end;
        tie(begin, end) = get_bonds(m_mol);
        return std::make_pair(mol_bond_iter(this, begin, end), mol_bond_iter(this, end, end));
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

      molecule_type& mol() const
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

      Index atomIndex(atom_type atom) const
      {
        return m_atomIndices[get_index(m_mol, atom)];
      }

      Index bondIndex(bond_type bond) const
      {
        return m_bondIndices[get_index(m_mol, bond)];
      }

    private:
      molecule_type &m_mol;
      Size m_numAtoms;
      Size m_numBonds;
      std::vector<bool> m_atoms;
      std::vector<bool> m_bonds;
      std::vector<Index> m_atomIndices;
      std::vector<Index> m_bondIndices;
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

  template<typename SubstructureType>
  Size num_atoms(const SubstructureType &mol)
  {
    return mol.numAtoms();
  }

  template<typename SubstructureType>
  std::pair<typename molecule_traits<SubstructureType>::mol_atom_iter, typename molecule_traits<SubstructureType>::mol_atom_iter>
  get_atoms(SubstructureType &mol)
  {
    return mol.atoms();
  }

  template<typename SubstructureType>
  typename molecule_traits<SubstructureType>::atom_type get_atom(const SubstructureType &mol, Index index)
  {
    return mol.atom(index);
  }


  /////////////////////////////
  //
  // Bonds
  //
  /////////////////////////////

  template<typename SubstructureType>
  Size num_bonds(const SubstructureType &mol)
  {
    return mol.numBonds();
  }

  template<typename SubstructureType>
  std::pair<typename molecule_traits<SubstructureType>::mol_bond_iter, typename molecule_traits<SubstructureType>::mol_bond_iter>
  get_bonds(SubstructureType &mol)
  {
    return mol.bonds();
  }

  template<typename SubstructureType>
  typename molecule_traits<SubstructureType>::bond_type get_bond(const SubstructureType &mol, Index index)
  {
    return mol.bond(index);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Atom
  //
  //////////////////////////////////////////////////////////////////////////////

  template<typename SubstructureType>
  Index get_index(const SubstructureType &mol, typename molecule_traits<SubstructureType>::atom_type atom)
  {
    return mol.atomIndex(atom);
  }

  template<typename SubstructureType>
  std::pair<typename molecule_traits<SubstructureType>::atom_bond_iter, typename molecule_traits<SubstructureType>::atom_bond_iter>
  get_bonds(const SubstructureType &mol, typename molecule_traits<SubstructureType>::atom_type atom)
  {
    molecule_traits<HeMol>::atom_bond_iter begin, end;
    tie(begin, end) = get_bonds(mol.mol(), atom);
    typedef typename molecule_traits<SubstructureType>::atom_bond_iter atom_bond_iter;
    return std::make_pair(atom_bond_iter(&mol, begin, end), atom_bond_iter(&mol, end, end));
  }

  template<typename SubstructureType>
  std::pair<typename molecule_traits<SubstructureType>::atom_atom_iter, typename molecule_traits<SubstructureType>::atom_atom_iter>
  get_nbrs(const SubstructureType &mol, typename molecule_traits<SubstructureType>::atom_type atom)
  {
    typename molecule_traits<SubstructureType>::atom_bond_iter begin, end;
    tie(begin, end) = get_bonds(mol, atom);
    typedef typename molecule_traits<SubstructureType>::atom_atom_iter atom_atom_iter;
    return std::make_pair(atom_atom_iter(atom, begin), atom_atom_iter(atom, end));
  }

  template<typename SubstructureType>
  bool is_aromatic(const SubstructureType &mol, typename molecule_traits<SubstructureType>::atom_type atom)
  {
    return atom->isAromatic();
  }

  template<typename SubstructureType>
  bool is_cyclic(const SubstructureType &mol, typename molecule_traits<SubstructureType>::atom_type atom)
  {
    return atom->isCyclic();
  }

  template<typename SubstructureType>
  int get_element(const SubstructureType &mol, typename molecule_traits<SubstructureType>::atom_type atom)
  {
    return atom->element();
  }

  template<typename SubstructureType>
  int get_mass(const SubstructureType &mol, typename molecule_traits<SubstructureType>::atom_type atom)
  {
    return atom->mass();
  }

  template<typename SubstructureType>
  int get_degree(const SubstructureType &mol, typename molecule_traits<SubstructureType>::atom_type atom)
  {
    typename molecule_traits<SubstructureType>::atom_bond_iter bond, end_bonds;
    tie(bond, end_bonds) = get_bonds(const_cast<SubstructureType&>(mol), atom);
    int d = 0;
    for (; bond != end_bonds; ++bond)
      ++d;
    return d;
  }

  template<typename SubstructureType>
  int num_hydrogens(const SubstructureType &mol, typename molecule_traits<SubstructureType>::atom_type atom)
  {
    return atom->hydrogens();
  }

  template<typename SubstructureType>
  int get_charge(const SubstructureType &mol, typename molecule_traits<SubstructureType>::atom_type atom)
  {
    return atom->charge();
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Bond
  //
  //////////////////////////////////////////////////////////////////////////////

  template<typename SubstructureType>
  Index get_index(const SubstructureType &mol, typename molecule_traits<SubstructureType>::bond_type bond)
  {
    return mol.bondIndex(bond);
  }

  template<typename SubstructureType>
  typename molecule_traits<SubstructureType>::atom_type get_source(const SubstructureType &mol, typename molecule_traits<SubstructureType>::bond_type bond)
  {
    return bond->source();
  }

  template<typename SubstructureType>
  typename molecule_traits<SubstructureType>::atom_type get_target(const SubstructureType &mol, typename molecule_traits<SubstructureType>::bond_type bond)
  {
    return bond->target();
  }
  
  template<typename SubstructureType>
  typename molecule_traits<SubstructureType>::atom_type get_other(const SubstructureType &mol, 
                                                     typename molecule_traits<SubstructureType>::bond_type bond, 
                                                     typename molecule_traits<SubstructureType>::atom_type atom)
  {
    return bond->other(atom);
  }

  template<typename SubstructureType>
  bool is_aromatic(const SubstructureType &mol, typename molecule_traits<SubstructureType>::bond_type bond)
  {
    return bond->isAromatic();
  }

  template<typename SubstructureType>
  bool is_cyclic(const SubstructureType &mol, typename molecule_traits<SubstructureType>::bond_type bond)
  {
    return bond->isCyclic();
  }

  template<typename SubstructureType>
  bool get_order(const SubstructureType &mol, typename molecule_traits<SubstructureType>::bond_type bond)
  {
    return bond->order();
  }

  template<typename SubstructureType>
  typename molecule_traits<SubstructureType>::bond_type get_bond(const SubstructureType &mol, typename molecule_traits<SubstructureType>::atom_type source,
                                                                             typename molecule_traits<SubstructureType>::atom_type target)
  {
    typename molecule_traits<SubstructureType>::atom_bond_iter bond, end_bonds;
    tie(bond, end_bonds) = get_bonds(mol, source);
    for (; bond != end_bonds; ++bond)
      if (get_other(mol.mol(), *bond, source) == target)
        return *bond;
    return 0;
  }

  //@endcond

}

#endif
