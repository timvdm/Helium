#ifndef HELIUM_SUBSTRUCTURE_H
#define HELIUM_SUBSTRUCTURE_H

#include <vector>
#include <algorithm>
#include <cassert>

#include "tie.h"
#include "hemol.h"

//#include <iostream>

namespace Helium {

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

  class Substructure
  {
    public:
      typedef HeMol molecule_type;

      typedef typename molecule_traits<molecule_type>::atom_type atom_type;
      typedef typename molecule_traits<molecule_type>::bond_type bond_type;

      typedef substructure_atom_iterator<Substructure> mol_atom_iter;
      typedef substructure_bond_iterator<Substructure> mol_bond_iter;
      typedef substructure_atom_bond_iterator<Substructure> atom_bond_iter;
      typedef substructure_nbr_iterator<Substructure> atom_atom_iter;

      typedef void atom_factory;
      typedef void bond_factory;

      Substructure(molecule_type *mol, const std::vector<bool> &atoms,                      
          const std::vector<bool> &bonds) : m_mol(mol),
          m_atoms(atoms), m_bonds(bonds)
      {
        assert(atoms.size() == num_atoms(mol));
        assert(bonds.size() == num_bonds(mol));
        m_numAtoms = std::count(atoms.begin(), atoms.end(), true);
        m_numBonds = std::count(bonds.begin(), bonds.end(), true);
        
        molecule_traits<molecule_type>::mol_atom_iter atom, end_atoms;
        tie(atom, end_atoms) = get_atoms(mol);
        unsigned int index = 0;
        for (; atom != end_atoms; ++atom)
          if (m_atoms[get_index(mol, *atom)])
            m_atomIndices.push_back(index++);
          else
            m_atomIndices.push_back(-1);
        
        molecule_traits<molecule_type>::mol_bond_iter bond, end_bonds;
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

      std::size_t numAtoms() const
      {
        return m_numAtoms;
      }

      std::size_t numBonds() const
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

      atom_type atom(unsigned int index) const
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

      bond_type bond(unsigned int index) const
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

      static unsigned int null_index()
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

      molecule_type* mol() const
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

      unsigned int atomIndex(atom_type atom) const
      {
        return m_atomIndices[get_index(m_mol, atom)];
      }

      unsigned int bondIndex(bond_type bond) const
      {
        return m_bondIndices[get_index(m_mol, bond)];
      }

    private:
      molecule_type *m_mol;
      unsigned int m_numAtoms;
      unsigned int m_numBonds;
      std::vector<bool> m_atoms;
      std::vector<bool> m_bonds;
      std::vector<unsigned int> m_atomIndices;
      std::vector<unsigned int> m_bondIndices;
  };



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

  template<>
  std::size_t num_atoms<Substructure>(const Substructure *mol)
  {
    return mol->numAtoms();
  }

  template<>
  std::pair<molecule_traits<Substructure>::mol_atom_iter, molecule_traits<Substructure>::mol_atom_iter>
  get_atoms<Substructure>(Substructure *mol)
  {
    return mol->atoms();
  }

  template<>
  typename molecule_traits<Substructure>::atom_type get_atom<Substructure>(const Substructure *mol, std::size_t index)
  {
    return mol->atom(index);
  }


  /////////////////////////////
  //
  // Bonds
  //
  /////////////////////////////

  template<>
  std::size_t num_bonds<Substructure>(const Substructure*mol)
  {
    return mol->numBonds();
  }

  template<>
  std::pair<molecule_traits<Substructure>::mol_bond_iter, molecule_traits<Substructure>::mol_bond_iter>
  get_bonds(Substructure*mol)
  {
    return mol->bonds();
  }

  template<>
  molecule_traits<Substructure>::bond_type get_bond(const Substructure*mol, std::size_t index)
  {
    return mol->bond(index);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Atom
  //
  //////////////////////////////////////////////////////////////////////////////

  template<>
  std::size_t get_index(const Substructure *mol, typename molecule_traits<Substructure>::atom_type atom)
  {
    return mol->atomIndex(atom);
  }

  template<>
  std::pair<molecule_traits<Substructure>::atom_bond_iter, molecule_traits<Substructure>::atom_bond_iter>
  get_bonds(const Substructure *mol, molecule_traits<Substructure>::atom_type atom)
  {
    molecule_traits<HeMol>::atom_bond_iter begin, end;
    tie(begin, end) = get_bonds(mol->mol(), atom);
    typedef molecule_traits<Substructure>::atom_bond_iter atom_bond_iter;
    return std::make_pair(atom_bond_iter(mol, begin, end), atom_bond_iter(mol, end, end));
  }

  template<>
  std::pair<typename molecule_traits<Substructure>::atom_atom_iter, typename molecule_traits<Substructure>::atom_atom_iter>
  get_nbrs(const Substructure *mol, typename molecule_traits<Substructure>::atom_type atom)
  {
    molecule_traits<Substructure>::atom_bond_iter begin, end;
    tie(begin, end) = get_bonds(mol, atom);
    typedef molecule_traits<Substructure>::atom_atom_iter atom_atom_iter;
    return std::make_pair(atom_atom_iter(atom, begin), atom_atom_iter(atom, end));
  }

  template<>
  bool is_aromatic(const Substructure *mol, typename molecule_traits<Substructure>::atom_type atom)
  {
    return atom->isAromatic();
  }

  template<>
  bool is_cyclic(const Substructure *mol, typename molecule_traits<Substructure>::atom_type atom)
  {
    return atom->isCyclic();
  }

  template<>
  int get_element(const Substructure *mol, typename molecule_traits<Substructure>::atom_type atom)
  {
    return atom->element();
  }

  template<>
  int get_mass(const Substructure *mol, typename molecule_traits<Substructure>::atom_type atom)
  {
    return atom->mass();
  }

  template<>
  int get_degree(const Substructure *mol, typename molecule_traits<Substructure>::atom_type atom)
  {
    molecule_traits<Substructure>::atom_bond_iter bond, end_bonds;
    tie(bond, end_bonds) = get_bonds(const_cast<Substructure*>(mol), atom);
    int d = 0;
    for (; bond != end_bonds; ++bond)
      ++d;
    return d;
  }

  template<>
  int num_hydrogens(const Substructure *mol, typename molecule_traits<Substructure>::atom_type atom)
  {
    return atom->hydrogens();
  }

  template<>
  int get_charge(const Substructure *mol, typename molecule_traits<Substructure>::atom_type atom)
  {
    return atom->charge();
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Bond
  //
  //////////////////////////////////////////////////////////////////////////////

  template<>
  std::size_t get_index(const Substructure *mol, typename molecule_traits<Substructure>::bond_type bond)
  {
    return mol->bondIndex(bond);
  }

  template<>
  HeAtom* get_source<Substructure>(const Substructure *mol, typename molecule_traits<Substructure>::bond_type bond)
  {
    return bond->source();
  }

  template<>
  HeAtom* get_target(const Substructure *mol, typename molecule_traits<Substructure>::bond_type bond)
  {
    return bond->target();
  }
  
  template<>
  molecule_traits<Substructure>::atom_type get_other(const Substructure*mol, 
                                                     molecule_traits<Substructure>::bond_type bond, 
                                                     molecule_traits<Substructure>::atom_type atom)
  {
    return bond->other(atom);
  }

  template<>
  bool is_aromatic(const Substructure *mol, typename molecule_traits<Substructure>::bond_type bond)
  {
    return bond->isAromatic();
  }

  template<>
  bool is_cyclic(const Substructure *mol, typename molecule_traits<Substructure>::bond_type bond)
  {
    return bond->isCyclic();
  }

  template<>
  bool get_order(const Substructure *mol, typename molecule_traits<Substructure>::bond_type bond)
  {
    return bond->order();
  }

  template<>
  molecule_traits<Substructure>::bond_type get_bond(const Substructure *mol, molecule_traits<Substructure>::atom_type source,
                                                                             molecule_traits<Substructure>::atom_type target)
  {
    molecule_traits<Substructure>::atom_bond_iter bond, end_bonds;
    tie(bond, end_bonds) = get_bonds(mol, source);
    for (; bond != end_bonds; ++bond)
      if (get_other(mol->mol(), *bond, source) == target)
        return *bond;
    return 0;
  }


}

#endif
