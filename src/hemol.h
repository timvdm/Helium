#ifndef HELIUM_HEMOL_H
#define HELIUM_HEMOL_H

#include "molecule.h"
#include "tie.h"

#include <vector>
#include <istream>
#include <algorithm>
#include <cassert>

namespace Helium {

  class HeBond;
  class HeMol;
 
  namespace impl {

    template<typename HeAtomType, typename HeBondType>
    class nbr_iterator
    {
      public:
        nbr_iterator()
        {
        }

        nbr_iterator(HeAtomType *atom, typename std::vector<HeBondType*>::iterator iter) : m_atom(atom), m_iter(iter)
        {
        }

        HeAtomType* operator*() const
        {
          return (*m_iter)->other(m_atom);
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
        HeAtomType *m_atom;
        typename std::vector<HeBondType*>::iterator m_iter;
    };

  }

  /**
   * @brief Class representing an atom in an HeMol.
   */
  class HeAtom
  {
    public:
      HeAtom(int index, bool aromatic, bool cyclic, int element, int mass,
           int hydrogens, int charge) : m_index(index), m_aromatic(aromatic),
           m_cyclic(cyclic), m_element(element), m_mass(mass),
           m_hydrogens(hydrogens), m_charge(charge)
      {
      }

      std::pair<std::vector<HeBond*>::iterator, std::vector<HeBond*>::iterator> bonds()
      {
        return std::make_pair(m_bonds.begin(), m_bonds.end());
      }

      std::pair<impl::nbr_iterator<HeAtom, HeBond>, impl::nbr_iterator<HeAtom, HeBond> > nbrs()
      {
        return std::make_pair(impl::nbr_iterator<HeAtom, HeBond>(this, m_bonds.begin()),
                              impl::nbr_iterator<HeAtom, HeBond>(this, m_bonds.end()));
      }

      Index index() const
      {
        return m_index;
      }

      bool isAromatic() const
      {
        return m_aromatic;
      }

      bool isCyclic() const
      {
        return m_cyclic;
      }

      int element() const
      {
        return m_element;
      }

      int mass() const
      {
        return m_mass;
      }

      int degree() const
      {
        return m_bonds.size();
      }

      int hydrogens() const
      {
        return m_hydrogens;
      }

      int charge() const
      {
        return m_charge;
      }

      void addBond(HeBond *bond)
      {
        m_bonds.push_back(bond);
      }

    private:
      friend class HeMol;

      void setIndex(Index index)
      {
        m_index = index;
      }

      std::vector<HeBond*> m_bonds;
      Index m_index;
      bool m_aromatic;
      bool m_cyclic;
      unsigned char m_element;
      unsigned char m_mass;
      unsigned char m_hydrogens;
      signed char m_charge;
  };

  class HeBond
  {
    public:
      HeBond(Index index, HeAtom *source, HeAtom *target, bool aromatic, bool cyclic,
           int order) : m_source(source), m_target(target), m_index(index),
           m_aromatic(aromatic), m_cyclic(cyclic), m_order(order)
      {
      }

      Index index() const
      {
        return m_index;
      }

      HeAtom* source() const
      {
        return m_source;
      }

      HeAtom* target() const
      {
        return m_target;
      }

      HeAtom* other(const HeAtom *atom) const
      {
        return m_source == atom ? m_target : m_source;
      }

      bool isAromatic() const
      {
        return m_aromatic;
      }

      bool isCyclic() const
      {
        return m_cyclic;
      }

      int order() const
      {
        return m_order;
      }

    private:
      HeAtom *m_source;
      HeAtom *m_target;
      Index m_index;
      bool m_aromatic;
      bool m_cyclic;
      unsigned char m_order;
  };

  template<typename MoleculeType>
  bool read_molecule(std::istream &is, MoleculeType &mol);

  /**
   * @brief Class representing a molecule.
   */
  class HeMol
  {
    public:
      typedef HeAtom* atom_type;
      typedef HeBond* bond_type;

      // iterators
      typedef std::vector<HeAtom*>::iterator mol_atom_iter;
      typedef std::vector<HeBond*>::iterator mol_bond_iter;
      typedef std::vector<HeBond*>::iterator atom_bond_iter;
      typedef impl::nbr_iterator<HeAtom, HeBond> atom_atom_iter;

      // constant iterators
      typedef std::vector<HeAtom*>::const_iterator const_atom_iter;
      typedef std::vector<HeBond*>::const_iterator const_bond_iter;

      ~HeMol()
      {
        for (std::size_t i = 0; i < m_atoms.size(); ++i)
          delete m_atoms[i];
        for (std::size_t i = 0; i < m_bonds.size(); ++i)
          delete m_bonds[i];
      }

      Size numHeAtoms() const
      {
        return m_atoms.size();
      }

      Size numHeBonds() const
      {
        return m_bonds.size();
      }

      std::pair<std::vector<HeAtom*>::iterator, std::vector<HeAtom*>::iterator> atoms()
      {
        return std::make_pair(m_atoms.begin(), m_atoms.end());
      }

      std::pair<const_atom_iter, const_atom_iter> atoms() const
      {
        return std::make_pair(m_atoms.begin(), m_atoms.end());
      }

      std::pair<std::vector<HeBond*>::iterator, std::vector<HeBond*>::iterator> bonds()
      {
        return std::make_pair(m_bonds.begin(), m_bonds.end());
      }

      std::pair<const_bond_iter, const_bond_iter> bonds() const
      {
        return std::make_pair(m_bonds.begin(), m_bonds.end());
      }

      HeAtom* atom(Index index) const
      {
        return m_atoms[index];
      }

      HeBond* bond(Index index) const
      {
        return m_bonds[index];
      }

      static Index null_index()
      {
        return -1;
      }

      static HeAtom* null_atom()
      {
        return 0;
      }

      static HeBond* null_bond()
      {
        return 0;
      }

      void renumberAtoms(const std::vector<Index> &permutation)
      {
        assert(permutation.size() == m_atoms.size());
        //assert(unique_elements(permutation) == m_atoms.size());
        //assert(*std::min_elment(permutations.begin(), permutation.end()) == 0);
        //assert(*std::max_elment(permutations.begin(), permutation.end()) == m_atoms.size() - 1);
        std::vector<HeAtom*> atoms;
        for (std::size_t i = 0; i < permutation.size(); ++i)
          atoms.push_back(m_atoms[permutation[i]]);
        m_atoms = atoms;

        for (std::size_t i = 0; i < m_atoms.size(); ++i)
          m_atoms[i]->setIndex(i);
      }

     private:
      template<typename MoleculeType> friend bool read_molecule(std::istream &is, MoleculeType &mol);

      std::vector<HeAtom*> m_atoms;
      std::vector<HeBond*> m_bonds;
  };

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
    return mol.numHeAtoms();
  }

  inline std::pair< molecule_traits<HeMol>::mol_atom_iter,  molecule_traits<HeMol>::mol_atom_iter>
  get_atoms(HeMol &mol)
  {
    return mol.atoms();
  }

  inline std::pair< molecule_traits<HeMol>::const_atom_iter,  molecule_traits<HeMol>::const_atom_iter>
  get_atoms(const HeMol &mol)
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
    return mol.numHeBonds();
  }

  inline std::pair< molecule_traits<HeMol>::mol_bond_iter,  molecule_traits<HeMol>::mol_bond_iter>
  get_bonds(HeMol &mol)
  {
    return mol.bonds();
  }

  inline std::pair< molecule_traits<HeMol>::const_bond_iter,  molecule_traits<HeMol>::const_bond_iter>
  get_bonds(const HeMol &mol)
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
    return atom->index();
  }

  inline std::pair< molecule_traits<HeMol>::atom_bond_iter,  molecule_traits<HeMol>::atom_bond_iter>
  get_bonds(const HeMol &mol,  molecule_traits<HeMol>::atom_type atom)
  {
    return atom->bonds();
  }

  inline std::pair< molecule_traits<HeMol>::atom_atom_iter,  molecule_traits<HeMol>::atom_atom_iter>
  get_nbrs(const HeMol &mol,  molecule_traits<HeMol>::atom_type atom)
  {
    return atom->nbrs();
  }

  inline bool is_aromatic(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom->isAromatic();
  }

  inline bool is_cyclic(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom->isCyclic();
  }

  inline int get_element(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom->element();
  }

  inline int get_mass(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom->mass();
  }

  inline int get_degree(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom->degree();
  }

  inline int num_hydrogens(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom->hydrogens();
  }

  inline int get_charge(const HeMol &mol, const  molecule_traits<HeMol>::atom_type atom)
  {
    return atom->charge();
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // HeBond
  //
  //////////////////////////////////////////////////////////////////////////////

  inline Index get_index(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond->index();
  }

  inline HeAtom* get_source(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond->source();
  }

  inline HeAtom* get_target(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond->target();
  }

  inline HeAtom* get_other(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond, const  molecule_traits<HeMol>::atom_type atom)
  {
    return bond->other(atom);
  }

  inline bool is_aromatic(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond->isAromatic();
  }

  inline bool is_cyclic(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond->isCyclic();
  }

  inline bool get_order(const HeMol &mol, const  molecule_traits<HeMol>::bond_type bond)
  {
    return bond->order();
  }

  inline  molecule_traits<HeMol>::bond_type get_bond(const HeMol &mol,  molecule_traits<HeMol>::atom_type source,
                                                                                molecule_traits<HeMol>::atom_type target)
  {
    molecule_traits<HeMol>::atom_bond_iter bond, end_bonds;
    tie(bond, end_bonds) = get_bonds(mol, source);
    for (; bond != end_bonds; ++bond)
      if (get_other(mol, *bond, source) == target)
        return *bond;
    return 0;
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
      os << "          \t" << (*atom)->index() << "\t" << (*atom)->element() << std::endl;
    
    os << "    Bonds:\tsource\ttarget\torder" << std::endl;
    molecule_traits<HeMol>::mol_bond_iter bond, end_bonds;
    tie(bond, end_bonds) = mol.bonds();
    for (; bond != end_bonds; ++bond)
      os << "          \t" << (*bond)->source()->index() << "\t" << (*bond)->target()->index() << "\t" << (*bond)->order() << std::endl;
    return os;
  }

}

#endif
