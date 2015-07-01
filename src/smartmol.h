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
#ifndef HELIUM_SMARTMOL_H
#define HELIUM_SMARTMOL_H

#include <Helium/hemol.h>
#include <Helium/util.h>
#include <stdexcept>
#include <map>

namespace Helium {

  /**
   * @file smartmol.h
   * @brief SmartMol implementation.
   */

  //@cond DEV

  namespace impl {

    template<typename MoleculeType>
    class incident_iterator : public std::iterator<std::forward_iterator_tag,
        typename molecule_traits<MoleculeType>::bond_type>
    {
      public:
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

        incident_iterator()
        {
        }

        incident_iterator(MoleculeType *mol, typename std::vector<HeMol::bond_type>::iterator iter)
          : m_mol(mol), m_iter(iter)
        {
        }

        bond_type& operator*() const
        {
          m_bond = bond_type(m_mol, (*m_iter).index());
          return m_bond;
        }

        incident_iterator<MoleculeType>& operator++()
        {
          ++m_iter;
          return *this;
        }

        incident_iterator<MoleculeType> operator++(int)
        {
          incident_iterator<MoleculeType> tmp = *this;
          ++m_iter;
          return tmp;
        }

        bool operator==(const incident_iterator<MoleculeType> &other) const
        {
          return m_iter == other.m_iter;
        }

        bool operator!=(const incident_iterator<MoleculeType> &other) const
        {
          return m_iter != other.m_iter;
        }

      private:
        MoleculeType *m_mol;
        typename std::vector<HeMol::bond_type>::iterator m_iter;
        mutable bond_type m_bond;
    };

    template<typename MoleculeType>
    class nbr_iterator_wrapper : public std::iterator<std::forward_iterator_tag,
        typename molecule_traits<MoleculeType>::atom_type>
    {
      public:
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

        nbr_iterator_wrapper()
        {
        }

        nbr_iterator_wrapper(MoleculeType *mol, const HeMol::atom_type &atom,
            typename std::vector<HeMol::bond_type>::iterator iter)
          : m_mol(mol), m_atom(atom), m_iter(iter)
        {
        }

        atom_type& operator*() const
        {
          HeMol::atom_type nbr = (*m_iter).other(m_atom);
          m_nbr = atom_type(m_mol, nbr.index());
          return m_nbr;
        }

        nbr_iterator_wrapper<MoleculeType>& operator++()
        {
          ++m_iter;
          return *this;
        }

        nbr_iterator_wrapper<MoleculeType> operator++(int)
        {
          nbr_iterator_wrapper<MoleculeType> tmp = *this;
          ++m_iter;
          return tmp;
        }

        bool operator==(const nbr_iterator_wrapper<MoleculeType> &other) const
        {
          return m_iter == other.m_iter;
        }

        bool operator!=(const nbr_iterator_wrapper<MoleculeType> &other) const
        {
          return m_iter != other.m_iter;
        }

      private:
        MoleculeType *m_mol;
        HeMol::atom_type m_atom;
        typename std::vector<HeMol::bond_type>::iterator m_iter;
        mutable atom_type m_nbr;
    };

  }

  template<typename MoleculeType>
  class SmartAtom
  {
    public:
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      SmartAtom(MoleculeType *mol = 0, Index index = -1) : m_mol(mol), m_index(index)
      {
      }

      iterator_pair<typename molecule_traits<MoleculeType>::incident_iter> bonds()
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        iterator_pair<HeMol::incident_iter> iters = get_bonds(m_mol->m_mol, atom);
        return make_iterator_pair(impl::incident_iterator<MoleculeType>(m_mol, iters.begin()),
                                  impl::incident_iterator<MoleculeType>(m_mol, iters.end()));
      }

      iterator_pair<impl::nbr_iterator_wrapper<MoleculeType> > nbrs()
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        iterator_pair<HeMol::incident_iter> iters = get_bonds(m_mol->m_mol, atom);
        return make_iterator_pair(impl::nbr_iterator_wrapper<MoleculeType>(m_mol, atom, iters.begin()),
                                  impl::nbr_iterator_wrapper<MoleculeType>(m_mol, atom, iters.end()));
      }

      MoleculeType* mol() const
      {
        return m_mol;
      }

      Index index() const
      {
        return m_index;
      }

      bool isAromatic() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return is_aromatic(m_mol->m_mol, atom);
      }

      void setAromatic(bool value)
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        set_aromatic(m_mol->m_mol, atom, value);
      }

      int element() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return get_element(m_mol->m_mol, atom);
      }

      void setElement(int value)
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        set_element(m_mol->m_mol, atom, value);
      }

      int mass() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return get_mass(m_mol->m_mol, atom);
      }

      void setMass(int value)
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        set_mass(m_mol->m_mol, atom, value);
      }

      int degree() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return get_degree(m_mol->m_mol, atom);
      }

      int hydrogens() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return get_hydrogens(m_mol->m_mol, atom);
      }

      void setHydrogens(int value)
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        set_hydrogens(m_mol->m_mol, atom, value);
      }

      int charge() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return get_charge(m_mol->m_mol, atom);
      }

      void setCharge(int value)
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        set_charge(m_mol->m_mol, atom, value);
      }

      bool isHydrogen() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return is_hydrogen(m_mol->m_mol, atom);
      }

      bool isCarbon() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return is_carbon(m_mol->m_mol, atom);
      }

      bool isNitrogen() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return is_nitrogen(m_mol->m_mol, atom);
      }

      bool isOxygen() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return is_oxygen(m_mol->m_mol, atom);
      }

      bool isPhosphorus() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return is_phosphorus(m_mol->m_mol, atom);
      }

      bool isSulfur() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return is_sulfur(m_mol->m_mol, atom);
      }

      int heavyDegree() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return get_heavy_degree(m_mol->m_mol, atom);
      }

      int boSum() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return get_bosum(m_mol->m_mol, atom);
      }

      int valence() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return get_valence(m_mol->m_mol, atom);
      }

      int connectivity() const
      {
        HeMol::atom_type atom = get_atom(m_mol->m_mol, m_index);
        return get_connectivity(m_mol->m_mol, atom);
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

    private:
      MoleculeType *m_mol;
      Index m_index;
  };

  template<typename MoleculeType>
  class SmartBond
  {
    public:
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      SmartBond(MoleculeType *mol = 0, Index index = -1) : m_mol(mol), m_index(index)
      {
      }

      Index index() const
      {
        return m_index;
      }

      MoleculeType* mol() const
      {
        return m_mol;
      }

      atom_type source() const
      {
        HeMol::bond_type bond = get_bond(m_mol->m_mol, m_index);
        HeMol::atom_type source = get_source(m_mol->m_mol, bond);
        return atom_type(m_mol, source.index());
      }

      atom_type target() const
      {
        HeMol::bond_type bond = get_bond(m_mol->m_mol, m_index);
        HeMol::atom_type target = get_target(m_mol->m_mol, bond);
        return atom_type(m_mol, target.index());
      }

      atom_type other(const atom_type &atom) const
      {
        HeMol::bond_type bond = get_bond(m_mol->m_mol, m_index);
        HeMol::atom_type other = get_other(m_mol->m_mol, bond, get_atom(m_mol->m_mol, atom.index()));
        return atom_type(m_mol, get_index(m_mol->m_mol, other));
      }

      bool isAromatic() const
      {
        HeMol::bond_type bond = get_bond(m_mol->m_mol, m_index);
        return is_aromatic(m_mol->m_mol, bond);
      }

      void setAromatic(bool value)
      {
        HeMol::bond_type bond = get_bond(m_mol->m_mol, m_index);
        set_aromatic(m_mol->m_mol, bond, value);
      }

      int order() const
      {
        HeMol::bond_type bond = get_bond(m_mol->m_mol, m_index);
        return get_order(m_mol->m_mol, bond);
      }

      void setOrder(int value)
      {
        HeMol::bond_type bond = get_bond(m_mol->m_mol, m_index);
        set_order(m_mol->m_mol, bond, value);
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

    private:
      MoleculeType *m_mol;
      Index m_index;
  };

  /*
  class Attribute
  {
    public:
      enum Type {
        Null,
        Bool,
        Integer,
        Real,
        String,
        Pointer
      };

      Attribute() : m_type(Null), m_integer(0)
      {
      }

      ~Attribute()
      {
        if (m_type == String)
          delete m_string;
        if (m_type == Pointer)
          delete m_ptr;
      }

      Type type() const
      {
        return m_type;
      }

      bool toBool() const
      {
        assert(m_type == Bool);
        return m_bool;
      }

      long toInteger() const
      {
        assert(m_type == Integer);
        return m_integer;
      }

      double toReal() const
      {
        assert(m_type == Real);
        return m_real;
      }

      std::string toString() const
      {
        assert(m_type == String);
        return *m_string;
      }

      void* toPointer() const
      {
        assert(m_type == Pointer);
        return m_ptr;
      }

      void setBool(bool value)
      {
        m_type = Bool;
        m_bool = value;
      }

      void setInteger(long value)
      {
        m_type = Integer;
        m_integer = value;
      }

      void setReal(double value)
      {
        m_type = Real;
        m_real = value;
      }

      void setString(const std::string &value)
      {
        m_type = String;
        m_string = new std::string(value);
      }

      void setPointer(void *value)
      {
        m_ptr = value;
      }

    private:
      Type m_type;
      union {
        bool m_bool;
        long m_integer;
        double m_real;
        std::string *m_string;
        void *m_ptr;
      };
  };
  */

  class SmartMol;

  class SmartAttribute
  {
    public:
      enum ChangeType {
        ChangeElement,
        ChangeMass,
        ChangeCharge,
        ChangeHydrogens,
        ChangeAromaticity,
        ChangeOrder
      };

      SmartAttribute(const SmartMol &mol) : m_mol(mol)
      {
      }

      virtual ~SmartAttribute()
      {
      }

      const SmartMol& molecule() const
      {
        return m_mol;
      }

      virtual std::string name() const = 0;

      virtual void addAtom(Index index)
      {
      }

      virtual void addBond(Index index)
      {
      }

      virtual void removeAtom(Index index)
      {
      }

      virtual void removeBond(Index index)
      {
      }

      virtual void atomChange(Index index, ChangeType type, int oldValue)
      {
      }

      virtual void bondChange(Index index, ChangeType type, int oldValue)
      {
      }

      virtual void clear()
      {
      }

    private:
      const SmartMol &m_mol;
  };

  class SmartMol
  {
    public:
      typedef SmartAtom<SmartMol> atom_type;
      typedef SmartBond<SmartMol> bond_type;

      typedef impl::atom_iterator<SmartMol> atom_iter;
      typedef impl::bond_iterator<SmartMol> bond_iter;
      typedef impl::incident_iterator<SmartMol> incident_iter;
      typedef impl::nbr_iterator_wrapper<SmartMol> nbr_iter;

      SmartMol()
      {
      }

      SmartMol(const SmartMol &other)
      {
        m_mol = other.m_mol;
        m_attributes = other.m_attributes;
      }

      SmartMol(const SmartMol &source, const std::vector<bool> &atoms,
          const std::vector<bool> &bonds, bool adjustHydrogens = true)
      {
        if (atoms.size() != source.numAtoms())
          throw std::runtime_error("atoms parameter does not have the correct size");
        if (bonds.size() != source.numBonds())
          throw std::runtime_error("bonds parameter does not have the correct size");

        make_substructure(m_mol, source, atoms, bonds);
        m_attributes = source.m_attributes;
      }

      Size numAtoms() const
      {
        return num_atoms(m_mol);
      }

      Size numBonds() const
      {
        return num_bonds(m_mol);
      }

      iterator_pair<atom_iter> atoms() const
      {
        return make_iterator_pair(atom_iter(this, 0), atom_iter(this, m_mol.numAtoms()));
      }

      iterator_pair<bond_iter> bonds() const
      {
        return make_iterator_pair(bond_iter(this, 0), bond_iter(this, m_mol.numBonds()));
      }

      atom_type atom(Index index) const
      {
        if (index >= num_atoms(m_mol))
          throw std::runtime_error("Invalid atom index");
        return atom_type(const_cast<SmartMol*>(this), index);
      }

      bond_type bond(Index index) const
      {
        if (index >= num_bonds(m_mol))
          throw std::runtime_error("Invalid bond index");
        return bond_type(const_cast<SmartMol*>(this), index);
      }

      inline bond_type bond(const atom_type &source, const atom_type &target) const;

      atom_type addAtom()
      {
        // add a new bond to the molecule
        molecule_traits<HeMol>::atom_type atom = add_atom(m_mol);

        // invoke lazy attribute callbacks
        std::map<std::string, SmartAttribute*>::iterator attr;
        for (attr = m_attributes.begin(); attr != m_attributes.end(); ++attr)
          attr->second->addAtom(atom.index());

        // return the newly created atom
        return atom_type(this, get_index(m_mol, atom));
      }

      void removeAtom(const atom_type &atom)
      {
        // invoke lazy attribute callbacks
        std::map<std::string, SmartAttribute*>::iterator attr;
        for (attr = m_attributes.begin(); attr != m_attributes.end(); ++attr)
          attr->second->removeAtom(atom.index());

        // actually remove the atom
        remove_atom(m_mol, get_atom(m_mol, atom.index()));
      }

      bond_type addBond(const atom_type &source, const atom_type &target)
      {
        // add a new bond to the molecule
        molecule_traits<HeMol>::bond_type bond = add_bond(m_mol, get_atom(m_mol, source.index()), get_atom(m_mol, target.index()));

        // invoke lazy attribute callbacks
        std::map<std::string, SmartAttribute*>::iterator attr;
        for (attr = m_attributes.begin(); attr != m_attributes.end(); ++attr)
          attr->second->addBond(bond.index());

        // return the newly created bond
        return bond_type(this, get_index(m_mol, bond));
      }

      void removeBond(const bond_type &bond)
      {
        // invoke lazy attribute callbacks
        std::map<std::string, SmartAttribute*>::iterator attr;
        for (attr = m_attributes.begin(); attr != m_attributes.end(); ++attr)
          attr->second->removeBond(bond.index());

        // actually remove the bond
        remove_bond(m_mol, get_bond(m_mol, bond.index()));
      }

      void clear()
      {
        // invoke lazy attribute callbacks
        std::map<std::string, SmartAttribute*>::iterator attr;
        for (attr = m_attributes.begin(); attr != m_attributes.end(); ++attr)
          attr->second->clear();

        // actually clear the molecule
        clear_molecule(m_mol);
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

      void addAttribute(SmartAttribute *attribute)
      {
        // check if there is already an attribute with this name
        std::map<std::string, SmartAttribute*>::const_iterator attr = m_attributes.find(attribute->name());
        if (attr != m_attributes.end())
          delete attr->second;
        // set the new attribute
        m_attributes[attribute->name()] = attribute;
      }

      SmartAttribute* attribute(const std::string &name) const
      {
        std::map<std::string, SmartAttribute*>::const_iterator attr = m_attributes.find(name);
        if (attr != m_attributes.end())
          return attr->second;
        return 0;
      }

    private:
      template<typename> friend class SmartAtom;
      template<typename> friend class SmartBond;

      HeMol m_mol;
      std::map<std::string, SmartAttribute*> m_attributes;
  };

  //////////////////////////////////////////////////////////////////////////////
  //
  // Molecule
  //
  //////////////////////////////////////////////////////////////////////////////

  template<>
  inline void clear_molecule<SmartMol>(SmartMol &mol)
  {
    mol.clear();
  }

  /////////////////////////////
  //
  // Atoms
  //
  /////////////////////////////

  template<>
  inline Size num_atoms<SmartMol>(const SmartMol &mol)
  {
    return mol.numAtoms();
  }

  template<>
  inline iterator_pair<molecule_traits<SmartMol>::atom_iter>
  get_atoms<SmartMol>(const SmartMol &mol)
  {
    return mol.atoms();
  }


  template<>
  inline molecule_traits<SmartMol>::atom_type
  get_atom<SmartMol>(const SmartMol &mol, Index index)
  {
    PRE(index < mol.numAtoms());
    return mol.atom(index);
  }

  template<>
  inline molecule_traits<SmartMol>::atom_type add_atom<SmartMol>(SmartMol &mol)
  {
    return mol.addAtom();
  }

  template<>
  inline void remove_atom<SmartMol>(SmartMol &mol, molecule_traits<SmartMol>::atom_type atom)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    mol.removeAtom(atom);
  }

  /////////////////////////////
  //
  // Bonds
  //
  /////////////////////////////

  template<>
  inline Size num_bonds<SmartMol>(const SmartMol &mol)
  {
    return mol.numBonds();
  }

  template<>
  inline iterator_pair<molecule_traits<SmartMol>::bond_iter>
  get_bonds<SmartMol>(const SmartMol &mol)
  {
    return mol.bonds();
  }

  template<>
  inline molecule_traits<SmartMol>::bond_type
  get_bond<SmartMol>(const SmartMol &mol, Index index)
  {
    PRE(index < mol.numBonds());
    return mol.bond(index);
  }

  template<>
  inline molecule_traits<SmartMol>::bond_type add_bond<SmartMol>(SmartMol &mol,
      molecule_traits<SmartMol>::atom_type source,
      molecule_traits<SmartMol>::atom_type target)
  {
    PRE(source != molecule_traits<SmartMol>::null_atom());
    PRE(target != molecule_traits<SmartMol>::null_atom());
    return mol.addBond(source, target);
  }

  template<>
  inline void remove_bond<SmartMol>(SmartMol &mol, molecule_traits<SmartMol>::bond_type bond)
  {
    PRE(bond != molecule_traits<SmartMol>::null_bond());
    mol.removeBond(bond);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Atom
  //
  //////////////////////////////////////////////////////////////////////////////

  template<>
  inline Index get_index<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::atom_type atom)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    return atom.index();
  }

  template<>
  inline iterator_pair<molecule_traits<SmartMol>::incident_iter>
  get_bonds<SmartMol>(const SmartMol &mol, molecule_traits<SmartMol>::atom_type atom)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    return atom.bonds();
  }

  template<>
  inline iterator_pair<molecule_traits<SmartMol>::nbr_iter>
  get_nbrs<SmartMol>(const SmartMol &mol, molecule_traits<SmartMol>::atom_type atom)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    return atom.nbrs();
  }

  template<>
  inline bool is_aromatic<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::atom_type atom)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    return atom.isAromatic();
  }

  template<>
  inline void set_aromatic<SmartMol>(SmartMol &mol,
      molecule_traits<SmartMol>::atom_type atom, bool value)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    atom.setAromatic(value);
  }

  template<>
  inline int get_element<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::atom_type atom)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    return atom.element();
  }

  template<>
  inline void set_element<SmartMol>(SmartMol &mol,
      molecule_traits<SmartMol>::atom_type atom, int value)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    atom.setElement(value);
  }

  template<>
  inline int get_mass<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::atom_type atom)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    return atom.mass();
  }

  template<>
  inline void set_mass<SmartMol>(SmartMol &mol,
      molecule_traits<SmartMol>::atom_type atom, int value)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    atom.setMass(value);
  }

  template<>
  inline int get_degree<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::atom_type atom)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    return atom.degree();
  }

  template<>
  inline int get_hydrogens<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::atom_type atom)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    return atom.hydrogens();
  }

  template<>
  inline void set_hydrogens<SmartMol>(SmartMol &mol,
      molecule_traits<SmartMol>::atom_type atom, int value)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    atom.setHydrogens(value);
  }

  template<>
  inline int get_charge<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::atom_type atom)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    return atom.charge();
  }

  template<>
  inline void set_charge<SmartMol>(SmartMol &mol,
      molecule_traits<SmartMol>::atom_type atom, int value)
  {
    PRE(atom != molecule_traits<SmartMol>::null_atom());
    atom.setCharge(value);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // HeBond
  //
  //////////////////////////////////////////////////////////////////////////////

  template<>
  inline Index get_index<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::bond_type bond)
  {
    PRE(bond != molecule_traits<SmartMol>::null_bond());
    return bond.index();
  }

  template<>
  inline typename molecule_traits<SmartMol>::atom_type
  get_source<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::bond_type bond)
  {
    PRE(bond != molecule_traits<SmartMol>::null_bond());
    return bond.source();
  }

  template<>
  inline typename molecule_traits<SmartMol>::atom_type
  get_target<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::bond_type bond)
  {
    PRE(bond != molecule_traits<SmartMol>::null_bond());
    return bond.target();
  }

  template<>
  inline typename molecule_traits<SmartMol>::atom_type
  get_other<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::bond_type bond,
      const molecule_traits<SmartMol>::atom_type atom)
  {
    PRE(bond != molecule_traits<SmartMol>::null_bond());
    PRE(atom == bond.source() || atom == bond.target());
    return bond.other(atom);
  }

  template<>
  inline bool is_aromatic<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::bond_type bond)
  {
    PRE(bond != molecule_traits<SmartMol>::null_bond());
    return bond.isAromatic();
  }

  template<>
  inline void set_aromatic<SmartMol>(SmartMol &mol,
      molecule_traits<SmartMol>::bond_type bond, bool value)
  {
    PRE(bond != molecule_traits<SmartMol>::null_bond());
    bond.setAromatic(value);
  }

  template<>
  inline int get_order<SmartMol>(const SmartMol &mol,
      const molecule_traits<SmartMol>::bond_type bond)
  {
    PRE(bond != molecule_traits<SmartMol>::null_bond());
    return bond.order();
  }

  template<>
  inline void set_order<SmartMol>(SmartMol &mol,
      molecule_traits<SmartMol>::bond_type bond, int value)
  {
    PRE(bond != molecule_traits<SmartMol>::null_bond());
    bond.setOrder(value);
  }

  template<>
  inline molecule_traits<SmartMol>::bond_type
  get_bond<SmartMol>(const SmartMol &mol,
      molecule_traits<SmartMol>::atom_type source,
      molecule_traits<SmartMol>::atom_type target)
  {
    PRE(source != molecule_traits<SmartMol>::null_atom());
    PRE(target != molecule_traits<SmartMol>::null_atom());
    for (auto &bond : get_bonds(mol, source))
      if (get_other(mol, bond, source) == target)
        return bond;
    return mol.null_bond();
  }

  SmartMol::bond_type SmartMol::bond(const SmartMol::atom_type &source,
      const SmartMol::atom_type &target) const
  {
    bond_type bond = get_bond(*this, source, target);
    if (bond == molecule_traits<SmartMol>::null_bond())
      throw std::runtime_error(make_string("There is no bond between atoms ",
            source.index(), " and ", target.index()));
    return bond;
  }

  /**
   * @brief STL output stream operator for SmartMol.
   */
  inline std::ostream& operator<<(std::ostream &os, SmartMol &mol)
  {
    os << "Molecule:" << std::endl;

    os << "    Atoms:\tindex\telement" << std::endl;
    for (auto &atom : get_atoms(mol))
      os << "          \t" << atom.index() << "\t" << atom.element() << std::endl;

    os << "    Bonds:\tsource\ttarget\torder" << std::endl;
    for (auto &bond : get_bonds(mol))
      os << "          \t" << bond.source().index() << "\t" << bond.target().index() << "\t" << bond.order() << std::endl;

    return os;
  }

  inline std::ostream& operator<<(std::ostream &os, const SmartMol::atom_type &atom)
  {
    os << "SmartAtom(" << atom.index() << ")";
    return os;
  }

  inline std::ostream& operator<<(std::ostream &os, const SmartMol::bond_type &bond)
  {
    os << "SmartBond(" << bond.index() << ")";
    return os;
  }

  //@endcond

}

#endif
