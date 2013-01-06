#ifndef HELIUM_MOLECULE_H
#define HELIUM_MOLECULE_H

#include <vector>
#include <istream>
#include <algorithm>

namespace Helium {

  class Bond;
  
  template<typename AtomType, typename BondType>
  class nbr_iterator
  {
    public:
      nbr_iterator()
      {
      }

      nbr_iterator(AtomType *atom, typename std::vector<BondType*>::iterator iter) : m_atom(atom), m_iter(iter)
      {
      }

      AtomType* operator*() const
      {
        return (*m_iter)->other(m_atom);
      }

      nbr_iterator<AtomType, BondType>& operator++()
      {
        ++m_iter;
        return *this;
      }

      nbr_iterator<AtomType, BondType> operator++(int)
      {
        nbr_iterator<AtomType, BondType> tmp = *this;
        ++m_iter;
        return tmp;
      }

      bool operator!=(const nbr_iterator<AtomType, BondType> &other)
      {
        return m_iter != other.m_iter;
      }

    private:
      AtomType *m_atom;
      typename std::vector<BondType*>::iterator m_iter;
  };

  class Atom
  {
    public:
      Atom(int index, bool aromatic, bool cyclic, int element, int mass,
           int hydrogens, int charge) : m_index(index), m_aromatic(aromatic),
           m_cyclic(cyclic), m_element(element), m_mass(mass),
           m_hydrogens(hydrogens), m_charge(charge)
      {
      }

      std::pair<std::vector<Bond*>::iterator, std::vector<Bond*>::iterator> bonds()
      {
        return std::make_pair(m_bonds.begin(), m_bonds.end());
      }

      std::pair<nbr_iterator<Atom, Bond>, nbr_iterator<Atom, Bond> > nbrs()
      {
        return std::make_pair(nbr_iterator<Atom, Bond>(this, m_bonds.begin()),
                              nbr_iterator<Atom, Bond>(this, m_bonds.end()));
      }

      unsigned int index() const
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

      void addBond(Bond *bond)
      {
        m_bonds.push_back(bond);
      }

    private:
      std::vector<Bond*> m_bonds;
      unsigned int m_index;
      bool m_aromatic;
      bool m_cyclic;
      unsigned char m_element;
      unsigned char m_mass;
      unsigned char m_hydrogens;
      signed char m_charge;
  };

  class Bond
  {
    public:
      Bond(unsigned int index, Atom *source, Atom *target, bool aromatic, bool cyclic,
           int order) : m_source(source), m_target(target), m_index(index),
           m_aromatic(aromatic), m_cyclic(cyclic), m_order(order)
      {
      }

      unsigned int index() const
      {
        return m_index;
      }

      Atom* source() const
      {
        return m_source;
      }

      Atom* target() const
      {
        return m_target;
      }

      Atom* other(const Atom *atom) const
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
      Atom *m_source;
      Atom *m_target;
      unsigned int m_index;
      bool m_aromatic;
      bool m_cyclic;
      unsigned char m_order;
  };

  template<typename AtomFactory, typename BondFactory>
  class HeliumMolecule;

  template<typename MoleculeType>
  bool read_molecule(std::istream &is, MoleculeType &mol);

  template<typename AtomFactory, typename BondFactory>
  class HeliumMolecule
  {
    public:
      typedef Atom* atom_type;
      typedef Bond* bond_type;

      typedef std::vector<Atom*>::iterator mol_atom_iter;
      typedef std::vector<Bond*>::iterator mol_bond_iter;
      typedef std::vector<Bond*>::iterator atom_bond_iter;
      typedef nbr_iterator<Atom, Bond> atom_atom_iter;

      typedef AtomFactory atom_factory;
      typedef BondFactory bond_factory;

      ~HeliumMolecule()
      {
        for (std::size_t i = 0; i < m_atoms.size(); ++i)
          atom_factory::instance()->destroy(m_atoms[i]);
        for (std::size_t i = 0; i < m_bonds.size(); ++i)
          bond_factory::instance()->destroy(m_bonds[i]);
      }

      std::size_t numAtoms() const
      {
        return m_atoms.size();
      }

      std::size_t numBonds() const
      {
        return m_bonds.size();
      }

      std::pair<std::vector<Atom*>::iterator, std::vector<Atom*>::iterator> atoms()
      {
        return std::make_pair(m_atoms.begin(), m_atoms.end());
      }

      std::pair<std::vector<Bond*>::iterator, std::vector<Bond*>::iterator> bonds()
      {
        return std::make_pair(m_bonds.begin(), m_bonds.end());
      }

      Atom* atom(unsigned int index) const
      {
        return m_atoms[index];
      }

      Bond* bond(unsigned int index) const
      {
        return m_bonds[index];
      }

      static unsigned int null_index()
      {
        return -1;
      }

      static Atom* null_atom()
      {
        return 0;
      }

      static Bond* null_bond()
      {
        return 0;
      }

     private:
      template<typename MoleculeType> friend bool read_molecule(std::istream &is, MoleculeType &mol);

      std::vector<Atom*> m_atoms;
      std::vector<Bond*> m_bonds;
  };

  /*
  */
  class AtomBondFactory
  {
    public:
      static AtomBondFactory* instance()
      {
        static AtomBondFactory *singleton = 0;
        if (!singleton)
          singleton = new AtomBondFactory;
        return singleton;
      }

      Atom* create(int index, bool aromatic, bool cyclic, int element, 
           int mass, int hydrogens, int charge)
      {
        return new Atom(index, aromatic, cyclic, element, mass, hydrogens, charge);
      }

      Bond* create(unsigned int index, Atom *source, Atom *target, bool aromatic, bool cyclic, int order)
      {
        return new Bond(index, source, target, aromatic, cyclic, order);
      }

      void destroy(Atom *atom)
      {
        delete atom;
      }

      void destroy(Bond *bond)
      {
        delete bond;
      }
 
  };


  /*
  class AtomBondFactory
  {
    public:
      static AtomBondFactory* instance()
      {
        static AtomBondFactory *singleton = 0;
        if (!singleton)
          singleton = new AtomBondFactory;
        return singleton;
      }

      Atom* create(int index, bool aromatic, bool cyclic, int element, 
           int mass, int hydrogens, int charge)
      {
        Atom *atom;
        if (!m_destroyedAtoms.empty()) {
          atom = new (&m_atoms[m_destroyedAtoms.back()]) Atom(index, aromatic, cyclic, element,
                                                              mass, hydrogens, charge);
          m_destroyedAtoms.pop_back();
        } else {
          m_atoms.push_back(Atom(index, aromatic, cyclic, element, mass,
                                 hydrogens, charge));
          atom = &m_atoms.back();
        }
        return atom;
      }

      Bond* create(unsigned int index, Atom *source, Atom *target, bool aromatic, bool cyclic, int order)
      {
        Bond *bond;
        if (!m_destroyedBonds.empty()) {
          bond = new (&m_bonds[m_destroyedBonds.back()]) Bond(index, source, target, aromatic, cyclic, order);
          m_destroyedBonds.pop_back();
        } else {
          m_bonds.push_back(Bond(index, source, target, aromatic, cyclic, order));
          bond = &m_bonds.back();
        }
        return bond;
      }

      void destroy(Atom *atom)
      {
        std::size_t index = atom - &m_atoms.front();
        m_destroyedAtoms.push_back(index);
      }

      void destroy(Bond *bond)
      {
        std::size_t index = bond - &m_bonds.front();
        m_destroyedBonds.push_back(index);
      }
 
    private:
      std::vector<Atom> m_atoms;
      std::vector<Bond> m_bonds;
      std::vector<unsigned int> m_destroyedAtoms;
      std::vector<unsigned int> m_destroyedBonds;
  };

  std::vector<Atom> AtomBondFactory::m_atoms;
  std::vector<Bond> AtomBondFactory::m_bonds;
  std::vector<unsigned int> AtomBondFactory::m_destroyedAtoms;
  std::vector<unsigned int> AtomBondFactory::m_destroyedBonds;
  */

  typedef HeliumMolecule<AtomBondFactory, AtomBondFactory> Molecule;


  template<typename MoleculeType>
  struct molecule_traits
  {
    typedef typename MoleculeType::atom_type atom_type;
    typedef typename MoleculeType::bond_type bond_type;

    typedef typename MoleculeType::mol_atom_iter mol_atom_iter;
    typedef typename MoleculeType::mol_bond_iter mol_bond_iter;
    typedef typename MoleculeType::atom_atom_iter atom_atom_iter;
    typedef typename MoleculeType::atom_bond_iter atom_bond_iter;

    typedef typename MoleculeType::atom_factory atom_factory;
    typedef typename MoleculeType::bond_factory bond_factory;

    static unsigned int null_index()
    {
      return MoleculeType::null_index();
    }

    static atom_type null_atom()
    {
      return MoleculeType::null_atom();
    }

    static bond_type null_bond()
    {
      return MoleculeType::null_bond();
    }
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

  template<typename MoleculeType>
  std::size_t num_atoms(const MoleculeType *mol)
  {
    return mol->numAtoms();
  }

  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::mol_atom_iter, typename molecule_traits<MoleculeType>::mol_atom_iter>
  get_atoms(MoleculeType *mol)
  {
    return mol->atoms();
  }

  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::atom_type get_atom(const MoleculeType *mol, std::size_t index)
  {
    return mol->atom(index);
  }


  /////////////////////////////
  //
  // Bonds
  //
  /////////////////////////////

  template<typename MoleculeType>
  std::size_t num_bonds(const MoleculeType *mol)
  {
    return mol->numBonds();
  }

  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::mol_bond_iter, typename molecule_traits<MoleculeType>::mol_bond_iter>
  get_bonds(MoleculeType *mol)
  {
    return mol->bonds();
  }

  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::bond_type get_bond(const MoleculeType *mol, std::size_t index)
  {
    return mol->bond(index);
  }

  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::bond_type get_bond(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type source,
                                                                                      typename molecule_traits<MoleculeType>::atom_type target)
  {
    typename molecule_traits<MoleculeType>::atom_bond_iter bond, end_bonds;
    tie(bond, end_bonds) = get_bonds(mol, source);
    for (; bond != end_bonds; ++bond)
      if (get_other(mol, *bond, source) == target)
        return *bond;
    return 0;
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Atom
  //
  //////////////////////////////////////////////////////////////////////////////

  template<typename MoleculeType>
  std::size_t get_index(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return atom->index();
  }

  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::atom_bond_iter, typename molecule_traits<MoleculeType>::atom_bond_iter>
  get_bonds(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return atom->bonds();
  }

  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::atom_atom_iter, typename molecule_traits<MoleculeType>::atom_atom_iter>
  get_nbrs(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return atom->nbrs();
  }

  template<typename MoleculeType>
  bool is_aromatic(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return atom->isAromatic();
  }

  template<typename MoleculeType>
  bool is_cyclic(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return atom->isCyclic();
  }

  template<typename MoleculeType>
  int get_element(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return atom->element();
  }

  template<typename MoleculeType>
  int get_mass(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return atom->mass();
  }

  template<typename MoleculeType>
  int get_degree(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return atom->degree();
  }

  template<typename MoleculeType>
  int num_hydrogens(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return atom->hydrogens();
  }

  template<typename MoleculeType>
  int get_charge(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return atom->charge();
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Bond
  //
  //////////////////////////////////////////////////////////////////////////////

  template<typename MoleculeType>
  std::size_t get_index(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond)
  {
    return bond->index();
  }

  template<typename MoleculeType>
  Atom* get_source(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond)
  {
    return bond->source();
  }

  template<typename MoleculeType>
  Atom* get_target(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond)
  {
    return bond->target();
  }

  template<typename MoleculeType>
  Atom* get_other(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return bond->other(atom);
  }

  template<typename MoleculeType>
  bool is_aromatic(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond)
  {
    return bond->isAromatic();
  }

  template<typename MoleculeType>
  bool is_cyclic(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond)
  {
    return bond->isCyclic();
  }

  template<typename MoleculeType>
  bool get_order(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond)
  {
    return bond->order();
  }


}

#endif
