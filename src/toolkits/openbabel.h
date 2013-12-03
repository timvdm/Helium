#ifndef HELIUM_OPENBAEL_H
#define HELIUM_OPENABEL_H

#include <Helium/molecule.h>
#include <openbabel/mol.h>

namespace Helium {

  namespace impl {

    class openbabel_nbr_iterator
    {
      public:
        openbabel_nbr_iterator()
        {
        }

        openbabel_nbr_iterator(OpenBabel::OBAtom *atom) : m_iter(atom)
        {
        }

        OpenBabel::OBAtom* operator*() const
        {
          return m_iter.operator->();
        }

        openbabel_nbr_iterator& operator++()
        {
          ++m_iter;
          return *this;
        }

        openbabel_nbr_iterator operator++(int)
        {
          openbabel_nbr_iterator tmp = *this;
          ++m_iter;
          return tmp;
        }

        bool operator!=(const openbabel_nbr_iterator &other)
        {
          return m_iter != other.m_iter;
        }

      private:
        OpenBabel::OBAtomAtomIter m_iter;
    };

  }

  template<>
  struct molecule_traits<OpenBabel::OBMol>
  {
    typedef OpenBabel::OBAtom* atom_type;
    typedef OpenBabel::OBBond* bond_type;

    typedef OpenBabel::OBAtomIterator atom_iter;
    typedef OpenBabel::OBBondIterator bond_iter;
    typedef impl::openbabel_nbr_iterator nbr_iter;
    typedef OpenBabel::OBBondIterator incident_iter;

    static Index null_index()
    {
      return -1;
    }

    static atom_type null_atom()
    {
      return 0;
    }

    static bond_type null_bond()
    {
      return 0;
    }
  };

  template<>
  void clear_molecule<OpenBabel::OBMol>(OpenBabel::OBMol &mol)
  {
    mol.Clear();
  }

  template<>
  Size num_atoms(const OpenBabel::OBMol &mol)
  {
    return mol.NumAtoms();
  }

  template<>
  std::pair<typename molecule_traits<OpenBabel::OBMol>::atom_iter,
            typename molecule_traits<OpenBabel::OBMol>::atom_iter>
  get_atoms(const OpenBabel::OBMol &mol)
  {
    return std::make_pair(const_cast<OpenBabel::OBMol&>(mol).BeginAtoms(),
                          const_cast<OpenBabel::OBMol&>(mol).EndAtoms());
  }

  template<>
  typename molecule_traits<OpenBabel::OBMol>::atom_type
  get_atom(const OpenBabel::OBMol &mol, Index index)
  {
    return mol.GetAtom(index + 1);
  }

  template<>
  typename molecule_traits<OpenBabel::OBMol>::atom_type
  add_atom(OpenBabel::OBMol &mol)
  {
    return mol.NewAtom();
  }

  template<>
  void remove_atom(OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom)
  {
    mol.DeleteAtom(atom);
  }

  template<>
  Size num_bonds(const OpenBabel::OBMol &mol)
  {
    return mol.NumBonds();
  }

  template<>
  std::pair<typename molecule_traits<OpenBabel::OBMol>::bond_iter,
            typename molecule_traits<OpenBabel::OBMol>::bond_iter>
  get_bonds(const OpenBabel::OBMol &mol)
  {
    return std::make_pair(const_cast<OpenBabel::OBMol&>(mol).BeginBonds(),
                          const_cast<OpenBabel::OBMol&>(mol).EndBonds());
  }

  template<>
  typename molecule_traits<OpenBabel::OBMol>::bond_type
  get_bond(const OpenBabel::OBMol &mol, Index index)
  {
    return mol.GetBond(index);
  }

  template<>
  typename molecule_traits<OpenBabel::OBMol>::bond_type
  get_bond(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type source,
      typename molecule_traits<OpenBabel::OBMol>::atom_type target)
  {
    return mol.GetBond(source, target);
  }

  template<>
  typename molecule_traits<OpenBabel::OBMol>::bond_type
  add_bond(OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type source,
      typename molecule_traits<OpenBabel::OBMol>::atom_type target)
  {
    mol.AddBond(source->GetIdx(), target->GetIdx(), 1);
    return mol.GetBond(mol.NumBonds() - 1);
  }

  template<>
  void remove_bond(OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::bond_type bond)
  {
    mol.DeleteBond(bond);
  }

  template<>
  Index get_index(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom)
  {
    return atom->GetIndex();
  }

  template<>
  std::pair<typename molecule_traits<OpenBabel::OBMol>::incident_iter,
            typename molecule_traits<OpenBabel::OBMol>::incident_iter>
  get_bonds(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom)
  {
    return std::make_pair(atom->BeginBonds(), atom->EndBonds());
  }

  template<>
  std::pair<typename molecule_traits<OpenBabel::OBMol>::nbr_iter,
    typename molecule_traits<OpenBabel::OBMol>::nbr_iter>
  get_nbrs(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom)
  {
    return std::make_pair(impl::openbabel_nbr_iterator(atom),
                          impl::openbabel_nbr_iterator());
  }

  template<>
  bool is_aromatic(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom)
  {
    return atom->IsAromatic();
  }

  template<>
  void set_aromatic(OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom,
      bool value)
  {
    if (value)
      atom->SetAromatic();
    else
      atom->UnsetAromatic();
  }

  template<>
  int get_element(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom)
  {
    return atom->GetAtomicNum();
  }

  template<>
  void set_element(OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom,
      int value)
  {
    atom->SetAtomicNum(value);
    atom->SetIsotope(Element::averageMass(value));
  }

  template<>
  int get_mass(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom)
  {
    return atom->GetIsotope();
  }

  template<>
  void set_mass(OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom,
      int value)
  {
    atom->SetIsotope(value);
  }

  template<>
  int get_degree(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom)
  {
    return atom->GetValence();
  }

  template<>
  int num_hydrogens(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom)
  {
    return atom->GetImplicitValence() - atom->GetValence();
  }

  template<>
  void set_hydrogens(OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom,
      int value)
  {
    atom->SetImplicitValence(value - atom->GetValence());
  }

  template<>
  int get_charge(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom)
  {
    return atom->GetFormalCharge();
  }

  template<>
  void set_charge(OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom,
      int value)
  {
    atom->SetFormalCharge(value);
  }

  template<>
  Index get_index(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::bond_type bond)
  {
    return bond->GetIdx();
  }

  template<>
  typename molecule_traits<OpenBabel::OBMol>::atom_type
  get_source(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::bond_type bond)
  {
    return bond->GetBeginAtom();
  }

  template<>
  typename molecule_traits<OpenBabel::OBMol>::atom_type
  get_target(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::bond_type bond)
  {
    return bond->GetEndAtom();
  }

  template<>
  typename molecule_traits<OpenBabel::OBMol>::atom_type
  get_other(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::bond_type bond,
      typename molecule_traits<OpenBabel::OBMol>::atom_type atom)
  {
    return bond->GetNbrAtom(atom);
  }

  template<>
  bool is_aromatic(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::bond_type bond)
  {
    return bond->IsAromatic();
  }

  template<>
  void set_aromatic(OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::bond_type bond,
      bool value)
  {
    if (value)
      bond->SetAromatic();
    else
      bond->UnsetAromatic();
  }

  template<>
  int get_order(const OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::bond_type bond)
  {
    return bond->GetBondOrder();
  }

  template<>
  void set_order(OpenBabel::OBMol &mol,
      typename molecule_traits<OpenBabel::OBMol>::bond_type bond,
      int value)
  {
    bond->SetBondOrder(value);
  }

}

#endif
