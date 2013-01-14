#ifndef HELIUM_MOLECULE_H
#define HELIUM_MOLECULE_H

#include <utility>

namespace Helium {

  typedef unsigned int Index;
  typedef unsigned int Size;

  template<typename MoleculeType>
  struct molecule_traits
  {
    typedef typename MoleculeType::atom_type atom_type;
    typedef typename MoleculeType::bond_type bond_type;

    typedef typename MoleculeType::mol_atom_iter mol_atom_iter;
    typedef typename MoleculeType::mol_bond_iter mol_bond_iter;
    typedef typename MoleculeType::atom_atom_iter atom_atom_iter;
    typedef typename MoleculeType::atom_bond_iter atom_bond_iter;

    static Index null_index()
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
  Size num_atoms(const MoleculeType *mol);

  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::mol_atom_iter, typename molecule_traits<MoleculeType>::mol_atom_iter>
  get_atoms(MoleculeType *mol);

  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::atom_type get_atom(const MoleculeType *mol, Index index);

  /////////////////////////////
  //
  // Bonds
  //
  /////////////////////////////

  template<typename MoleculeType>
  Size num_bonds(const MoleculeType *mol);

  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::mol_bond_iter, typename molecule_traits<MoleculeType>::mol_bond_iter>
  get_bonds(MoleculeType *mol);

  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::bond_type get_bond(const MoleculeType *mol, Index index);

  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::bond_type get_bond(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type source,
                                                                                      typename molecule_traits<MoleculeType>::atom_type target);

  //////////////////////////////////////////////////////////////////////////////
  //
  // Atom
  //
  //////////////////////////////////////////////////////////////////////////////

  template<typename MoleculeType>
  Index get_index(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom);

  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::atom_bond_iter, typename molecule_traits<MoleculeType>::atom_bond_iter>
  get_bonds(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom);

  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::atom_atom_iter, typename molecule_traits<MoleculeType>::atom_atom_iter>
  get_nbrs(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom);

  template<typename MoleculeType>
  bool is_aromatic(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom);

  template<typename MoleculeType>
  bool is_cyclic(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom);

  template<typename MoleculeType>
  int get_element(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom);

  template<typename MoleculeType>
  int get_mass(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom);

  template<typename MoleculeType>
  int get_degree(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom);

  template<typename MoleculeType>
  int num_hydrogens(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom);

  template<typename MoleculeType>
  int get_charge(const MoleculeType *mol, typename molecule_traits<MoleculeType>::atom_type atom);

  //////////////////////////////////////////////////////////////////////////////
  //
  // Bond
  //
  //////////////////////////////////////////////////////////////////////////////

  template<typename MoleculeType>
  Index get_index(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond);

  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::atom_type get_source(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond);

  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::atom_type get_target(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond);

  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::atom_type get_other(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond,
                                                                                       typename molecule_traits<MoleculeType>::atom_type atom);

  template<typename MoleculeType>
  bool is_aromatic(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond);

  template<typename MoleculeType>
  bool is_cyclic(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond);

  template<typename MoleculeType>
  bool get_order(const MoleculeType *mol, typename molecule_traits<MoleculeType>::bond_type bond);

}

#endif
