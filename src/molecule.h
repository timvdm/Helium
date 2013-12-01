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
#ifndef HELIUM_MOLECULE_H
#define HELIUM_MOLECULE_H

#include <utility>
#include <vector>
#include <sstream>

/**
 * @brief Iterate over all the atoms in a molecule.
 *
 * @param atom The atom.
 * @param mol The molecule.
 * @param MoleculeType The type of the molecule.
 */
#define FOREACH_ATOM(atom, mol, MoleculeType) \
  for (typename Helium::impl::ForeachAtom<MoleculeType> atom(mol); atom.begin != atom.end; ++atom.begin)

/**
 * @brief Iterate over all the bonds in a molecule.
 *
 * @param bond The bond.
 * @param mol The molecule.
 * @param MoleculeType The type of the molecule.
 */
#define FOREACH_BOND(bond, mol, MoleculeType) \
  for (typename Helium::impl::ForeachBond<MoleculeType> bond(mol); bond.begin != bond.end; ++bond.begin)

/**
 * @brief Iterate over all the neighbors of an atom.
 *
 * @param nbr The neighbor atom.
 * @param atom The center atom.
 * @param mol The molecule.
 * @param MoleculeType The type of the molecule.
 */
#define FOREACH_NBR(nbr, atom, mol, MoleculeType) \
  for (typename Helium::impl::ForeachNbr<MoleculeType> nbr(mol, atom); nbr.begin != nbr.end; ++nbr.begin)

/**
 * @brief Iterate over all the incident bonds of an atom.
 *
 * @param bond The incident bond.
 * @param atom The center atom.
 * @param mol The molecule.
 * @param MoleculeType The type of the molecule.
 */
#define FOREACH_INCIDENT(bond, atom, mol, MoleculeType) \
  for (typename Helium::impl::ForeachIncident<MoleculeType> bond(mol, atom); bond.begin != bond.end; ++bond.begin)

namespace Helium {

  /**
   * @defgroup molecule_group Molecule Concept
   * @{
   */

  /**
   * @file molecule.h
   * @brief Type traits and functions for the Molecule concept.
   */

  /**
   * @brief Type used for atom and bond indices.
   */
  typedef unsigned int Index;
  /**
   * @brief Type used for various sizes.
   *
   * Examples include number of atoms/bonds, number of neighbors, ...
   */
  typedef unsigned int Size;

  /**
   * @brief Molecule type traits struct.
   *
   * The molecule_traits struct template is used to associate various types
   * with a model of the Molecule concept. When writing a new class to model
   * the Molecule concept, the default typedefs and function names can be
   * added to the new class and in these cases no molecule_traits specialization
   * has to be provided. In other cases where an existing class is to be made
   * a model of the Molecule concept, a specialization of molecule_traits
   * has to be provided for the specific type.
   */
  template<typename MoleculeType>
  struct molecule_traits
  {
    /**
     * @brief The type of an atom.
     *
     * By default, the atom_type typedef inside MoleculeType is used.
     */
    typedef typename MoleculeType::atom_type atom_type;
    /**
     * @brief The type of a bond.
     *
     * By default, the bond_type typedef inside MoleculeType is used.
     */
    typedef typename MoleculeType::bond_type bond_type;

    /**
     * @brief The type of an iterator over a molecule's atoms.
     *
     * By default, the atom_iter typedef inside MoleculeType is used.
     */
    typedef typename MoleculeType::atom_iter atom_iter;
    /**
   * @brief The type of an iterator over a molecule's bonds.
     *
     * By default, the bond_iter typedef inside MoleculeType is used.
     */
    typedef typename MoleculeType::bond_iter bond_iter;
    /**
     * @brief The type of an iterator over an atom's neighboring atoms.
     *
     * By default, the nbr_iter typedef inside MoleculeType is used.
     */
    typedef typename MoleculeType::nbr_iter nbr_iter;
    /**
     * @brief The type of an iterator over an atom's incident bonds.
     *
     * By default, the incident_iter typedef inside MoleculeType is used.
     */
    typedef typename MoleculeType::incident_iter incident_iter;

    /**
     * @brief Get the null value for an atom/bond index.
     *
     * This index value indicates an invalid index and has no atom/bond
     * associated with it. By default, the static null_index() function
     * inside MoleculeType is used.
     */
    static Index null_index()
    {
      return MoleculeType::null_index();
    }

    /**
     * @brief Get a null atom (i.e. an invalid atom)
     *
     * A null atom is an atom that is not valid. This could be the result after
     * creating an atom using the default constructor or simply a 0 if the
     * atom_type is MyAtom* for example. By default, the static null_atom()
     * function inside MoleculeType is used.
     */
    static atom_type null_atom()
    {
      return MoleculeType::null_atom();
    }

    /**
     * @brief Get a null bond (i.e. an invalid bond)
     *
     * A null bond is a bond that is not valid. This could be the result after
     * creating a bond using the default constructor or simply a 0 if the
     * bond_type is MyBond* for example. By default, the static null_bond()
     * function inside MoleculeType is used.
     */
    static bond_type null_bond()
    {
      return MoleculeType::null_bond();
    }
  };

  namespace impl {

    template<typename MoleculeType>
    struct ForeachAtom
    {
      ForeachAtom(const MoleculeType &mol)
      {
        TIE(begin, end) = get_atoms(mol);
      }

      typename molecule_traits<MoleculeType>::atom_type operator*()
      {
        return *begin;
      }

      typename molecule_traits<MoleculeType>::atom_iter begin, end;
    };

    template<typename MoleculeType>
    struct ForeachBond
    {
      ForeachBond(const MoleculeType &mol)
      {
        TIE(begin, end) = get_bonds(mol);
      }

      typename molecule_traits<MoleculeType>::bond_type operator*()
      {
        return *begin;
      }

      typename molecule_traits<MoleculeType>::bond_iter begin, end;
    };

    template<typename MoleculeType>
    struct ForeachNbr
    {
      ForeachNbr(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
      {
        TIE(begin, end) = get_nbrs(mol, atom);
      }

      typename molecule_traits<MoleculeType>::atom_type operator*()
      {
        return *begin;
      }

      typename molecule_traits<MoleculeType>::nbr_iter begin, end;
    };

    template<typename MoleculeType>
    struct ForeachIncident
    {
      ForeachIncident(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
      {
        TIE(begin, end) = get_bonds(mol, atom);
      }

      typename molecule_traits<MoleculeType>::bond_type operator*()
      {
        return *begin;
      }

      typename molecule_traits<MoleculeType>::incident_iter begin, end;
    };

  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Molecule
  //
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Clear the molecule (i.e. remove all atoms/bonds).
   *
   * @param mol The molecule.
   */
  template<typename EditableMoleculeType>
  void clear_molecule(EditableMoleculeType &mol);

  /////////////////////////////
  //
  // Atoms
  //
  /////////////////////////////

  /**
   * @brief Get the number of atoms in a molecule.
   *
   * @param mol The molecule.
   *
   * @return The number of atoms.
   */
  template<typename MoleculeType>
  Size num_atoms(const MoleculeType &mol);

  /**
   * @brief Get an iterator pair to iterator over the atoms inside a molecule.
   *
   * @param mol The molecule.
   *
   * @return The iterator pair to iterate over the atoms inside a molecule.
   */
  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::atom_iter,
    typename molecule_traits<MoleculeType>::atom_iter>
  get_atoms(const MoleculeType &mol);

  /**
   * @brief Get the atom with the specified index.
   *
   * @pre The index must be valid (i.e. in the range [0,n) where n is the number
   *      of atoms)
   *
   * @param mol The molecule.
   * @param index The index of the atom to get.
   *
   * @return The atom with the specified @p index.
   */
  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::atom_type
  get_atom(const MoleculeType &mol, Index index);

  /**
   * @brief Add an atom to the molecule.
   *
   * @param mol The molecule.
   *
   * @return The newly created atom.
   */
  template<typename EditableMoleculeType>
  typename molecule_traits<EditableMoleculeType>::atom_type
  add_atom(EditableMoleculeType &mol);

  /**
   * @brief Remove an atom from the molecule.
   *
   * @param mol The molecule.
   * @param atom The atom to remove.
   */
  template<typename EditableMoleculeType>
  void remove_atom(EditableMoleculeType &mol,
      typename molecule_traits<EditableMoleculeType>::atom_type atom);

  /////////////////////////////
  //
  // Bonds
  //
  /////////////////////////////

  /**
   * @brief Get the number of bonds in a molecule.
   *
   * @param mol The molecule.
   *
   * @return The number of bonds.
   */
  template<typename MoleculeType>
  Size num_bonds(const MoleculeType &mol);

  /**
   * @brief Get an iterator pair to iterator over the bonds inside a molecule.
   *
   * @param mol The molecule.
   *
   * @return The iterator pair to iterate over the bonds inside a molecule.
   */
  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::bond_iter,
            typename molecule_traits<MoleculeType>::bond_iter>
  get_bonds(const MoleculeType &mol);

  /**
   * @brief Get the bond with the specified index.
   *
   * @pre The index must be valid (i.e. in the range [0,n) where n is the number
   *      of bonds)
   *
   * @param mol The molecule.
   * @param index The index of the bond to get.
   *
   * @return The bond with the specified @p index.
   */
  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::bond_type
  get_bond(const MoleculeType &mol, Index index);

  /**
   * @brief Get the bond between the specified source and target atoms.
   *
   * @pre The source and target atoms must be valid (i.e. not equal to
   *       molecule_traits<MoleculeType>::null_atom())
   *
   * @param mol The molecule.
   * @param source One atom of the bond.
   * @param target Another atom of the bond.
   *
   * @return The bond between the specified atom or
   *         molecule_traits<MoleculeType>::null_bond() if there is no such bond.
   */
  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::bond_type
  get_bond(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type source,
      typename molecule_traits<MoleculeType>::atom_type target);

  /**
   * @brief Add a bond to the molecule.
   *
   * @param mol The molecule.
   *
   * @return The newly created bond.
   */
  template<typename EditableMoleculeType>
  typename molecule_traits<EditableMoleculeType>::bond_type
  add_bond(EditableMoleculeType &mol,
      typename molecule_traits<EditableMoleculeType>::atom_type source,
      typename molecule_traits<EditableMoleculeType>::atom_type target);

  /**
   * @brief Remove a bond from the molecule.
   *
   * @param mol The molecule.
   * @param bond The bond to remove.
   */
  template<typename EditableMoleculeType>
  void remove_bond(EditableMoleculeType &mol,
      typename molecule_traits<EditableMoleculeType>::bond_type bond);

  //////////////////////////////////////////////////////////////////////////////
  //
  // Atom
  //
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Get the atom's index.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   *
   * @return The atom index.
   */
  template<typename MoleculeType>
  Index get_index(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom);

  /**
   * @brief Get the iterator pair to iterate over an atom's incident bonds.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   *
   * @return The iterator pair to iterate over an atom's incident bonds.
   */
  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::incident_iter,
            typename molecule_traits<MoleculeType>::incident_iter>
  get_bonds(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom);

  /**
   * @brief Get the iterator pair to iterate over an atom's neighboring atoms.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   *
   * @return The iterator pair to iterate over an atom's neighboring atoms.
   */
  template<typename MoleculeType>
  std::pair<typename molecule_traits<MoleculeType>::nbr_iter,
    typename molecule_traits<MoleculeType>::nbr_iter>
  get_nbrs(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom);

  /**
   * @brief Get the atom's aromaticity.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   *
   * @return True if the atom is aromatic, false otherwise.
   */
  template<typename MoleculeType>
  bool is_aromatic(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom);

  /**
   * @brief Set the atom's aromaticity.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   * @param value The new value.
   */
  template<typename EditableMoleculeType>
  void set_aromatic(EditableMoleculeType &mol,
      typename molecule_traits<EditableMoleculeType>::atom_type atom,
      bool value);

  /**
   * @brief Get the atom's chemical element number.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   *
   * @return The atom's chemical element number.
   */
  template<typename MoleculeType>
  int get_element(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom);

  /**
   * @brief Set the atom's element.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   * @param value The new value.
   */
  template<typename EditableMoleculeType>
  void set_element(EditableMoleculeType &mol,
      typename molecule_traits<EditableMoleculeType>::atom_type atom,
      int value);

  /**
   * @brief Get the atom's mass number.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   *
   * @return The atom's mass number.
   */
  template<typename MoleculeType>
  int get_mass(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom);

  /**
   * @brief Set the atom's mass.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   * @param value The new value.
   */
  template<typename EditableMoleculeType>
  void set_mass(EditableMoleculeType &mol,
      typename molecule_traits<EditableMoleculeType>::atom_type atom,
      int value);

  /**
   * @brief Get the atom's degree.
   *
   * The degree is the number of neighboring atoms (or incident bonds).
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   *
   * @return The atom's degree.
   */
  template<typename MoleculeType>
  int get_degree(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom);

  /**
   * @brief Get the atom's number of connected hydrogens.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   *
   * @return The atom's number of connected hydrogens.
   */
  template<typename MoleculeType>
  int num_hydrogens(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom);

  /**
   * @brief Set the atom's number of hydrogens.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   * @param value The new value.
   */
  template<typename EditableMoleculeType>
  void set_hydrogens(EditableMoleculeType &mol,
      typename molecule_traits<EditableMoleculeType>::atom_type atom,
      int value);

  /**
   * @brief Get the atom's charge.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   *
   * @return The atom's charge.
   */
  template<typename MoleculeType>
  int get_charge(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom);

  /**
   * @brief Set the atom's charge.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   * @param value The new value.
   */
  template<typename EditableMoleculeType>
  void set_charge(EditableMoleculeType &mol,
      typename molecule_traits<EditableMoleculeType>::atom_type atom,
      int value);

  //////////////////////////////////////////////////////////////////////////////
  //
  // Bond
  //
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Get the bond's index.
   *
   * @pre The bond must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_bond()).
   *
   * @param mol The molecule.
   * @param bond The bond.
   *
   * @return The bond index.
   */
  template<typename MoleculeType>
  Index get_index(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond);

  /**
   * @brief Get the bond's source atom.
   *
   * @pre The bond must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_bond()).
   *
   * @param mol The molecule.
   * @param bond The bond.
   *
   * @return The bond's source atom.
   */
  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::atom_type get_source(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond);

  /**
   * @brief Get the bond's target atom.
   *
   * @pre The bond must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_bond()).
   *
   * @param mol The molecule.
   * @param bond The bond.
   *
   * @return The bond's target atom.
   */
  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::atom_type get_target(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond);

  /**
   * @brief Get the bond's other atom.
   *
   * @pre The bond must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_bond()).
   *      The specified @p atom must be one of the bond's atoms, if this is not the case
   *      any of the returned value is undefined.
   *
   * @param mol The molecule.
   * @param bond The bond.
   * @param atom One of the bond atoms, the other bond atom is returned.
   *
   * @return The bond's other atom (i.e. The bond atom that is not equal to @p atom).
   */
  template<typename MoleculeType>
  typename molecule_traits<MoleculeType>::atom_type get_other(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond,
                                                                                       typename molecule_traits<MoleculeType>::atom_type atom);

  /**
   * @brief Get the bond's aromaticity.
   *
   * @pre The bond must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_bond()).
   *
   * @param mol The molecule.
   * @param bond The bond.
   *
   * @return True if the bond is aromatic, false otherwise.
   */
  template<typename MoleculeType>
  bool is_aromatic(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond);

  /**
   * @brief Set the bond's aromaticity.
   *
   * @pre The bond must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_bond()).
   *
   * @param mol The molecule.
   * @param bond The bond.
   * @param value The new value.
   */
  template<typename EditableMoleculeType>
  void set_aromatic(EditableMoleculeType &mol,
      typename molecule_traits<EditableMoleculeType>::bond_type bond,
      bool value);

  /**
   * @brief Get the bond's order.
   *
   * @pre The bond must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_bond()).
   *
   * @param mol The molecule.
   * @param bond The bond.
   *
   * @return The bond's order.
   */
  template<typename MoleculeType>
  int get_order(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond);

  /**
   * @brief Set the bond's order.
   *
   * @pre The bond must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_bond()).
   *
   * @param mol The molecule.
   * @param bond The bond.
   * @param value The new value.
   */
  template<typename EditableMoleculeType>
  void set_order(EditableMoleculeType &mol,
      typename molecule_traits<EditableMoleculeType>::bond_type bond,
      int value);

  //////////////////////////////////////////////////////////////////////////////
  //
  // Utility Functions
  //
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Check if an atom is a hydrogen atom.
   *
   * @param mol The molecule.
   * @param atom The atom to check.
   *
   * @return True if the atom is a hydrogen  atom.
   */
  template<typename MoleculeType>
  bool is_hydrogen(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return get_element(mol, atom) == 1;
  }

  /**
   * @brief Check if an atom is a carbon atom.
   *
   * @param mol The molecule.
   * @param atom The atom to check.
   *
   * @return True if the atom is a carbon atom.
   */
  template<typename MoleculeType>
  bool is_carbon(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return get_element(mol, atom) == 6;
  }

  /**
   * @brief Check if an atom is a nitrogen atom.
   *
   * @param mol The molecule.
   * @param atom The atom to check.
   *
   * @return True if the atom is a nitrogen atom.
   */
  template<typename MoleculeType>
  bool is_nitrogen(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return get_element(mol, atom) == 7;
  }

  /**
   * @brief Check if an atom is a oxygen atom.
   *
   * @param mol The molecule.
   * @param atom The atom to check.
   *
   * @return True if the atom is a oxygen atom.
   */
  template<typename MoleculeType>
  bool is_oxygen(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return get_element(mol, atom) == 8;
  }

  /**
   * @brief Check if an atom is a phosphorus atom.
   *
   * @param mol The molecule.
   * @param atom The atom to check.
   *
   * @return True if the atom is a phosphorus atom.
   */
  template<typename MoleculeType>
  bool is_phosphorus(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return get_element(mol, atom) == 15;
  }

  /**
   * @brief Check if an atom is a sulfur atom.
   *
   * @param mol The molecule.
   * @param atom The atom to check.
   *
   * @return True if the atom is a sulfur atom.
   */
  template<typename MoleculeType>
  bool is_sulfur(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return get_element(mol, atom) == 16;
  }

  /**
   * @brief Get the number of attached heavy atoms.
   *
   * All atoms except hydrogen are heavy atoms.
   *
   * @param mol The molecule.
   * @param atom The atom to check.
   *
   * @return The number of attached heavy atoms.
   */
  template<typename MoleculeType>
  int get_heavy_degree(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    int degree = 0;
    FOREACH_NBR (nbr, atom, mol, MoleculeType)
      if (get_element(mol, *nbr) > 1)
        ++degree;
    return degree;
  }

  /**
   * @brief Get the valence atom.
   *
   * The valence is the bond order sum of the explicit bonds + the number of
   * implicit/explicit hydrogens.
   *
   * @param mol The molecule.
   * @param atom The atom to check.
   *
   * @return The valence for the atom.
   */
  template<typename MoleculeType>
  int get_valence(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    double val = 0;
    FOREACH_INCIDENT (bond, atom, mol, MoleculeType) {
      int order = get_order(mol, *bond);
      if (order == 5)
        val += 1.5;
      else
        val += order;
    }
    return val;
  }

  /**
   * @brief Get the atom connectivity.
   *
   * The connectivity is the number of attached heavy atoms + the number of
   * implicit/explicit hydrogens.
   *
   * @param mol The molecule.
   * @param atom The atom to check.
   *
   * @return The connectivity for the atom.
   */
  template<typename MoleculeType>
  int get_connectivity(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
  {
    return get_valence(mol, atom) + num_hydrogens(mol, atom);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Atom Predicates
  //
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Base class for atom predicates.
   *
   * Atom predicates are used in combination with a number of functions
   * (e.g. molecule_has_atom(), atom_has_nbr(), bond_has_atom(), ...).
   * An atom predicate is a functor that returns true if an atom matches
   * the description of the predicate. There are a number of predicates
   * predefined (e.g. ElementPredicate, ChargePredicate, ...). Additional
   * predicates can be added by inheriting AtomPredicate and implementing
   * the pure virtual function operator.
   */
  template<typename MoleculeType>
  class AtomPredicate
  {
    public:
      /**
       * @brief The the atom type.
       */
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

      /**
       * @brief Destructor.
       */
      virtual ~AtomPredicate()
      {
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param atom The atom to test.
       *
       * @return True if the atom matches the predicate.
       */
      virtual bool operator()(const MoleculeType &mol, atom_type atom) const = 0;
  };

  /**
   * @brief Class for combining atom predicates using a locgical and.
   */
  template<typename MoleculeType>
  class AtomAndPredicates : public AtomPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The atom type.
       */
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

      /**
       * @brief Constructor.
       *
       * @param predicate1 The first atom predicate.
       * @param predicate2 The second atom predicate.
       */
      template<typename Predicate1, typename Predicate2>
      AtomAndPredicates(const Predicate1 &predicate1, const Predicate2 &predicate2)
      {
        m_predicates.push_back(new Predicate1(predicate1));
        m_predicates.push_back(new Predicate2(predicate2));
      }

      /**
       * @brief Constructor.
       *
       * @param predicate1 The first atom predicate.
       * @param predicate2 The second atom predicate.
       * @param predicate3 The 3th atom predicate.
       */
      template<typename Predicate1, typename Predicate2, typename Predicate3>
      AtomAndPredicates(const Predicate1 &predicate1, const Predicate2 &predicate2,
          const Predicate3 &predicate3)
      {
        m_predicates.push_back(new Predicate1(predicate1));
        m_predicates.push_back(new Predicate2(predicate2));
        m_predicates.push_back(new Predicate3(predicate3));
      }

      /**
       * @brief Destructor.
       */
      ~AtomAndPredicates()
      {
        for (std::size_t i = 0; i < m_predicates.size(); ++i)
          delete m_predicates[i];
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param atom The atom to test.
       *
       * @return True if all predicates return true, otherwise false is returned.
       */
      bool operator()(const MoleculeType &mol, atom_type atom) const
      {
        for (std::size_t i = 0; i < m_predicates.size(); ++i)
          if (!m_predicates[i]->operator()(mol, atom))
            return false;
        return true;
      }

    private:
      std::vector<AtomPredicate<MoleculeType>*> m_predicates; //!< The atom predicates.
  };

  /**
   * @brief Utility function to create an AtomAndPredicate.
   *
   * @param mol The molecule (only used for type deduction).
   * @param predicate1 The 1st atom predicate.
   * @param predicate2 The 2nd atom predicate.
   *
   * @return The AtomAndPredicate.
   */
  template<typename MoleculeType, typename Predicate1, typename Predicate2>
  AtomAndPredicates<MoleculeType> atom_and_predicates(const MoleculeType &mol,
      const Predicate1 &predicate1, const Predicate2 &predicate2)
  {
    return AtomAndPredicates<MoleculeType>(predicate1, predicate2);
  }

  /**
   * @brief Utility function to create an AtomAndPredicate.
   *
   * @param mol The molecule (only used for type deduction).
   * @param predicate1 The 1st atom predicate.
   * @param predicate2 The 2nd atom predicate.
   * @param predicate3 The 3th atom predicate.
   *
   * @return The AtomAndPredicate.
   */
  template<typename MoleculeType, typename Predicate1, typename Predicate2, typename Predicate3>
  AtomAndPredicates<MoleculeType> atom_and_predicates(const MoleculeType &mol,
      const Predicate1 &predicate1, const Predicate2 &predicate2,
      const Predicate3 &predicate3)
  {
    return AtomAndPredicates<MoleculeType>(predicate1, predicate2, predicate3);
  }

  /**
   * @brief Class for combining atom predicates using a locgical or.
   */
  template<typename MoleculeType>
  class AtomOrPredicates : public AtomPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The atom type.
       */
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

      /**
       * @brief Constructor.
       *
       * @param predicate1 The first atom predicate.
       * @param predicate2 The second atom predicate.
       */
      template<typename Predicate1, typename Predicate2>
      AtomOrPredicates(const Predicate1 &predicate1, const Predicate2 &predicate2)
      {
        m_predicates.push_back(new Predicate1(predicate1));
        m_predicates.push_back(new Predicate2(predicate2));
      }

      /**
       * @brief Constructor.
       *
       * @param predicate1 The first atom predicate.
       * @param predicate2 The second atom predicate.
       * @param predicate3 The 3th atom predicate.
       */
      template<typename Predicate1, typename Predicate2, typename Predicate3>
      AtomOrPredicates(const Predicate1 &predicate1, const Predicate2 &predicate2,
          const Predicate3 &predicate3)
      {
        m_predicates.push_back(new Predicate1(predicate1));
        m_predicates.push_back(new Predicate2(predicate2));
        m_predicates.push_back(new Predicate3(predicate3));
      }

      /**
       * @brief Destructor.
       */
      ~AtomOrPredicates()
      {
        for (std::size_t i = 0; i < m_predicates.size(); ++i)
          delete m_predicates[i];
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param atom The atom to test.
       *
       * @return True if once a predicate return true, otherwise false is
       *         returned if all predicates return false.
       */
      bool operator()(const MoleculeType &mol, atom_type atom) const
      {
        for (std::size_t i = 0; i < m_predicates.size(); ++i)
          if (m_predicates[i]->operator()(mol, atom))
            return true;
        return false;
      }

    private:
      std::vector<AtomPredicate<MoleculeType>*> m_predicates; //!< The atom predicates.
  };

  /**
   * @brief Utility function to create an AtomOrPredicate.
   *
   * @param mol The molecule (only used for type deduction).
   * @param predicate1 The 1st atom predicate.
   * @param predicate2 The 2nd atom predicate.
   *
   * @return The AtomOrPredicate.
   */
  template<typename MoleculeType, typename Predicate1, typename Predicate2>
  AtomOrPredicates<MoleculeType> atom_or_predicates(const MoleculeType &mol,
      const Predicate1 &predicate1, const Predicate2 &predicate2)
  {
    return AtomOrPredicates<MoleculeType>(predicate1, predicate2);
  }

  /**
   * @brief Utility function to create an AtomOrPredicate.
   *
   * @param mol The molecule (only used for type deduction).
   * @param predicate1 The 1st atom predicate.
   * @param predicate2 The 2nd atom predicate.
   * @param predicate3 The 3th atom predicate.
   *
   * @return The AtomOrPredicate.
   */
  template<typename MoleculeType, typename Predicate1, typename Predicate2, typename Predicate3>
  AtomOrPredicates<MoleculeType> atom_and_predicates(const MoleculeType &mol,
      const Predicate1 &predicate1, const Predicate2 &predicate2,
      const Predicate3 &predicate3)
  {
    return AtomOrPredicates<MoleculeType>(predicate1, predicate2, predicate3);
  }

  /**
   * @brief Element atom predicate.
   *
   * This AtomPredicate is used for comparing the atom's element.
   */
  template<typename MoleculeType, template<typename> class Compare = std::equal_to>
  class ElementPredicate : public AtomPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The atom type.
       */
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

      /**
       * @brief Constructor.
       *
       * @param element The element number to compare with.
       */
      ElementPredicate(int element) : m_element(element)
      {
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param atom The atom to test.
       *
       * @return The result of comparing the atom's element with the element
       *         given in the constructor.
       */
      bool operator()(const MoleculeType &mol, atom_type atom) const
      {
        Compare<int> compare;
        return compare(get_element(mol, atom), m_element);
      }

    private:
      int m_element; //!< The element to compare with.
  };

  /**
   * @brief Utility function to create an ElementPredicate with std::equal_to compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param element The element to compare with.
   *
   * @return The ElementPredicate.
   */
  template<typename MoleculeType>
  ElementPredicate<MoleculeType, std::equal_to> element_eq_predicate(const MoleculeType &mol, int element)
  {
    return ElementPredicate<MoleculeType, std::equal_to>(element);
  }

  /**
   * @brief Utility function to create an ElementPredicate with std::less compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param element The element to compare with.
   *
   * @return The ElementPredicate.
   */
  template<typename MoleculeType>
  ElementPredicate<MoleculeType, std::less> element_lt_predicate(const MoleculeType &mol, int element)
  {
    return ElementPredicate<MoleculeType, std::less>(element);
  }

  /**
   * @brief Utility function to create an ElementPredicate with std::greater compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param element The element to compare with.
   *
   * @return The ElementPredicate.
   */
  template<typename MoleculeType>
  ElementPredicate<MoleculeType, std::greater> element_gt_predicate(const MoleculeType &mol, int element)
  {
    return ElementPredicate<MoleculeType, std::greater>(element);
  }

  /**
   * @brief Utility function to create an ElementPredicate with std::less_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param element The element to compare with.
   *
   * @return The ElementPredicate.
   */
  template<typename MoleculeType>
  ElementPredicate<MoleculeType, std::less_equal> element_leq_predicate(const MoleculeType &mol, int element)
  {
    return ElementPredicate<MoleculeType, std::less_equal>(element);
  }

  /**
   * @brief Utility function to create an ElementPredicate with std::greater_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param element The element to compare with.
   *
   * @return The ElementPredicate.
   */
  template<typename MoleculeType>
  ElementPredicate<MoleculeType, std::greater_equal> element_geq_predicate(const MoleculeType &mol, int element)
  {
    return ElementPredicate<MoleculeType, std::greater_equal>(element);
  }

  /**
   * @brief Mass atom predicate.
   *
   * This AtomPredicate is used for comparing the atom's mass.
   */
  template<typename MoleculeType, template<typename> class Compare = std::equal_to>
  class MassPredicate : public AtomPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The atom type.
       */
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

      /**
       * @brief Constructor.
       *
       * @param mass The mass number to compare with.
       */
      MassPredicate(int mass) : m_mass(mass)
      {
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param atom The atom to test.
       *
       * @return The result of comparing the atom's mass with the mass
       *         given in the constructor.
       */
      bool operator()(const MoleculeType &mol, atom_type atom) const
      {
        Compare<int> compare;
        return compare(get_mass(mol, atom), m_mass);
      }

    private:
      int m_mass; //!< The mass to compare with.
  };

  /**
   * @brief Utility function to create an MassPredicate with std::equal_to compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param mass The mass to compare with.
   *
   * @return The MassPredicate.
   */
  template<typename MoleculeType>
  MassPredicate<MoleculeType, std::equal_to> mass_eq_predicate(const MoleculeType &mol, int mass)
  {
    return MassPredicate<MoleculeType, std::equal_to>(mass);
  }

  /**
   * @brief Utility function to create an MassPredicate with std::less compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param mass The mass to compare with.
   *
   * @return The MassPredicate.
   */
  template<typename MoleculeType>
  MassPredicate<MoleculeType, std::less> mass_lt_predicate(const MoleculeType &mol, int mass)
  {
    return MassPredicate<MoleculeType, std::less>(mass);
  }

  /**
   * @brief Utility function to create an MassPredicate with std::greater compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param mass The mass to compare with.
   *
   * @return The MassPredicate.
   */
  template<typename MoleculeType>
  MassPredicate<MoleculeType, std::greater> mass_gt_predicate(const MoleculeType &mol, int mass)
  {
    return MassPredicate<MoleculeType, std::greater>(mass);
  }

  /**
   * @brief Utility function to create an MassPredicate with std::less_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param mass The mass to compare with.
   *
   * @return The MassPredicate.
   */
  template<typename MoleculeType>
  MassPredicate<MoleculeType, std::less_equal> mass_leq_predicate(const MoleculeType &mol, int mass)
  {
    return MassPredicate<MoleculeType, std::less_equal>(mass);
  }

  /**
   * @brief Utility function to create an MassPredicate with std::greater_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param mass The mass to compare with.
   *
   * @return The MassPredicate.
   */
  template<typename MoleculeType>
  MassPredicate<MoleculeType, std::greater_equal> mass_geq_predicate(const MoleculeType &mol, int mass)
  {
    return MassPredicate<MoleculeType, std::greater_equal>(mass);
  }

  /**
   * @brief Charge atom predicate.
   *
   * This AtomPredicate is used for comparing the atom's charge.
   */
  template<typename MoleculeType, template<typename> class Compare = std::equal_to>
  class ChargePredicate : public AtomPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The atom type.
       */
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

      /**
       * @brief Constructor.
       *
       * @param charge The charge to compare with.
       */
      ChargePredicate(int charge) : m_charge(charge)
      {
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param atom The atom to test.
       *
       * @return The result of comparing the atom's charge with the charge
       *         given in the constructor.
       */
      bool operator()(const MoleculeType &mol, atom_type atom) const
      {
        Compare<int> compare;
        return compare(get_charge(mol, atom), m_charge);
      }

    private:
      int m_charge; //!< The charge to compare with.
  };

  /**
   * @brief Utility function to create an ChargePredicate with std::equal_to compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param charge The charge to compare with.
   *
   * @return The ChargePredicate.
   */
  template<typename MoleculeType>
  ChargePredicate<MoleculeType, std::equal_to> charge_eq_predicate(const MoleculeType &mol, int charge)
  {
    return ChargePredicate<MoleculeType, std::equal_to>(charge);
  }

  /**
   * @brief Utility function to create an ChargePredicate with std::less compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param charge The charge to compare with.
   *
   * @return The ChargePredicate.
   */
  template<typename MoleculeType>
  ChargePredicate<MoleculeType, std::less> charge_lt_predicate(const MoleculeType &mol, int charge)
  {
    return ChargePredicate<MoleculeType, std::less>(charge);
  }

  /**
   * @brief Utility function to create an ChargePredicate with std::greater compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param charge The charge to compare with.
   *
   * @return The ChargePredicate.
   */
  template<typename MoleculeType>
  ChargePredicate<MoleculeType, std::greater> charge_gt_predicate(const MoleculeType &mol, int charge)
  {
    return ChargePredicate<MoleculeType, std::greater>(charge);
  }

  /**
   * @brief Utility function to create an ChargePredicate with std::less_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param charge The charge to compare with.
   *
   * @return The ChargePredicate.
   */
  template<typename MoleculeType>
  ChargePredicate<MoleculeType, std::less_equal> charge_leq_predicate(const MoleculeType &mol, int charge)
  {
    return ChargePredicate<MoleculeType, std::less_equal>(charge);
  }

  /**
   * @brief Utility function to create an ChargePredicate with std::greater_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param charge The charge to compare with.
   *
   * @return The ChargePredicate.
   */
  template<typename MoleculeType>
  ChargePredicate<MoleculeType, std::greater_equal> charge_geq_predicate(const MoleculeType &mol, int charge)
  {
    return ChargePredicate<MoleculeType, std::greater_equal>(charge);
  }

  /**
   * @brief Number of hydrogens atom predicate.
   *
   * This AtomPredicate is used for comparing the atom's number of hydrogens.
   */
  template<typename MoleculeType, template<typename> class Compare = std::equal_to>
  class NumHydrogensPredicate : public AtomPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The atom type.
       */
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

      /**
       * @brief Constructor.
       *
       * @param numHydrogens The number of hydrogens to compare with.
       */
      NumHydrogensPredicate(int numHydrogens) : m_numHydrogens(numHydrogens)
      {
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param atom The atom to test.
       *
       * @return The result of comparing the atom's number of hydrogens with
       *         the number of hydrogens given in the constructor.
       */
      bool operator()(const MoleculeType &mol, atom_type atom) const
      {
        Compare<int> compare;
        return compare(num_hydrogens(mol, atom), m_numHydrogens);
      }

    private:
      int m_numHydrogens; //!< The number of hydrogens to compare with.
  };

  /**
   * @brief Utility function to create an NumHydrogensPredicate with std::equal_to compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param numHydrogens The number of hydrogens to compare with.
   *
   * @return The NumHydrogensPredicate.
   */
  template<typename MoleculeType>
  NumHydrogensPredicate<MoleculeType, std::equal_to> num_hydrogens_eq_predicate(const MoleculeType &mol, int numHydrogens)
  {
    return NumHydrogensPredicate<MoleculeType, std::equal_to>(numHydrogens);
  }

  /**
   * @brief Utility function to create an NumHydrogensPredicate with std::less compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param numHydrogens The number of hydrogens to compare with.
   *
   * @return The NumHydrogensPredicate.
   */
  template<typename MoleculeType>
  NumHydrogensPredicate<MoleculeType, std::less> num_hydrogens_lt_predicate(const MoleculeType &mol, int numHydrogens)
  {
    return NumHydrogensPredicate<MoleculeType, std::less>(numHydrogens);
  }

  /**
   * @brief Utility function to create an NumHydrogensPredicate with std::greater compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param numHydrogens The number of hydrogens to compare with.
   *
   * @return The NumHydrogensPredicate.
   */
  template<typename MoleculeType>
  NumHydrogensPredicate<MoleculeType, std::greater> num_hydrogens_gt_predicate(const MoleculeType &mol, int numHydrogens)
  {
    return NumHydrogensPredicate<MoleculeType, std::greater>(numHydrogens);
  }

  /**
   * @brief Utility function to create an NumHydrogensPredicate with std::less_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param numHydrogens The number of hydrogens to compare with.
   *
   * @return The NumHydrogensPredicate.
   */
  template<typename MoleculeType>
  NumHydrogensPredicate<MoleculeType, std::less_equal> num_hydrogens_leq_predicate(const MoleculeType &mol, int numHydrogens)
  {
    return NumHydrogensPredicate<MoleculeType, std::less_equal>(numHydrogens);
  }

  /**
   * @brief Utility function to create an NumHydrogensPredicate with std::greater_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param numHydrogens The number of hydrogens to compare with.
   *
   * @return The NumHydrogensPredicate.
   */
  template<typename MoleculeType>
  NumHydrogensPredicate<MoleculeType, std::greater_equal> num_hydrogens_geq_predicate(const MoleculeType &mol, int numHydrogens)
  {
    return NumHydrogensPredicate<MoleculeType, std::greater_equal>(numHydrogens);
  }

  /**
   * @brief Degree atom predicate.
   *
   * This AtomPredicate is used for comparing the atom's degree.
   */
  template<typename MoleculeType, template<typename> class Compare = std::equal_to>
  class DegreePredicate : public AtomPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The atom type.
       */
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

      /**
       * @brief Constructor.
       *
       * @param degree The degree to compare with.
       */
      DegreePredicate(int degree) : m_degree(degree)
      {
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param atom The atom to test.
       *
       * @return The result of comparing the atom's degree with the degree
       *         given in the constructor.
       */
      bool operator()(const MoleculeType &mol, atom_type atom) const
      {
        Compare<int> compare;
        return compare(get_degree(mol, atom), m_degree);
      }

    private:
      int m_degree; //!< The degree to compare with.
  };

  /**
   * @brief Utility function to create an DegreePredicate with std::equal_to compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param degree The degree to compare with.
   *
   * @return The DegreePredicate.
   */
  template<typename MoleculeType>
  DegreePredicate<MoleculeType, std::equal_to> degree_eq_predicate(const MoleculeType &mol, int degree)
  {
    return DegreePredicate<MoleculeType, std::equal_to>(degree);
  }

  /**
   * @brief Utility function to create an DegreePredicate with std::less compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param degree The degree to compare with.
   *
   * @return The DegreePredicate.
   */
  template<typename MoleculeType>
  DegreePredicate<MoleculeType, std::less> degree_lt_predicate(const MoleculeType &mol, int degree)
  {
    return DegreePredicate<MoleculeType, std::less>(degree);
  }

  /**
   * @brief Utility function to create an DegreePredicate with std::greater compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param degree The degree to compare with.
   *
   * @return The DegreePredicate.
   */
  template<typename MoleculeType>
  DegreePredicate<MoleculeType, std::greater> degree_gt_predicate(const MoleculeType &mol, int degree)
  {
    return DegreePredicate<MoleculeType, std::greater>(degree);
  }

  /**
   * @brief Utility function to create an DegreePredicate with std::less_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param degree The degree to compare with.
   *
   * @return The DegreePredicate.
   */
  template<typename MoleculeType>
  DegreePredicate<MoleculeType, std::less_equal> degree_leq_predicate(const MoleculeType &mol, int degree)
  {
    return DegreePredicate<MoleculeType, std::less_equal>(degree);
  }

  /**
   * @brief Utility function to create an DegreePredicate with std::greater_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param degree The degree to compare with.
   *
   * @return The DegreePredicate.
   */
  template<typename MoleculeType>
  DegreePredicate<MoleculeType, std::greater_equal> degree_geq_predicate(const MoleculeType &mol, int degree)
  {
    return DegreePredicate<MoleculeType, std::greater_equal>(degree);
  }

  /**
   * @brief Aromatic atom predicate.
   *
   * This AtomPredicate is used for comparing the atom's aromaticity.
   */
  template<typename MoleculeType>
  class AromaticAtomPredicate : public AtomPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The atom type.
       */
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

      /**
       * @brief Constructor.
       *
       * @param aromatic The aromaticity to compare with.
       */
      AromaticAtomPredicate(bool aromatic) : m_aromatic(aromatic)
      {
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param atom The atom to test.
       *
       * @return The result of comparing the atom's aromaticity with the
       *         aromaticity given in the constructor.
       */
      bool operator()(const MoleculeType &mol, atom_type atom) const
      {
        return is_aromatic(mol, atom) == m_aromatic;
      }

    private:
      bool m_aromatic; //!< The aromaticity to compare with.
  };

  /**
   * @brief Utility function to create an AromaticAtomPredicate for aromatic atoms.
   *
   * @param mol The molecule (only used for type deduction).
   *
   * @return The AromaticAtomPredicate.
   */
  template<typename MoleculeType>
  AromaticAtomPredicate<MoleculeType> aromatic_atom_predicate(const MoleculeType &mol)
  {
    return AromaticAtomPredicate<MoleculeType>(true);
  }

  /**
   * @brief Utility function to create an AromaticAtomPredicate for non-aromatic atoms.
   *
   * @param mol The molecule (only used for type deduction).
   *
   * @return The AromaticAtomPredicate.
   */
  template<typename MoleculeType>
  AromaticAtomPredicate<MoleculeType> not_aromatic_atom_predicate(const MoleculeType &mol)
  {
    return AromaticAtomPredicate<MoleculeType>(false);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Bond Predicates
  //
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Base class for bond predicates.
   *
   * Bond predicates are used in combination with a number of functions
   * (e.g. molecule_has_bond(), atom_has_bond(), ...).
   * A bond predicate is a functor that returns true if a bond matches
   * the description of the predicate. There are a number of predicates
   * predefined (e.g. OrderPredicate, AromaticBondPredicate, ...). Additional
   * predicates can be added by inheriting BondPredicate and implementing
   * the pure virtual function operator.
   */
  template<typename MoleculeType>
  class BondPredicate
  {
    public:
      /**
       * @brief The bond type.
       */
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      /**
       * @brief Destructor.
       */
      virtual ~BondPredicate()
      {
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param bond The bond to test.
       *
       * @return True if the bond matches the predicate.
       */
      virtual bool operator()(const MoleculeType &mol, bond_type bond) const = 0;
  };



  /**
   * @brief Class for combining bond predicates using a locgical and.
   */
  template<typename MoleculeType>
  class BondAndPredicates : public BondPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The bond type.
       */
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      /**
       * @brief Constructor.
       *
       * @param predicate1 The 1st predicate.
       * @param predicate2 The 2nd predicate.
       */
      template<typename Predicate1, typename Predicate2>
      BondAndPredicates(const Predicate1 &predicate1, const Predicate2 &predicate2)
      {
        m_predicates.push_back(new Predicate1(predicate1));
        m_predicates.push_back(new Predicate2(predicate2));
      }

      /**
       * @brief Constructor.
       *
       * @param predicate1 The 1st bond predicate.
       * @param predicate2 The 2nd bond predicate.
       * @param predicate3 The 3th bond predicate.
       */
      template<typename Predicate1, typename Predicate2, typename Predicate3>
      BondAndPredicates(const Predicate1 &predicate1, const Predicate2 &predicate2,
          const Predicate3 &predicate3)
      {
        m_predicates.push_back(new Predicate1(predicate1));
        m_predicates.push_back(new Predicate2(predicate2));
        m_predicates.push_back(new Predicate3(predicate3));
      }

      /**
       * @brief Destructor.
       */
      ~BondAndPredicates()
      {
        for (std::size_t i = 0; i < m_predicates.size(); ++i)
          delete m_predicates[i];
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param bond The bond to test.
       *
       * @return True if a predicates return true, otherwise false is returned.
       */
      bool operator()(const MoleculeType &mol, bond_type bond) const
      {
        for (std::size_t i = 0; i < m_predicates.size(); ++i)
          if (!m_predicates->operator()(mol, bond))
            return false;
        return true;
      }

    private:
      std::vector<BondPredicate<MoleculeType>*> m_predicates; //!< The bond predicates
  };

  /**
   * @brief Utility function to create an BondAndPredicate.
   *
   * @param mol The molecule (only used for type deduction).
   * @param predicate1 The 1st bond predicate.
   * @param predicate2 The 2nd bond predicate.
   *
   * @return The BondAndPredicate.
   */
  template<typename MoleculeType, typename Predicate1, typename Predicate2>
  BondAndPredicates<MoleculeType> bond_and_predicates(const MoleculeType &mol,
      const Predicate1 &predicate1, const Predicate2 &predicate2)
  {
    return BondAndPredicates<MoleculeType>(predicate1, predicate2);
  }

  /**
   * @brief Utility function to create an BondAndPredicate.
   *
   * @param mol The molecule (only used for type deduction).
   * @param predicate1 The 1st bond predicate.
   * @param predicate2 The 2nd bond predicate.
   * @param predicate3 The 3th bond predicate.
   *
   * @return The BondAndPredicate.
   */
  template<typename MoleculeType, typename Predicate1, typename Predicate2, typename Predicate3>
  BondAndPredicates<MoleculeType> bond_and_predicates(const MoleculeType &mol,
      const Predicate1 &predicate1, const Predicate2 &predicate2,
      const Predicate3 &predicate3)
  {
    return BondAndPredicates<MoleculeType>(predicate1, predicate2, predicate3);
  }

  /**
   * @brief Class for combining bond predicates using a locgical or.
   */
  template<typename MoleculeType>
  class BondOrPredicates : public BondPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The bond type.
       */
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      /**
       * @brief Constructor.
       *
       * @param predicate1 The first bond predicate.
       * @param predicate2 The second bond predicate.
       */
      template<typename Predicate1, typename Predicate2>
      BondOrPredicates(const Predicate1 &predicate1, const Predicate2 &predicate2)
      {
        m_predicates.push_back(new Predicate1(predicate1));
        m_predicates.push_back(new Predicate2(predicate2));
      }

      /**
       * @brief Constructor.
       *
       * @param predicate1 The first bond predicate.
       * @param predicate2 The second bond predicate.
       * @param predicate3 The 3th bond predicate.
       */
      template<typename Predicate1, typename Predicate2, typename Predicate3>
      BondOrPredicates(const Predicate1 &predicate1, const Predicate2 &predicate2,
          const Predicate3 &predicate3)
      {
        m_predicates.push_back(new Predicate1(predicate1));
        m_predicates.push_back(new Predicate2(predicate2));
        m_predicates.push_back(new Predicate3(predicate3));
      }

      /**
       * @brief Destructor.
       */
      ~BondOrPredicates()
      {
        for (std::size_t i = 0; i < m_predicates.size(); ++i)
          delete m_predicates[i];
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param bond The bond to test.
       *
       * @return True if once a predicate return true, otherwise false is
       *         returned if all predicates return false.
       */
      bool operator()(const MoleculeType &mol, bond_type bond) const
      {
        for (std::size_t i = 0; i < m_predicates.size(); ++i)
          if (m_predicates[i]->operator()(mol, bond))
            return true;
        return false;
      }

    private:
      std::vector<BondPredicate<MoleculeType>*> m_predicates; //!< The bond predicates.
  };

  /**
   * @brief Utility function to create an BondOrPredicate.
   *
   * @param mol The molecule (only used for type deduction).
   * @param predicate1 The 1st bond predicate.
   * @param predicate2 The 2nd bond predicate.
   *
   * @return The BondOrPredicate.
   */
  template<typename MoleculeType, typename Predicate1, typename Predicate2>
  BondOrPredicates<MoleculeType> bond_or_predicates(const MoleculeType &mol,
      const Predicate1 &predicate1, const Predicate2 &predicate2)
  {
    return BondOrPredicates<MoleculeType>(predicate1, predicate2);
  }

  /**
   * @brief Utility function to create an BondOrPredicate.
   *
   * @param mol The molecule (only used for type deduction).
   * @param predicate1 The 1st bond predicate.
   * @param predicate2 The 2nd bond predicate.
   * @param predicate3 The 3th bond predicate.
   *
   * @return The BondOrPredicate.
   */
  template<typename MoleculeType, typename Predicate1, typename Predicate2, typename Predicate3>
  BondOrPredicates<MoleculeType> bond_and_predicates(const MoleculeType &mol,
      const Predicate1 &predicate1, const Predicate2 &predicate2,
      const Predicate3 &predicate3)
  {
    return BondOrPredicates<MoleculeType>(predicate1, predicate2, predicate3);
  }

  /**
   * @brief Order bond predicate.
   *
   * This BondPredicate is used for comparing the bond's order.
   */
  template<typename MoleculeType, template<typename> class Compare = std::equal_to>
  class OrderPredicate : public BondPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The bond type.
       */
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      /**
       * @brief Constructor.
       *
       * @param order The bond order to compare with.
       */
      OrderPredicate(int order) : m_order(order)
      {
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param bond The bond to test.
       *
       * @return The result of comparing the bond's order with the order
       *         given in the constructor.
       */
      bool operator()(const MoleculeType &mol, bond_type bond) const
      {
        Compare<int> compare;
        return compare(get_order(mol, bond), m_order);
      }

    private:
      int m_order; //!< The order to compare with.
  };

  /**
   * @brief Utility function to create an OrderPredicate with std::equal_to compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param order The bond order to compare with.
   *
   * @return The OrderPredicate.
   */
  template<typename MoleculeType>
  OrderPredicate<MoleculeType> order_eq_predicate(const MoleculeType &mol, int order)
  {
    return OrderPredicate<MoleculeType, std::equal_to>(order);
  }

  /**
   * @brief Utility function to create an OrderPredicate with std::less compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param order The bond order to compare with.
   *
   * @return The OrderPredicate.
   */
  template<typename MoleculeType>
  OrderPredicate<MoleculeType> order_lt_predicate(const MoleculeType &mol, int order)
  {
    return OrderPredicate<MoleculeType, std::less>(order);
  }

  /**
   * @brief Utility function to create an OrderPredicate with std::greater compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param order The bond order to compare with.
   *
   * @return The OrderPredicate.
   */
  template<typename MoleculeType>
  OrderPredicate<MoleculeType> order_gt_predicate(const MoleculeType &mol, int order)
  {
    return OrderPredicate<MoleculeType, std::greater>(order);
  }

  /**
   * @brief Utility function to create an OrderPredicate with std::less_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param order The bond order to compare with.
   *
   * @return The OrderPredicate.
   */
  template<typename MoleculeType>
  OrderPredicate<MoleculeType> order_leq_predicate(const MoleculeType &mol, int order)
  {
    return OrderPredicate<MoleculeType, std::less_equal>(order);
  }

  /**
   * @brief Utility function to create an OrderPredicate with std::greater_equal compare.
   *
   * @param mol The molecule (only used for type deduction).
   * @param order The bond order to compare with.
   *
   * @return The OrderPredicate.
   */
  template<typename MoleculeType>
  OrderPredicate<MoleculeType> order_geq_predicate(const MoleculeType &mol, int order)
  {
    return OrderPredicate<MoleculeType, std::greater_equal>(order);
  }

  /**
   * @brief Aromatic bond predicate.
   *
   * This BondPredicate is used for comparing the bond's aromaticity.
   */
  template<typename MoleculeType>
  class AromaticBondPredicate : public BondPredicate<MoleculeType>
  {
    public:
      /**
       * @brief The bond type.
       */
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      /**
       * @brief Constructor.
       *
       * @param aromatic The aromaticity to compare with.
       */
      AromaticBondPredicate(bool aromatic) : m_aromatic(aromatic)
      {
      }

      /**
       * @brief The predicate function.
       *
       * @param mol The molecule.
       * @param bond The bond to test.
       *
       * @return The result of comparing the bond's aromaticity with the
       *         aromaticity given in the constructor.
       */
      bool operator()(const MoleculeType &mol, bond_type bond) const
      {
        return is_aromatic(mol, bond) == m_aromatic;
      }

    private:
      bool m_aromatic; //!< The aromaticity to compare with.
  };

  /**
   * @brief Utility function to create an AromaticBondPredicate for aromatic bonds.
   *
   * @param mol The molecule (only used for type deduction).
   *
   * @return The AromaticBondPredicate.
   */
  template<typename MoleculeType>
  AromaticBondPredicate<MoleculeType> aromatic_bond_predicate(const MoleculeType &mol)
  {
    return AromaticBondPredicate<MoleculeType>(true);
  }

  /**
   * @brief Utility function to create an AromaticBondPredicate for non-aromatic bonds.
   *
   * @param mol The molecule (only used for type deduction).
   *
   * @return The AromaticBondPredicate.
   */
  template<typename MoleculeType>
  AromaticBondPredicate<MoleculeType> not_aromatic_bond_predicate(const MoleculeType &mol)
  {
    return AromaticBondPredicate<MoleculeType>(false);
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // Applying Predicates
  //
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Apply an atom predicate to each atom in a molecule.
   *
   * @param mol The molecule.
   * @param predicate The atom predicate.
   *
   * @return True once the predicate returns true for an atom.
   */
  template<typename MoleculeType, typename AtomPredicateType>
  bool molecule_has_atom(const MoleculeType &mol, const AtomPredicateType &predicate)
  {
    FOREACH_ATOM (atom, mol, MoleculeType)
      if (predicate(mol, *atom))
        return true;
    return false;
  }

  /**
   * @brief Get all atoms in a molecule for which the predicate returns true.
   *
   * @param mol The molecule.
   * @param predicate The atom predicate.
   *
   * @return A list of all atoms for which the predicate returns true.
   */
  template<typename MoleculeType, typename AtomPredicateType>
  std::vector<typename molecule_traits<MoleculeType>::atom_type>
  molecule_get_atoms(const MoleculeType &mol, const AtomPredicateType &predicate)
  {
    std::vector<typename molecule_traits<MoleculeType>::atom_type> atoms;
    FOREACH_ATOM (atom, mol, MoleculeType)
      if (predicate(mol, *atom))
        atoms.push_back(*atom);
    return atoms;
  }

  /**
   * @brief Apply a bond predicate to each bond in a molecule.
   *
   * @param mol The molecule.
   * @param predicate The bond predicate.
   *
   * @return True once the predicate returns true for a bond.
   */
  template<typename MoleculeType, typename BondPredicateType>
  bool molecule_has_bond(const MoleculeType &mol, const BondPredicateType &predicate)
  {
    FOREACH_BOND (bond, mol, MoleculeType)
      if (predicate(mol, *bond))
        return true;
    return false;
  }

  /**
   * @brief Get all bonds in a molecule for which the predicate returns true.
   *
   * @param mol The molecule.
   * @param predicate The bond predicate.
   *
   * @return A list of all bonds for which the predicate returns true.
   */
  template<typename MoleculeType, typename BondPredicateType>
  std::vector<typename molecule_traits<MoleculeType>::bond_type>
  molecule_get_bonds(const MoleculeType &mol, const BondPredicateType &predicate)
  {
    std::vector<typename molecule_traits<MoleculeType>::bond_type> bonds;
    FOREACH_BOND (bond, mol, MoleculeType)
      if (predicate(mol, *bond))
        bonds.push_back(*bond);
    return bonds;
  }

  /**
   * @brief Apply an atom predicate to each neighbor of an atom.
   *
   * @param mol The molecule.
   * @param atom The central atom.
   * @param predicate The atom predicate.
   *
   * @return True once the predicate returns true for a neighbor atom.
   */
  template<typename MoleculeType, typename AtomPredicateType>
  bool atom_has_nbr(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom, const AtomPredicateType &predicate)
  {
    FOREACH_NBR (nbr, atom, mol, MoleculeType)
      if (predicate(mol, *nbr))
        return true;
    return false;
  }

  /**
   * @brief Get all neighbor atoms for which the predicate returns true.
   *
   * @param mol The molecule.
   * @param atom The central atom.
   * @param predicate The atom predicate.
   *
   * @return A list of all neighbor atoms for which the predicate returns true.
   */
  template<typename MoleculeType, typename AtomPredicateType>
  std::vector<typename molecule_traits<MoleculeType>::atom_type>
  atom_get_nbrs(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom, const AtomPredicateType &predicate)
  {
    std::vector<typename molecule_traits<MoleculeType>::atom_type> atoms;
    FOREACH_NBR (nbr, atom, mol, MoleculeType)
      if (predicate(mol, *atom))
        atoms.push_back(*atom);
    return atoms;
  }

  /**
   * @brief Apply an atom predicate to each incident bond of an atom.
   *
   * @param mol The molecule.
   * @param atom The central atom.
   * @param predicate The atom predicate.
   *
   * @return True once the predicate returns true for an incident bond.
   */
  template<typename MoleculeType, typename BondPredicateType>
  bool atom_has_bond(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom, const BondPredicateType &predicate)
  {
    FOREACH_INCIDENT (bond, atom, mol, MoleculeType)
      if (predicate(mol, *bond))
        return true;
    return false;
  }

  /**
   * @brief Get all incident bonds for which the predicate returns true.
   *
   * @param mol The molecule.
   * @param atom The central atom.
   * @param predicate The bond predicate.
   *
   * @return A list of all incident bonds for which the predicate returns true.
   */
  template<typename MoleculeType, typename BondPredicateType>
  std::vector<typename molecule_traits<MoleculeType>::bond_type>
  atom_get_bonds(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom, const BondPredicateType &predicate)
  {
    std::vector<typename molecule_traits<MoleculeType>::bond_type> bonds;
    FOREACH_INCIDENT (bond, atom, mol, MoleculeType)
      if (predicate(mol, *bond))
        bonds.push_back(*bond);
    return bonds;
  }





  //////////////////////////////////////////////////////////////////////////////
  //
  // Lists of atoms and bonds
  //
  //////////////////////////////////////////////////////////////////////////////

  /**
   * @brief Check if a list of atoms/bonds contains an atom/bond with the specified index.
   *
   * @param mol The molecule.
   * @param objects The atom or bond list.
   * @param index The index to compare with.
   *
   * @return True if an atom/bond with the specified index is found in the list.
   */
  template<typename MoleculeType, typename AtomBondType>
  bool contains_index(const MoleculeType &mol, const std::vector<AtomBondType> &objects, Index index)
  {
    for (std::size_t i = 0; i < objects.size(); ++i)
      if (get_index(mol, objects[i]) == index)
        return true;
    return false;
  }

  /**
   * @brief Check if a list of atoms/bonds contains an atom/bond with a smaller index.
   *
   * @param mol The molecule.
   * @param objects The atom or bond list.
   * @param index The index to compare with.
   *
   * @return True if an atom/bond with a smaller index is found in the list.
   */
  template<typename MoleculeType, typename AtomBondType>
  bool contains_smaller_index(const MoleculeType &mol, const std::vector<AtomBondType> &objects, Index index)
  {
    for (std::size_t i = 0; i < objects.size(); ++i)
      if (get_index(mol, objects[i]) < index)
        return true;
    return false;
  }

  /**
   * @brief Check if a list of atoms/bonds contains an atom/bond with a larger index.
   *
   * @param mol The molecule.
   * @param objects The atom or bond list.
   * @param index The index to compare with.
   *
   * @return True if an atom/bond with a larger index is found in the list.
   */
  template<typename MoleculeType, typename AtomBondType>
  bool contains_larger_index(const MoleculeType &mol, const std::vector<AtomBondType> &objects, Index index)
  {
    for (std::size_t i = 0; i < objects.size(); ++i)
      if (get_index(mol, objects[i]) > index)
        return true;
    return false;
  }

  /**
   * @brief Get a string with the indices of the atoms/bonds in a list.
   *
   * @param mol The molecule.
   * @param objects The atom or bond list.
   *
   * @return A string with the indices of the atoms/bonds in a list.
   */
  template<typename MoleculeType, typename AtomBondType>
  std::string index_string(const MoleculeType &mol, const std::vector<AtomBondType> &objects)
  {
    std::stringstream ss;
    for (std::size_t i = 0; i < objects.size(); ++i)
      ss << get_index(mol, objects[i]) << " ";
    return ss.str();
  }

  //@}


}

#endif
