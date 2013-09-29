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

#define FOREACH_ATOM(atom, mol, molecule_type) \
  typename molecule_traits<molecule_type>::atom_iter atom, end_##atom##__; \
  TIE(atom, end_##atom##__) = get_atoms(mol); \
  for (; atom != end_##atom##__; ++atom)

#define FOREACH_BOND(bond, mol, molecule_type) \
  typename molecule_traits<molecule_type>::bond_iter bond, end_##bond##__; \
  TIE(bond, end_##bond##__) = get_bonds(mol); \
  for (; bond != end_##bond##__; ++bond)

#define FOREACH_NBR(nbr, atom, mol, molecule_type) \
  typename molecule_traits<molecule_type>::nbr_iter nbr, end_##nbr##__; \
  TIE(nbr, end_##nbr##__) = get_nbrs(mol, atom); \
  for (; nbr != end_##nbr##__; ++nbr)

#define FOREACH_INCIDENT(bond, atom, mol, molecule_type) \
  typename molecule_traits<molecule_type>::incident_iter bond, end_##bond##__; \
  TIE(bond, end_##bond##__) = get_bonds(mol, atom); \
  for (; bond != end_##bond##__; ++bond)


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
  std::pair<typename molecule_traits<MoleculeType>::atom_iter, typename molecule_traits<MoleculeType>::atom_iter>
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
  typename molecule_traits<MoleculeType>::atom_type get_atom(const MoleculeType &mol, Index index);

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
  std::pair<typename molecule_traits<MoleculeType>::bond_iter, typename molecule_traits<MoleculeType>::bond_iter>
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
  typename molecule_traits<MoleculeType>::bond_type get_bond(const MoleculeType &mol, Index index);

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
  typename molecule_traits<MoleculeType>::bond_type get_bond(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type source,
                                                                                      typename molecule_traits<MoleculeType>::atom_type target);

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
  Index get_index(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom);

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
  std::pair<typename molecule_traits<MoleculeType>::incident_iter, typename molecule_traits<MoleculeType>::incident_iter>
  get_bonds(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom);

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
  std::pair<typename molecule_traits<MoleculeType>::nbr_iter, typename molecule_traits<MoleculeType>::nbr_iter>
  get_nbrs(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom);

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
  bool is_aromatic(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom);

  /**
   * @brief Get the atom's cyclic flag.
   *
   * @pre The atom must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_atom()).
   *
   * @param mol The molecule.
   * @param atom The atom.
   *
   * @return True if the atom is cyclic, false otherwise.
   */
  template<typename MoleculeType>
  bool is_cyclic(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom);

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
  int get_element(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom);

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
  int get_mass(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom);

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
  int get_degree(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom);

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
  int num_hydrogens(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom);

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
  int get_charge(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom);

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
   * @brief Get the bond's cyclic flag.
   *
   * @pre The bond must be valid (i.e. not equal to molecule_traits<MoleculeType>::null_bond()).
   *
   * @param mol The molecule.
   * @param bond The bond.
   *
   * @return True if the bond is cyclic, false otherwise.
   */
  template<typename MoleculeType>
  bool is_cyclic(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond);

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
  bool get_order(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond);

  //@}

}

#endif
