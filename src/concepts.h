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
#ifndef HELIUM_CONCEPTS_H
#define HELIUM_CONCEPTS_H

#include <Helium/molecule.h>

#include <boost/concept_check.hpp>

namespace Helium {

  template<typename MoleculeType>
  class AtomConcept
  {
    public:
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      typedef typename molecule_traits<MoleculeType>::incident_iter incident_iter;
      typedef typename molecule_traits<MoleculeType>::nbr_iter nbr_iter;

      BOOST_CONCEPT_ASSERT((boost::EqualityComparable<atom_type>));
      BOOST_CONCEPT_ASSERT((boost::LessThanComparable<atom_type>));

      BOOST_CONCEPT_USAGE(AtomConcept)
      {
        index = get_index(mol, atom);
        bonds = get_bonds(mol, atom);
        nbrs = get_nbrs(mol, atom);
        aromatic = is_aromatic(mol, atom);
        integer = get_element(mol, atom);
        integer = get_mass(mol, atom);
        integer = get_degree(mol, atom);
        integer = get_hydrogens(mol, atom);
        integer = get_charge(mol, atom);

        same_type(molecule_traits<MoleculeType>::null_bond(), *bonds.begin());
        same_type(molecule_traits<MoleculeType>::null_atom(), *nbrs.begin());
      }

      template<typename T>
      void same_type(T, T);

    private:
      MoleculeType mol;
      atom_type atom;
      iterator_pair<incident_iter> bonds;
      iterator_pair<nbr_iter> nbrs;
      Index index;
      bool aromatic;
      int integer;
  };

  template<typename MoleculeType>
  class BondConcept
  {
    public:
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      BOOST_CONCEPT_ASSERT((boost::EqualityComparable<bond_type>));
      BOOST_CONCEPT_ASSERT((boost::LessThanComparable<bond_type>));

      BOOST_CONCEPT_USAGE(BondConcept)
      {
        index = get_index(mol, bond);
        atom = get_source(mol, bond);
        atom = get_target(mol, bond);
        atom = get_other(mol, bond, atom);
        aromatic = is_aromatic(mol, bond);
        integer = get_order(mol, bond);
      }

    private:
      MoleculeType mol;
      atom_type atom;
      bond_type bond;
      Index index;
      bool aromatic;
      int integer;
  };

  template<typename MoleculeType>
  class MoleculeConcept /*: boost::Assignable<MoleculeType>*/
  {
    public:
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;
      typedef typename molecule_traits<MoleculeType>::bond_iter bond_iter;
      typedef typename molecule_traits<MoleculeType>::incident_iter incident_iter;
      typedef typename molecule_traits<MoleculeType>::nbr_iter nbr_iter;

      BOOST_CONCEPT_ASSERT((AtomConcept<MoleculeType>));
      BOOST_CONCEPT_ASSERT((BondConcept<MoleculeType>));

      BOOST_CONCEPT_ASSERT((boost::ForwardIterator<atom_iter>));
      BOOST_CONCEPT_ASSERT((boost::ForwardIterator<bond_iter>));
      BOOST_CONCEPT_ASSERT((boost::ForwardIterator<incident_iter>));
      BOOST_CONCEPT_ASSERT((boost::ForwardIterator<nbr_iter>));

      BOOST_CONCEPT_USAGE(MoleculeConcept)
      {
        atom = molecule_traits<MoleculeType>::null_atom();
        bond = molecule_traits<MoleculeType>::null_bond();
        index = molecule_traits<MoleculeType>::null_index();

        numAtoms = num_atoms(mol);
        numBonds = num_bonds(mol);

        atoms = get_atoms(mol);
        bonds = get_bonds(mol);

        atom = get_atom(mol, 0);
        bond = get_bond(mol, 0);

        same_type(molecule_traits<MoleculeType>::null_atom(), *atoms.begin());
        same_type(molecule_traits<MoleculeType>::null_bond(), *bonds.begin());
      }

      template<typename T>
      void same_type(T, T);

    private:
      MoleculeType mol;
      atom_type atom;
      bond_type bond;
      iterator_pair<atom_iter> atoms;
      iterator_pair<bond_iter> bonds;
      Size numAtoms;
      Size numBonds;
      Index index;
  };

  template<typename MoleculeType>
  void check_molecule_concept(const MoleculeType &mol)
  {
    BOOST_CONCEPT_ASSERT((MoleculeConcept<MoleculeType>));

    // the functions below are forward declared and are called to ensure
    // they are defined
    num_atoms(mol);
    num_bonds(mol);
    get_atoms(mol);
    get_bonds(mol);
    get_atom(mol, 0);
    get_bond(mol, 0);

    typename molecule_traits<MoleculeType>::atom_type atom;
    get_index(mol, atom);
    get_bonds(mol, atom);
    get_nbrs(mol, atom);
    is_aromatic(mol, atom);
    get_element(mol, atom);
    get_mass(mol, atom);
    get_degree(mol, atom);
    get_hydrogens(mol, atom);
    get_charge(mol, atom);

    typename molecule_traits<MoleculeType>::bond_type bond;
    get_index(mol, bond);
    get_source(mol, bond);
    get_target(mol, bond);
    get_other(mol, bond, atom);
    is_aromatic(mol, bond);
    get_order(mol, bond);
  }

  template<typename MoleculeType>
  class EditableAtomConcept : public AtomConcept<MoleculeType>
  {
    public:
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

      BOOST_CONCEPT_USAGE(EditableAtomConcept)
      {
        set_aromatic(mol, atom, true);
        set_element(mol, atom, 6);
        set_mass(mol, atom, 12);
        set_hydrogens(mol, atom, 4);
        set_charge(mol, atom, 0);
      }

    private:
      MoleculeType mol;
      atom_type atom;
  };

  template<typename MoleculeType>
  class EditableBondConcept : public BondConcept<MoleculeType>
  {
    public:
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      BOOST_CONCEPT_USAGE(EditableBondConcept)
      {
        set_aromatic(mol, bond, true);
        set_order(mol, bond, 1);
      }

    private:
      MoleculeType mol;
      bond_type bond;
  };

  template<typename MoleculeType>
  class EditableMoleculeConcept : public MoleculeConcept<MoleculeType>
  {
    public:
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      BOOST_CONCEPT_ASSERT((EditableAtomConcept<MoleculeType>));
      BOOST_CONCEPT_ASSERT((EditableBondConcept<MoleculeType>));

      BOOST_CONCEPT_USAGE(EditableMoleculeConcept)
      {
        atom = add_atom(mol);
        bond = add_bond(mol, atom, atom);

        remove_atom(mol, atom);
        remove_bond(mol, bond);

        clear_molecule(mol);
      }

    private:
      MoleculeType mol;
      atom_type atom;
      bond_type bond;
  };

  template<typename MoleculeType>
  void check_editable_molecule_concept(MoleculeType &mol)
  {
    BOOST_CONCEPT_ASSERT((EditableMoleculeConcept<MoleculeType>));
    check_molecule_concept(mol);

    // the functions below are forward declared and are called to ensure
    // they are defined
    typename molecule_traits<MoleculeType>::atom_type atom = add_atom(mol);
    typename molecule_traits<MoleculeType>::atom_type source = add_atom(mol);
    typename molecule_traits<MoleculeType>::atom_type target = add_atom(mol);
    typename molecule_traits<MoleculeType>::bond_type bond = add_bond(mol, source, target);

    set_aromatic(mol, atom, true);
    set_element(mol, atom, 6);
    set_mass(mol, atom, 12);
    set_hydrogens(mol, atom, 4);
    set_charge(mol, atom, 0);

    set_aromatic(mol, bond, true);
    set_order(mol, bond, 1);

    remove_bond(mol, bond);
    remove_atom(mol, atom);
    clear_molecule(mol);
  }

}

#endif
