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
#ifndef HELIUM_RING_H
#define HELIUM_RING_H

#include <Helium/config.h>
#include <Helium/molecule.h>
#include <Helium/tie.h>
#include <Helium/util/vector.h>

#include <vector>
#include <istream>
#include <algorithm>
#include <cassert>

namespace Helium {

  /**
   * @brief Class representing a ring in a molecule.
   */
  template<typename MoleculeType>
  class Ring
  {
    public:
      /**
       * @brief The atom type.
       */
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      /**
       * @brief The bond type.
       */
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      /**
       * @brief Default constructor.
       */
      Ring() : m_mol(0)
      {
      }

      /**
       * @brief Constructor.
       *
       * The ring bonds will be determined based on the ring atoms.
       *
       * @param mol The molecule.
       * @param atoms The ring atoms.
       */
      Ring(const MoleculeType &mol, const std::vector<atom_type> &atoms)
        : m_mol(&mol), m_atoms(atoms)
      {
        bondsFromAtoms();
      }

      /**
       * @brief Constructor.
       *
       * The ring atoms will be determined based on the ring bonds.
       *
       * @param mol The molecule.
       * @param bonds The ring bonds.
       */
      Ring(const MoleculeType &mol, const std::vector<bond_type> &bonds)
        : m_mol(&mol), m_bonds(bonds)
      {
        atomsFromBonds();
      }

      /**
       * @brief Constructor.
       *
       * @param mol The molecule.
       * @param atoms The ring atoms.
       * @param bonds The ring bonds.
       */
      Ring(const MoleculeType &mol, const std::vector<atom_type> &atoms,
          const std::vector<bond_type> &bonds) : m_mol(&mol), m_atoms(atoms),
          m_bonds(bonds)
      {
      }

      /**
       * @brief Get the ring size.
       *
       * This is the number of atoms or bonds in the ring.
       *
       * @return The ring size.
       */
      std::size_t size() const
      {
        return m_atoms.size();
      }

      /**
       * @brief Get the ring atoms.
       *
       * @return The ring atoms.
       */
      const std::vector<atom_type>& atoms() const
      {
        return m_atoms;
      }

      /**
       * @brief Get the ring atom with the specified index.
       *
       * @param index The index of the atom to get.
       *
       * @return The ring atoms.
       */
      atom_type atom(std::size_t index) const
      {
        return m_atoms[index];
      }

      /**
       * @brief Get the ring bonds.
       *
       * @return The ring bonds.
       */
      const std::vector<bond_type>& bonds() const
      {
        return m_bonds;
      }

      /**
       * @brief Get the ring bond with the specified index.
       *
       * @param index The index of the bond to get.
       *
       * @return The ring bonds.
       */
      bond_type bond(std::size_t index) const
      {
        return m_bonds[index];
      }

      /**
       * @brief Check if the ring contains the specified atom.
       *
       * @param atom The atom to check for.
       *
       * @return True if the ring contains the specified atom.
       */
      bool containsAtom(atom_type atom) const
      {
        assert(m_mol);
        return contains_index(*m_mol, m_atoms, get_index(*m_mol, atom));
      }

      /**
       * @brief Check if the ring contains the specified bond.
       *
       * @param bond The bond to check for.
       *
       * @return True if the ring contains the specified bond.
       */
      bool containsBond(bond_type bond) const
      {
        assert(m_mol);
        return contains_index(*m_mol, m_bonds, get_index(*m_mol, bond));
      }

    private:
      /**
       * @brief Determine ring bonds based on ring atoms.
       */
      void bondsFromAtoms()
      {
        for (std::size_t i = 1; i < m_atoms.size(); ++i) {
          bond_type bond = get_bond(*m_mol, m_atoms[i - 1], m_atoms[i]);
          assert(bond != molecule_traits<MoleculeType>::null_bond());
          m_bonds.push_back(bond);
        }
        bond_type bond = get_bond(*m_mol, m_atoms[m_atoms.size() - 1], m_atoms[0]);
        assert(bond != molecule_traits<MoleculeType>::null_bond());
        m_bonds.push_back(bond);
      }

      /**
       * @brief Determine ring atoms based on ring bonds.
       */
      void atomsFromBonds()
      {
        for (std::size_t i = 0; i < m_bonds.size(); ++i) {
          atom_type source = get_source(*m_mol, m_bonds[i]);
          if (!contains_index(*m_mol, m_atoms, source))
            m_atoms.push_back(source);
          atom_type target = get_target(*m_mol, m_bonds[i]);
          if (!contains_index(*m_mol, m_atoms, target))
            m_atoms.push_back(target);
        }
      }

      const MoleculeType *m_mol; //!< The molecule.
      std::vector<atom_type> m_atoms; //!< The ring atoms.
      std::vector<bond_type> m_bonds; //!< The ring bonds.
  };

  /**
   * @brief Class representing a set of rings in a molecule.
   */
  template<typename MoleculeType>
  class RingSet
  {
    public:
      /**
       * @brief The atom type.
       */
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      /**
       * @brief The bond type.
       */
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
      /**
       * @brief The ring type.
       */
      typedef Ring<MoleculeType> ring_type;
      /**
       * @brief The ring iterator type.
       */
      typedef typename std::vector<ring_type>::iterator ring_iter;
      /**
       * @brief The ring constant iterator type.
       */
      typedef typename std::vector<ring_type>::const_iterator const_ring_iter;

      /**
       * @brief Constructor.
       *
       * @param mol The molecule.
       */
      RingSet(const MoleculeType &mol) : m_mol(mol)
      {
        m_cyclicAtoms.resize(num_atoms(m_mol));
        m_cyclicBonds.resize(num_bonds(m_mol));
      }

      /**
       * @brief Get the number of rings in the set.
       *
       * @return The number of rings in the set.
       */
      std::size_t size() const
      {
        return m_rings.size();
      }

      /**
       * @brief Get the rings in the set.
       *
       * @return The rings in the set.
       */
      const std::vector<ring_type>& rings() const
      {
        return m_rings;
      }

      /**
       * @brief Get the ring with the specified index.
       *
       * @param index The index of the ring to get.
       *
       * @return The ring with the specified index.
       */
      const ring_type& ring(std::size_t index) const
      {
        return m_rings[index];
      }

      /**
       * @brief Check if the specified atom is in a ring.
       *
       * @param atom The atom to check.
       *
       * @return True if the atom is in a ring.
       */
      bool isAtomInRing(atom_type atom) const
      {
        return m_cyclicAtoms[get_index(m_mol, atom)];
      }

      /**
       * @brief Check if the specified bond is in a ring.
       *
       * @param bond The bond to check.
       *
       * @return True if the bond is in a ring.
       */
      bool isBondInRing(bond_type bond) const
      {
        return m_cyclicBonds[get_index(m_mol, bond)];
      }

      /**
       * @brief Check if the specified atom is in a ring of a certain size.
       *
       * @param atom The atom to check.
       * @param ringSize The ring size.
       *
       * @return True if the atom is in a ring of size @p ringSize.
       */
      bool isAtomInRingSize(atom_type atom, int ringSize) const
      {
        for (const_ring_iter r = m_rings.begin(); r != m_rings.end(); ++r)
          if (r->size() == ringSize && r->containsAtom(atom))
            return true;
        return false;
      }

      /**
       * @brief Check if the specified bond is in a ring of a certain size.
       *
       * @param bond The bond to check.
       * @param ringSize The ring size.
       *
       * @return True if the bond is in a ring of size @p ringSize.
       */
      bool isBondInRingSize(bond_type bond, int ringSize) const
      {
        for (const_ring_iter r = m_rings.begin(); r != m_rings.end(); ++r)
          if (r->size() == ringSize && r->containsBond(bond))
            return true;
        return false;
      }

      /**
       * @brief Get the number of incident ring bonds an atom has.
       *
       * @param atom The atom.
       *
       * @return The number of incident ring bonds.
       */
      int numRingBonds(atom_type atom) const
      {
        int result = 0;
        FOREACH_INCIDENT (bond, atom, m_mol, MoleculeType)
          if (isBondInRing(*bond))
            ++result;
        return result;
      }

      /**
       * @brief Get the number of neighbors an atom has that are in a ring.
       *
       * @param atom The atom.
       *
       * @return The number of ring neighbors.
       */
      int numRingNbrs(atom_type atom) const
      {
        int result = 0;
        FOREACH_NBR (nbr, atom, m_mol, MoleculeType)
          if (isAtomInRing(*nbr))
            ++result;
        return result;
      }

      /**
       * @brief Get the number of rings that the specified atom is part of.
       *
       * @param atom The atom.
       *
       * @return The number of rings that contain the atom.
       */
      int numRings(atom_type atom) const
      {
        int count = 0;
        for (std::size_t i = 0; i < m_rings.size(); ++i)
          if (m_rings[i].containsAtom(atom))
            ++count;
        return count;
      }

      /**
       * @brief Add a ring to the set.
       *
       * @param ring The ring to add.
       */
      void addRing(const ring_type &ring)
      {
        m_rings.push_back(ring);
        // update cyclic atoms/bonds
        for (std::size_t i = 0; i < ring.size(); ++i) {
          m_cyclicAtoms[get_index(m_mol, ring.atom(i))] = true;
          m_cyclicBonds[get_index(m_mol, ring.bond(i))] = true;
        }
      }

    private:
      const MoleculeType &m_mol; //!< The molecule.
      std::vector<ring_type> m_rings; //!< The rings.
      std::vector<bool> m_cyclicAtoms; //!< The cyclic atoms.
      std::vector<bool> m_cyclicBonds; //!< The cyclic bonds.
  };

}

#endif
