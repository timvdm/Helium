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

  template<typename MoleculeType>
  class Ring
  {
    public:
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      Ring() : m_mol(0)
      {
      }

      Ring(const MoleculeType &mol, const std::vector<atom_type> &atoms)
        : m_mol(&mol), m_atoms(atoms)
      {
        bondsFromAtoms();
      }

      Ring(const MoleculeType &mol, const std::vector<bond_type> &bonds)
        : m_mol(&mol), m_bonds(bonds)
      {
        atomsFromBonds();
      }

      Ring(const MoleculeType &mol, const std::vector<atom_type> &atoms,
          const std::vector<bond_type> &bonds) : m_mol(&mol), m_atoms(atoms),
          m_bonds(bonds)
      {
      }

      std::size_t size() const
      {
        return m_atoms.size();
      }

      const std::vector<atom_type>& atoms() const
      {
        return m_atoms;
      }

      atom_type atom(std::size_t index) const
      {
        return m_atoms[index];
      }

      const std::vector<bond_type>& bonds() const
      {
        return m_bonds;
      }

      bond_type bond(std::size_t index) const
      {
        return m_bonds[index];
      }

      bool containsAtom(atom_type atom) const
      {
        assert(m_mol);
        return contains_index(*m_mol, m_atoms, get_index(*m_mol, atom));
      }

      bool containsBond(bond_type bond) const
      {
        assert(m_mol);
        return contains_index(*m_mol, m_bonds, get_index(*m_mol, bond));
      }


    private:
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

      const MoleculeType *m_mol;
      std::vector<atom_type> m_atoms;
      std::vector<bond_type> m_bonds;
  };

  template<typename MoleculeType>
  class RingSet
  {
    public:
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
      typedef Ring<MoleculeType> ring_type;
      typedef typename std::vector<ring_type>::iterator ring_iter;
      typedef typename std::vector<ring_type>::const_iterator const_ring_iter;

      RingSet(const MoleculeType &mol) : m_mol(mol)
      {
        m_cyclicAtoms.resize(num_atoms(m_mol));
        m_cyclicBonds.resize(num_bonds(m_mol));
      }

      std::size_t size() const
      {
        return m_rings.size();
      }

      const std::vector<ring_type>& rings() const
      {
        return m_rings;
      }

      const ring_type& ring(std::size_t index) const
      {
        return m_rings[index];
      }

      bool isAtomInRing(atom_type atom) const
      {
        return m_cyclicAtoms[get_index(m_mol, atom)];
      }

      bool isBondInRing(bond_type bond) const
      {
        return m_cyclicBonds[get_index(m_mol, bond)];
      }

      bool isAtomInRingSize(atom_type atom, int ringSize) const
      {
        for (const_ring_iter r = m_rings.begin(); r != m_rings.end(); ++r)
          if (r->size() == ringSize && r->containsAtom(atom))
            return true;
        return false;
      }

      bool isBondInRingSize(bond_type bond, int ringSize) const
      {
        for (const_ring_iter r = m_rings.begin(); r != m_rings.end(); ++r)
          if (r->size() == ringSize && r->containsBond(bond))
            return true;
        return false;
      }

      int numRingBonds(atom_type atom) const
      {
        int result = 0;
        FOREACH_INCIDENT (bond, atom, m_mol, MoleculeType)
          if (isBondInRing(*bond))
            ++result;
        return result;
      }

      int numRingNbrs(atom_type atom) const
      {
        int result = 0;
        FOREACH_NBR (nbr, atom, m_mol, MoleculeType)
          if (isAtomInRing(*nbr))
            ++result;
        return result;
      }

      int numRings(atom_type atom) const
      {
        int count = 0;
        for (std::size_t i = 0; i < m_rings.size(); ++i)
          if (m_rings[i].containsAtom(atom))
            ++count;
        return count;
      }

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
      const MoleculeType &m_mol;
      std::vector<ring_type> m_rings;
      std::vector<bool> m_cyclicAtoms;
      std::vector<bool> m_cyclicBonds;
  };

}

#endif
