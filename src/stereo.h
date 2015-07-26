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
#ifndef HELIUM_STEREOUNIT_H
#define HELIUM_STEREOUNIT_H

#include <Helium/molecule.h>
#include <Helium/ring.h>

#include <vector>
#include <limits>
#include <cassert>

namespace Helium {

  struct Stereo
  {
    typedef unsigned int Ref;

    enum Type {
      None,
      Tetrahedral, //!< 4 ligands around center atom
      TrigonalBipyramidal, //!< 5 ligands, two axial with the remaning on a plane
      Octahedral, //!< 6 ligands, 3 axes
      CisTrans, //!< cis/trans double bond stereochemistry
      SquarePlanar, //!< 4 ligands on a plane
      Allene, //!< 4 ligands at the end of even number of double bonds
    };

    enum Class {
      TH1, TH2,
      AL1, AL2,
      SP1, SP2, SP3,
      TB1, TB2, TB3, TB4, TB5, TB6, TB7, TB8, TB9,
      TB10, TB11, TB12, TB13, TB14, TB15, TB16, TB17, TB18, TB19, TB20,
      OH1, OH2, OH3, OH4, OH5, OH6, OH7, OH8, OH9, OH10,
      OH11, OH12, OH13, OH14, OH15, OH16, OH17, OH18, OH19, OH20,
      OH21, OH22, OH23, OH24, OH25, OH26, OH27, OH28, OH29, OH30
    };

    static Ref nullRef()
    {
      return std::numeric_limits<Ref>::max();
    }

    static Ref implRef()
    {
      return std::numeric_limits<Ref>::max() - 1;
    }
  };

  /**
   *
   * Tetrahedral
   *   1 center, 4 refs: a b c d, viewing from ref a, refs b c d are clockwise
   *
   *        b  c
   *        | /
   *   a -- c1
   *         \
   *          d
   *
   * CisTrans
   *   2 centers, 4 refs: a b c d
   *
   *   a          d
   *    \        /
   *     c1 == c2
   *    /        \
   *   b          c
   *
   * SquarePlanar
   *   1 center, 4 refs: a b c d
   *
   *   a     d
   *    \   /
   *     c1
   *    /  \
   *   b    c
   *
   *
   *
   */
  class StereoStorage
  {
    public:
      StereoStorage(enum Stereo::Type type = Stereo::None) : m_type(type)
      {
        m_center = Stereo::nullRef();
        std::fill(m_refs, m_refs + 6, Stereo::nullRef());
      }

      ~StereoStorage()
      {
      }

      bool isValid() const
      {
        // check type
        if (m_type == Stereo::None)
          return false;

        // check center(s)
        switch (m_type) {
          default:
            if (m_center == Stereo::nullRef() || m_center == Stereo::implRef())
              return false;
        }

        // check refs
        int nullCount = std::count(m_refs, m_refs + numRefs(), Stereo::nullRef());
        if (nullCount)
          return false;

        int implCount = std::count(m_refs, m_refs + numRefs(), Stereo::implRef());
        switch (m_type) {
          case Stereo::Tetrahedral:
            if (implCount > 1)
              return false;
            break;
          case Stereo::Allene:
            break;
          case Stereo::CisTrans:
            // each double bond atom needs at least one explicit ref
            if (m_refs[0] == Stereo::implRef() && m_refs[1] == Stereo::implRef())
              return false;
            if (m_refs[2] == Stereo::implRef() && m_refs[3] == Stereo::implRef())
              return false;
            break;
          case Stereo::SquarePlanar:
            if (implCount > 3)
              return false;
            break;
          case Stereo::TrigonalBipyramidal:
          case Stereo::Octahedral:
            // invalid if both axis atoms are hydrogens
            //if (m_refs[0] == Stereo::implRef() && m_refs[4] == Stereo::implRef())
            //  return false;
            // invalid if 2 of the 3 plane atoms are hydrogens
            //if (std::count(m_refs + 1, m_refs + 4, Stereo::implRef()) > 1)
            //  return false;
            break;
          default:
            break;
        }

        return true;
      }

      enum Stereo::Type type() const
      {
        return m_type;
      }

      int numCenters() const
      {
        switch (m_type) {
          case Stereo::CisTrans:
          case Stereo::Allene:
            return 2;
          default:
            return 1;
        }
      }

      int numRefs() const
      {
        switch (m_type) {
          case Stereo::TrigonalBipyramidal:
            return 5;
          case Stereo::Octahedral:
            return 6;
          default:
            return 4;
        }
      }

      Stereo::Ref center() const
      {
        return m_center;
      }

      Stereo::Ref ref(int index) const
      {
        assert(index < numRefs());
        return m_refs[index];
      }

      const Stereo::Ref* refs() const
      {
        return m_refs;
      }

      template<typename RefIter>
      static StereoStorage create(enum Stereo::Type type, Stereo::Ref center,
          RefIter beginRefs, RefIter endRefs)
      {
        StereoStorage stereo(type);
        stereo.m_center = center;
        std::copy(beginRefs, endRefs, stereo.m_refs);
        return stereo;
      }

    private:
      enum Stereo::Type m_type; //!< The stereo type
      Stereo::Ref m_center; //!< The stereo center atom or bond
      Stereo::Ref m_refs[6]; //!< The 4 stereo
  };

  std::ostream& operator<<(std::ostream &os, const StereoStorage &stereo);

  class Stereochemistry
  {
    public:
      Stereochemistry()
      {
      }

      Stereochemistry(const Stereochemistry &other)
      {
        copy(other);
      }

      std::size_t numStereo() const
      {
        return m_stereo.size();
      }

      std::size_t numTetrahedral() const
      {
        return numType(Stereo::Tetrahedral);
      }

      std::size_t numAllene() const
      {
        return numType(Stereo::Allene);
      }

      std::size_t numCisTrans() const
      {
        return numType(Stereo::CisTrans);
      }

      std::size_t numSquarePlanar() const
      {
        return numType(Stereo::SquarePlanar);
      }

      std::size_t numTrigonalBipyramidal() const
      {
        return numType(Stereo::TrigonalBipyramidal);
      }

      std::size_t numOctahedral() const
      {
        return numType(Stereo::Octahedral);
      }

      bool isTetrahedral(Stereo::Ref center) const
      {
        return isType(Stereo::Tetrahedral, center);
      }

      bool isAllene(Stereo::Ref center) const
      {
        return isType(Stereo::Allene, center);
      }

      bool isCisTrans(Stereo::Ref bond) const
      {
        return isType(Stereo::CisTrans, bond);
      }

      bool isSquarePlanar(Stereo::Ref center) const
      {
        return isType(Stereo::SquarePlanar, center);
      }

      bool isTrigonalBipyramidal(Stereo::Ref center) const
      {
        return isType(Stereo::TrigonalBipyramidal, center);
      }

      bool isOctahedral(Stereo::Ref center) const
      {
        return isType(Stereo::Octahedral, center);
      }

      const StereoStorage& tetrahedral(Stereo::Ref center) const
      {
        return getType(Stereo::Tetrahedral, center);
      }

      const StereoStorage& allene(Stereo::Ref center) const
      {
        return getType(Stereo::Allene, center);
      }

      const StereoStorage& cisTrans(Stereo::Ref bond) const
      {
        return getType(Stereo::CisTrans, bond);
      }

      const StereoStorage& squarePlanar(Stereo::Ref center) const
      {
        return getType(Stereo::SquarePlanar, center);
      }

      const StereoStorage& trigonalBipyramidal(Stereo::Ref center) const
      {
        return getType(Stereo::TrigonalBipyramidal, center);
      }

      const StereoStorage& octahedral(Stereo::Ref center) const
      {
        return getType(Stereo::Octahedral, center);
      }

      const std::vector<StereoStorage>& allStereo() const
      {
        return m_stereo;
      }

      void clear()
      {
        m_stereo.clear();
      }

      void add(const StereoStorage &stereo)
      {
        m_stereo.push_back(stereo);
      }

      Stereochemistry& operator=(const Stereochemistry &other)
      {
        copy(other);
        return *this;
      }


    private:
      void copy(const Stereochemistry &other)
      {
        m_invalid = other.m_invalid;
        m_stereo = other.m_stereo;
      }

      std::size_t numType(enum Stereo::Type type) const
      {
        std::size_t result = 0;
        for (auto stereo : m_stereo)
          if (stereo.type() == type)
            ++result;
        return result;
      }

      bool isType(enum Stereo::Type type, Stereo::Ref center) const
      {
        for (auto stereo : m_stereo)
          if (stereo.type() == type && stereo.center() == center)
            return true;
        return false;
      }

      const StereoStorage& getType(enum Stereo::Type type, Stereo::Ref center) const
      {
        for (auto &stereo : m_stereo)
          if (stereo.type() == type && stereo.center() == center)
            return stereo;
        return m_invalid;
      }

      StereoStorage m_invalid;
      std::vector<StereoStorage> m_stereo;
  };

  Stereo::Class tetrahedral_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4);

  Stereo::Class squareplanar_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4);

  Stereo::Class trigonalbipyramidal_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4, Stereo::Ref ref5);

  Stereo::Class octahedral_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4, Stereo::Ref ref5, Stereo::Ref ref6);

  template<typename MoleculeType>
  class CisTransHelper
  {
    public:
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      CisTransHelper(const MoleculeType &mol, const StereoStorage &stereo) : m_mol(mol), m_stereo(stereo)
      {
        assert(m_stereo.center() < num_bonds(m_mol));
        m_bond = get_bond(mol, m_stereo.center());
        auto source = get_source(mol, m_bond);
        auto target = get_target(mol, m_bond);

        if (m_stereo.ref(0) != Stereo::implRef()) {
          auto a = refAtom(0);
          if (get_bond(mol, source, a) == molecule_traits<MoleculeType>::null_bond())
            std::swap(source, target);
        } else {
          assert(m_stereo.ref(1) != Stereo::implRef());
          auto b = refAtom(1);
          if (get_bond(mol, source, b) == molecule_traits<MoleculeType>::null_bond())
            std::swap(source, target);
        }

        m_source = source;
        m_target = target;
      }

      atom_type refAtom(int index) const
      {
        assert(index < 4);
        if (m_stereo.ref(index) == Stereo::implRef())
          return molecule_traits<MoleculeType>::null_atom();
        assert(m_stereo.ref(index) < num_atoms(m_mol));
        return get_atom(m_mol, m_stereo.ref(index));
      }

      atom_type source() const
      {
        return m_source;
      }

      atom_type target() const
      {
        return m_target;
      }

      atom_type atomAboveSource() const
      {
        return refAtom(0);
      }

      atom_type atomBelowSource() const
      {
        return refAtom(1);
      }

      atom_type atomAboveTarget() const
      {
        return refAtom(3);
      }

      atom_type atomBelowTarget() const
      {
        return refAtom(2);
      }

    private:
      const MoleculeType &m_mol;
      const StereoStorage &m_stereo;
      bond_type m_bond; // the double bond
      atom_type m_source; // ref nbrs 0 & 1
      atom_type m_target; // ref nbrs 2 & 3
  };

}

#endif
