/*
 * Copyright (c) 2015, Tim Vandermeersch
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
#ifndef HELIUM_STEREO_H
#define HELIUM_STEREO_H

#include <Helium/molecule.h>
#include <Helium/ring.h>

#include <vector>
#include <limits>
#include <cassert>

namespace Helium {

  /**
   * @brief Placeholder for various stereochemistry enums and values.
   */
  struct Stereo
  {
    typedef unsigned int Ref; //!< The type for stereochemistry reference points.

    /**
     * @brief The stereochemistry types.
     */
    enum Type {
      None, //!< No stereochmistry.
      Tetrahedral, //!< Tetrahedral stereochmistry.
      Allene, //!< Allene stereochemistry.
      SquarePlanar, //!< Square-planar stereochemistry.
      TrigonalBipyramidal, //!< Trigonal-bipyramidal stereochmistry.
      Octahedral, //!< Octahedral stereochemistry.
      CisTrans, //!< cis/trans double bond stereochemistry
    };

    enum Class {
      TH1, //!< Tetrahedral: view from 0, (1 2 3) are anti-clockwise.
      TH2, //!< Tetrahedral: view from 0, (1 2 3) are clockwise.
      AL1, //!< Allene: view from 0, (1 2 3) are anti-clockwise.
      AL2, //!< Allene: view from 0, (1 2 3) are clockwise.
      SP1, //!< Square-planar: (1 2 3 4) are layed out in a U shape.
      SP2, //!< Square-planar: (1 2 3 4) are layed out in a Z shape.
      SP3, //!< Square-planar: (1 2 3 4) are layed out in a 4 shape.
      TB1, //!< Trigonal-bipyramidal: view from 0 towards 4: remaining refs are anti-clockwise.
      TB2, //!< Trigonal-bipyramidal: view from 0 towards 4: remaining refs are clockwise.
      TB3, //!< Trigonal-bipyramidal: view from 0 towards 3: remaining refs are anti-clockwise.
      TB4, //!< Trigonal-bipyramidal: view from 0 towards 3: remaining refs are clockwise.
      TB5, //!< Trigonal-bipyramidal: view from 0 towards 2: remaining refs are anti-clockwise.
      TB6, //!< Trigonal-bipyramidal: view from 0 towards 2: remaining refs are clockwise.
      TB7, //!< Trigonal-bipyramidal: view from 0 towards 1: remaining refs are anti-clockwise.
      TB8, //!< Trigonal-bipyramidal: view from 0 towards 1: remaining refs are clockwise.
      TB9, //!< Trigonal-bipyramidal: view from 1 towards 4: remaining refs are anti-clockwise.
      TB10, //!< Trigonal-bipyramidal: view from 1 towards 3: remaining refs are anti-clockwise.
      TB11, //!< Trigonal-bipyramidal: view from 1 towards 4: remaining refs are clockwise.
      TB12, //!< Trigonal-bipyramidal: view from 1 towards 3: remaining refs are clockwise.
      TB13, //!< Trigonal-bipyramidal: view from 1 towards 2: remaining refs are anti-clockwise.
      TB14, //!< Trigonal-bipyramidal: view from 1 towards 2: remaining refs are clockwise.
      TB15, //!< Trigonal-bipyramidal: view from 2 towards 4: remaining refs are anti-clockwise.
      TB16, //!< Trigonal-bipyramidal: view from 2 towards 3: remaining refs are anti-clockwise.
      TB17, //!< Trigonal-bipyramidal: view from 3 towards 4: remaining refs are anti-clockwise.
      TB18, //!< Trigonal-bipyramidal: view from 3 towards 4: remaining refs are clockwise.
      TB19, //!< Trigonal-bipyramidal: view from 2 towards 3: remaining refs are clockwise.
      TB20, //!< Trigonal-bipyramidal: view from 2 towards 4: remaining refs are clockwise.
      OH1,
      OH2,
      OH3,
      OH4,
      OH5,
      OH6,
      OH7,
      OH8,
      OH9,
      OH10,
      OH11,
      OH12,
      OH13,
      OH14,
      OH15,
      OH16,
      OH17,
      OH18,
      OH19,
      OH20,
      OH21,
      OH22,
      OH23,
      OH24,
      OH25,
      OH26,
      OH27,
      OH28,
      OH29,
      OH30
    };

    /**
     * @brief Get the null reference point value.
     *
     * The value returned by this function is not a valid reference point.
     *
     * @return The null reference point value.
     */
    static Ref nullRef()
    {
      return std::numeric_limits<Ref>::max();
    }

    /**
     * @brief Get the implicit reference point value.
     *
     * @return The implicit reference point value.
     */
    static Ref implRef()
    {
      return std::numeric_limits<Ref>::max() - 1;
    }
  };

  /**
   * @brief Class for storing stereochemistry information.
   */
  class StereoStorage
  {
    public:
      /**
       * @brief Constructor.
       *
       * @param type The type of stereochemistry to store.
       */
      StereoStorage(enum Stereo::Type type = Stereo::None) : m_type(type)
      {
        m_center = Stereo::nullRef();
        std::fill(m_refs, m_refs + 6, Stereo::nullRef());
      }

      /**
       * @brief Constructor.
       *
       * @param type The type of stereochemistry to store.
       * @param center The stereogenic atom/bond.
       * @param beginRefs Begin iterator to the reference points.
       * @param endRefs End iterator to the reference points.
       */
      template<typename RefIter>
      StereoStorage(enum Stereo::Type type, Stereo::Ref center,
          RefIter beginRefs, RefIter endRefs) : m_type(type), m_center(center)
      {
        std::copy(beginRefs, endRefs, m_refs);
      }

      /**
       * @brief Check if the stored stereochmistry is valid.
       *
       * This function considers a the stored stereochemistry valid if the
       * following conditions are met:
       * - The type is not Stereo::None
       * - The stereo center is not Stereo::nullRef() or Stereo::implRef()
       * - Tetrahedral: at most 1 Stereo::implRef()
       * - Allene: at most 1 Stereo::implRef() per side
       * - CisTrans: at most 1 Stereo::implRef() per side
       * - Square-planar: at most 2 Stereo::implRef()
       * - Trigonal-Bipyramidal: at most 2 Stereo::implRef()
       * - Octahedral: at most 4 Stereo::implRef()
       *
       * @return True if the stored stereochemistry is valid.
       */
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

      /**
       * @brief Get the stereochemistry type.
       *
       * @return The stereochemistry type.
       */
      enum Stereo::Type type() const
      {
        return m_type;
      }

      /**
       * @brief Get the number of reference points.
       *
       * @return The number of reference points.
       */
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

      /**
       * @brief Get the stereochemistry center.
       *
       * This is the stereogenic atom for most stereochemistry types or the double
       * bond for cis/trans stereochemistry.
       *
       * @return The stereochemistry center.
       */
      Stereo::Ref center() const
      {
        return m_center;
      }

      /**
       * @brief Get a reference point.
       *
       * @param index Index of the reference point to get.
       *
       * @return The reference point with the specified index.
       */
      Stereo::Ref ref(int index) const
      {
        assert(index < numRefs());
        return m_refs[index];
      }

      /**
       * @brief Get a pointer to the reference points.
       *
       * @return A pointer to the reference points.
       */
      const Stereo::Ref* refs() const
      {
        return m_refs;
      }

    private:
      enum Stereo::Type m_type; //!< The stereochemistry type.
      Stereo::Ref m_center; //!< The stereogenic center atom or bond.
      Stereo::Ref m_refs[6]; //!< The reference points.
  };

  /**
   * @brief STL output stream operator for StereoStorage.
   *
   * @param os The STL output stream.
   * @param stereo The StereoStorage object.
   *
   * @return The STL output stream.
   */
  std::ostream& operator<<(std::ostream &os, const StereoStorage &stereo);

  /**
   * @brief Class for storing stereochemistry of a molecule.
   */
  class Stereochemistry
  {
    public:
      /**
       * @brief Get the number of stereogenic units.
       *
       * @return The number of stereogenic units.
       */
      std::size_t numStereo() const
      {
        return m_stereo.size();
      }

      /**
       * @brief Get the number of tetrahedral stereogenic atoms.
       *
       * @return The number of tetrahedral stereogenic atoms.
       */
      std::size_t numTetrahedral() const
      {
        return numType(Stereo::Tetrahedral);
      }

      /**
       * @brief Get the number of allene stereogenic atoms.
       *
       * @return The number of allene stereogenic atoms.
       */
      std::size_t numAllene() const
      {
        return numType(Stereo::Allene);
      }

      /**
       * @brief Get the number of cis/trans stereogenic atoms.
       *
       * @return The number of cis/trans stereogenic atoms.
       */
      std::size_t numCisTrans() const
      {
        return numType(Stereo::CisTrans);
      }

      /**
       * @brief Get the number of square-planar stereogenic atoms.
       *
       * @return The number of square-planar stereogenic atoms.
       */
      std::size_t numSquarePlanar() const
      {
        return numType(Stereo::SquarePlanar);
      }

      /**
       * @brief Get the number of trigonal-bipyramidal stereogenic atoms.
       *
       * @return The number of trigonal-bipyramidal stereogenic atoms.
       */
      std::size_t numTrigonalBipyramidal() const
      {
        return numType(Stereo::TrigonalBipyramidal);
      }

      /**
       * @brief Get the number of octahedral stereogenic atoms.
       *
       * @return The number of octahedral stereogenic atoms.
       */
      std::size_t numOctahedral() const
      {
        return numType(Stereo::Octahedral);
      }

      /**
       * @brief Check if an atom is tetrahedral.
       *
       * @param center The atom index to check.
       *
       * @return True if the atom is a tetrahedral stereo center.
       */
      bool isTetrahedral(Stereo::Ref center) const
      {
        return isType(Stereo::Tetrahedral, center);
      }

      /**
       * @brief Check if an atom is allene.
       *
       * @param center The atom index to check.
       *
       * @return True if the atom is a allene stereo center.
       */
      bool isAllene(Stereo::Ref center) const
      {
        return isType(Stereo::Allene, center);
      }

      /**
       * @brief Check if a bond is cis/trans.
       *
       * @param bond The bond index to check.
       *
       * @return True if the bond is a cis/trans stereo center.
       */
      bool isCisTrans(Stereo::Ref bond) const
      {
        return isType(Stereo::CisTrans, bond);
      }

      /**
       * @brief Check if an atom is square-planar.
       *
       * @param center The atom index to check.
       *
       * @return True if the atom is a square-planar stereo center.
       */
      bool isSquarePlanar(Stereo::Ref center) const
      {
        return isType(Stereo::SquarePlanar, center);
      }

      /**
       * @brief Check if an atom is trigonal-bipyramidal.
       *
       * @param center The atom index to check.
       *
       * @return True if the atom is a trigonal-bipyramidal stereo center.
       */
      bool isTrigonalBipyramidal(Stereo::Ref center) const
      {
        return isType(Stereo::TrigonalBipyramidal, center);
      }

      /**
       * @brief Check if an atom is octahedral.
       *
       * @param center The atom index to check.
       *
       * @return True if the atom is a octahedral stereo center.
       */
      bool isOctahedral(Stereo::Ref center) const
      {
        return isType(Stereo::Octahedral, center);
      }

      /**
       * @brief Get a tetrahedral StereoStorage object.
       *
       * @param center Index of the center atom.
       *
       * @return The tetrahedral StereoStorage object.
       */
      const StereoStorage& tetrahedral(Stereo::Ref center) const
      {
        return getType(Stereo::Tetrahedral, center);
      }

      /**
       * @brief Get a allene StereoStorage object.
       *
       * @param center Index of the center atom.
       *
       * @return The allene StereoStorage object.
       */
      const StereoStorage& allene(Stereo::Ref center) const
      {
        return getType(Stereo::Allene, center);
      }

      /**
       * @brief Get a cis/trans StereoStorage object.
       *
       * @param bond Index of the center bond.
       *
       * @return The cis/trans StereoStorage object.
       */
      const StereoStorage& cisTrans(Stereo::Ref bond) const
      {
        return getType(Stereo::CisTrans, bond);
      }

      /**
       * @brief Get a square-planar StereoStorage object.
       *
       * @param center Index of the center atom.
       *
       * @return The square-planar StereoStorage object.
       */
      const StereoStorage& squarePlanar(Stereo::Ref center) const
      {
        return getType(Stereo::SquarePlanar, center);
      }

      /**
       * @brief Get a trigonal-bipyramidal StereoStorage object.
       *
       * @param center Index of the center atom.
       *
       * @return The trigonal-bipyramidal StereoStorage object.
       */
      const StereoStorage& trigonalBipyramidal(Stereo::Ref center) const
      {
        return getType(Stereo::TrigonalBipyramidal, center);
      }

      /**
       * @brief Get a octahedral StereoStorage object.
       *
       * @param center Index of the center atom.
       *
       * @return The octahedral StereoStorage object.
       */
      const StereoStorage& octahedral(Stereo::Ref center) const
      {
        return getType(Stereo::Octahedral, center);
      }

      /**
       * @brief Get all StereoStorage objects.
       *
       * @return All StereoStorage objects.
       */
      const std::vector<StereoStorage>& allStereo() const
      {
        return m_stereo;
      }

      /**
       * @brief Clear all stored stereochemistry.
       */
      void clear()
      {
        m_stereo.clear();
      }

      /**
       * @brief Add a StereoStorage object.
       *
       * @param stereo The StereoStorage object to add.
       */
      void add(const StereoStorage &stereo)
      {
        m_stereo.push_back(stereo);
      }

    private:
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

      StereoStorage m_invalid; //!< An invalid StereoStorage object.
      std::vector<StereoStorage> m_stereo; //!< The StereoStorage objects.
  };

  /**
   * @brief STL output stream operator for Stereochemistry.
   *
   * @param os The STL output stream.
   * @param stereo The stereochemistry.
   *
   * @return The STL output stream.
   */
  std::ostream& operator<<(std::ostream &os, const Stereochemistry &stereo);

  /**
   * @brief Get the tetrahedral stereochemistry class.
   *
   * This function computes the correct class based on the stereochemistry stored
   * in the StereoStorage object and the order given by ref1 to ref4.
   *
   * @param storage The stored stereochemistry.
   * @param ref1 Reference point 1.
   * @param ref2 Reference point 2.
   * @param ref3 Reference point 3.
   * @param ref4 Reference point 4.
   *
   * @return The tetrahedral stereochemistry class.
   */
  Stereo::Class tetrahedral_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4);

  /**
   * @brief Get the square-planar stereochemistry class.
   *
   * This function computes the correct class based on the stereochemistry stored
   * in the StereoStorage object and the order given by ref1 to ref4.
   *
   * @param storage The stored stereochemistry.
   * @param ref1 Reference point 1.
   * @param ref2 Reference point 2.
   * @param ref3 Reference point 3.
   * @param ref4 Reference point 4.
   *
   * @return The square-planar stereochemistry class.
   */
  Stereo::Class squareplanar_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4);

  /**
   * @brief Get the trigonal-bipyramidal stereochemistry class.
   *
   * This function computes the correct class based on the stereochemistry stored
   * in the StereoStorage object and the order given by ref1 to ref5.
   *
   * @param storage The stored stereochemistry.
   * @param ref1 Reference point 1.
   * @param ref2 Reference point 2.
   * @param ref3 Reference point 3.
   * @param ref4 Reference point 4.
   * @param ref5 Reference point 5.
   *
   * @return The trigonal-bipyramidal stereochemistry class.
   */
  Stereo::Class trigonalbipyramidal_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4, Stereo::Ref ref5);

  /**
   * @brief Get the octahedral stereochemistry class.
   *
   * This function computes the correct class based on the stereochemistry stored
   * in the StereoStorage object and the order given by ref1 to ref6.
   *
   * @param storage The stored stereochemistry.
   * @param ref1 Reference point 1.
   * @param ref2 Reference point 2.
   * @param ref3 Reference point 3.
   * @param ref4 Reference point 4.
   * @param ref5 Reference point 5.
   * @param ref6 Reference point 6.
   *
   * @return The octahedral stereochemistry class.
   */
  Stereo::Class octahedral_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4, Stereo::Ref ref5, Stereo::Ref ref6);

  /**
   * @brief Helper class for working with cis/trans StereoStorage objects.
   */
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
