/*
 * Copyright (c) 2013-2015, Tim Vandermeersch
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
#ifndef HELIUM_SMILES_H
#define HELIUM_SMILES_H

#include <Helium/molecule.h>
#include <Helium/algorithms/dfs.h>
#include <Helium/algorithms/canonical.h>
#include <Helium/algorithms/components.h>
#include <Helium/algorithms/extendedconnectivities.h>
#include <Helium/stereo.h>
#include <Helium/element.h>
#include <Helium/error.h>
#include <Helium/smiley.h>

#include <sstream>


namespace Helium {

  /**
   * @file smiles.h
   * @brief Read/write SMILES.
   */

  /**
   * @class Smiles smiles.h <Helium/smiles.h>
   * @brief Class for reading SMILES.
   */
  class Smiles
  {
    public:
      /**
       * @brief SMILES features to write.
       */
      enum Flags {
        None = 0,
        Mass = 1,
        Charge = 2,
        Hydrogens = 4,
        Order = 8,
        All = Mass | Charge | Hydrogens | Order
      };

      /**
       * @brief Parse a SMILES string.
       *
       * @tparam EditableMoleculeType The type of the molecule, this must be a model
       *         of the EdiatbleMolecule concept.
       *
       * @param smiles The SMILES string to read.
       * @param mol The output molecule.
       * @param stereo The output stereochemistry.
       *
       * @return True if the SMILES was read correctly.
       */
      template<typename EditableMoleculeType>
      bool read(const std::string &smiles, EditableMoleculeType &mol,
          Stereochemistry &stereo);

      /**
       * @overload
       */
      template<typename EditableMoleculeType>
      bool read(const std::string &smiles, EditableMoleculeType &mol)
      {
        Stereochemistry stereo;
        return read(smiles, mol, stereo);
      }

      /**
       * @brief Write a SMILES string for the molecule.
       *
       * @param mol The molecule.
       * @param stereo The stereochemistry.
       * @param atomClasses The atom classes.
       * @param flags The SMILES features to write (see Flags).
       *
       * @return The SMILES string.
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, const Stereochemistry &stereo,
          const std::map<Index, int> &atomClasses, int flags = All);

      /**
       * @overload
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, const std::map<Index, int> &atomClasses,
          int flags = All)
      {
        Stereochemistry stereo;
        return write(mol, stereo, atomClasses, flags);
      }

      /**
       * @overload
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, const Stereochemistry &stereo,
          int flags = All)
      {
        std::map<Index, int> atomClasses;
        return write(mol, stereo, atomClasses, flags);
      }

      /**
       * @overload
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, int flags = All)
      {
        Stereochemistry stereo;
        std::map<Index, int> atomClasses;
        return write(mol, stereo, atomClasses, flags);
      }

      /**
       * @brief Write a SMILES string for the molecule using a specified order.
       *
       * This version of write() can be used to write canonical SMILES if a
       * canonical atom order is used.
       *
       * @param mol The molecule.
       * @param stereo The stereochemistry.
       * @param order The atom order.
       * @param atomClasses The atom classes.
       * @param flags The SMILES features to write (see Flags).
       *
       * @return The SMILES string.
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, const Stereochemistry &stereo,
          const std::vector<Index> &order, const std::map<Index, int> &atomClasses,
          int flags = All);

      /**
       * @overload
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, const std::vector<Index> &order,
          const std::map<Index, int> &atomClasses, int flags = All)
      {
        Stereochemistry stereo;
        return write(mol, stereo, order, atomClasses, flags);
      }

      /**
       * @overload
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, const Stereochemistry &stereo,
          const std::vector<Index> &order, int flags = All)
      {
        std::map<Index, int> atomClasses;
        return write(mol, stereo, order, atomClasses, flags);
      }

      /**
       * @overload
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, const std::vector<Index> &order,
          int flags = All)
      {
        Stereochemistry stereo;
        std::map<Index, int> atomClasses;
        return write(mol, stereo, order, atomClasses, flags);
      }


      /**
       * @brief Write a canonical SMILES string for the molecule.
       *
       * @param mol The molecule.
       * @param stereo The stereochemistry.
       * @param atomClasses The atom classes.
       * @param flags The SMILES features to write (see Flags).
       *
       * @return The SMILES string.
       */
      template<typename MoleculeType>
      std::string writeCanonical(const MoleculeType &mol, const Stereochemistry &stereo,
          const std::map<Index, int> &atomClasses, int flags = All)
      {
        std::pair<std::vector<Index>, std::vector<unsigned long> > canon = canonicalize(mol,
            extended_connectivities(mol, DefaultAtomInvariant(DefaultAtomInvariant::Element)),
            DefaultAtomInvariant(DefaultAtomInvariant::All), DefaultBondInvariant(DefaultBondInvariant::All),
            connected_atom_components(mol), connected_bond_components(mol));
        return write(mol, stereo, canon.first, atomClasses, flags);
      }

      /**
       * @overload
       */
      template<typename MoleculeType>
      std::string writeCanonical(const MoleculeType &mol, const Stereochemistry &stereo,
          int flags = All)
      {
        std::map<Index, int> atomClasses;
        return writeCanonical(mol, stereo, atomClasses, flags);
      }

      /**
       * @overload
       */
      template<typename MoleculeType>
      std::string writeCanonical(const MoleculeType &mol, const std::map<Index, int> &atomClasses,
          int flags = All)
      {
        Stereochemistry stereo;
        return writeCanonical(mol, stereo, atomClasses, flags);
      }

      /**
       * @overload
       */
      template<typename MoleculeType>
      std::string writeCanonical(const MoleculeType &mol, int flags = All)
      {
        Stereochemistry stereo;
        std::map<Index, int> atomClasses;
        return writeCanonical(mol, stereo, atomClasses, flags);
      }

      /**
       * @brief Get the error resulting from calling read().
       *
       * @return The parse error.
       */
      const Error& error() const
      {
        return m_error;
      }

    private:
      Error m_error; //!< The last error.
  };

  namespace impl {

    template<typename EditableMoleculeType>
    struct SmileyCallback : public Smiley::CallbackBase
    {
      SmileyCallback(EditableMoleculeType &mol_, Stereochemistry &stereo_)
          : mol(mol_), stereo(stereo_)
      {
      }

      void clear()
      {
        clear_molecule(mol);
        stereo.clear();
      }

      void addAtom(int element, bool aromatic, int isotope, int hCount, int charge, int atomClass)
      {
        typename molecule_traits<EditableMoleculeType>::atom_type atom = add_atom(mol);
        set_element(mol, atom, element);
        set_aromatic(mol, atom, aromatic);
        if (isotope != -1)
          set_mass(mol, atom, isotope);
        else
          set_mass(mol, atom, Element::averageMass(element));
        if (hCount != -1)
          set_hydrogens(mol, atom, hCount);
        else {
          switch (element) {
            case 5: // B
            case 6: // C
            case 7: // N
            case 8: // O
            case 15: // P
            case 16: // S
            case 9: // F
            case 17: // Cl
            case 35: // Br
            case 53: // I
              set_hydrogens(mol, atom, 99);
              break;
            default:
              break;
          }
        }
        set_charge(mol, atom, charge);
      }

      void addBond(int source, int target, int order, bool isUp, bool isDown)
      {
        //std::cout << "addBond(" << source << ", " << target << ", " << order << ", " << isUp << ", " << isDown << ")" << std::endl;

        auto bond = add_bond(mol, get_atom(mol, source), get_atom(mol, target));
        if (order == 5)
          set_aromatic(mol, bond, true);
        set_order(mol, bond, order);

        if (isUp && !isDown)
          upBonds.insert(get_index(mol, bond));
        if (!isUp && isDown)
          downBonds.insert(get_index(mol, bond));
      }

      void replaceImplicitHydrogens(Stereo::Ref *refs, int size)
      {
        for (std::size_t i = 0; i < size; ++i)
          if (refs[i] == Smiley::implicitHydrogen())
            refs[i] = Stereo::implRef();
      }

      void setChiral(int index, Smiley::Chirality chirality, const std::vector<int> &chiralNbrs)
      {
        //std::cout << "setChiral(index: " << index << ", type: " << chirality << ", chiralNbrs: " << chiralNbrs << ")" << std::endl;

        Stereo::Ref refs[6];
        Stereo::Type type = Stereo::None;
        int numRefs = 0;

        switch (chirality) {
          case Smiley::AntiClockwise:
          case Smiley::Clockwise:
          case Smiley::TH1:
          case Smiley::TH2:
            type = Stereo::Tetrahedral;
            numRefs = 4;
            break;
          case Smiley::AL1:
          case Smiley::AL2:
            type = Stereo::Allene;
            numRefs = 4;
            break;
          case Smiley::SP1:
          case Smiley::SP2:
          case Smiley::SP3:
            type = Stereo::SquarePlanar;
            numRefs = 4;
            break;
          default:
            if (chirality >= Smiley::TB1 && chirality <= Smiley::TB20) {
              type = Stereo::TrigonalBipyramidal;
              numRefs = 5;
            } else if (chirality >= Smiley::OH1 && chirality <= Smiley::OH30) {
              type = Stereo::Octahedral;
              numRefs = 6;
            }
            break;
        }

        assert(type != Stereo::None);
        assert(numRefs);

        if (chiralNbrs.size() != numRefs) {
          // FIXME: report invalid stereochemistry
          return;
        }

        std::copy(chiralNbrs.begin(), chiralNbrs.end(), refs);
        replaceImplicitHydrogens(refs, numRefs);

        switch (chirality) {
          case Smiley::NotChiral:
            break;
          case Smiley::AntiClockwise:
          case Smiley::TH1:
            break;
          case Smiley::Clockwise:
          case Smiley::TH2:
            type = Stereo::Tetrahedral;
            refs[0] = chiralNbrs[0];
            std::copy(chiralNbrs.rbegin(), chiralNbrs.rbegin() + 3, refs + 1);
            replaceImplicitHydrogens(refs, 4);
            break;
          case Smiley::AL1:
            break;
          case Smiley::AL2:
            type = Stereo::Allene;
            refs[0] = chiralNbrs[0];
            std::copy(chiralNbrs.rbegin(), chiralNbrs.rbegin() + 3, refs + 1);
            replaceImplicitHydrogens(refs, 4);
            break;
          case Smiley::SP1:
            // U shape
            break;
          case Smiley::SP2:
            // 4 shape -> U shape
            std::swap(refs[1], refs[2]);
            break;
          case Smiley::SP3:
            // Z shape -> U shape
            std::swap(refs[0], refs[1]);
            break;
          case Smiley::TB1: // axis a-e, other @
            break;
          case Smiley::TB2: // axis a-e, other @@
            std::swap(refs[0], refs[4]); // abcde -> ebcda
            break;
          case Smiley::TB3: // axis a-d, other @
            std::swap(refs[3], refs[4]); // abcde -> abced
            break;
          case Smiley::TB4: // axis a-d, other @@
            std::swap(refs[3], refs[4]); // abcde -> abced
            std::swap(refs[0], refs[4]); // abcde -> dbcea
            break;
          case Smiley::TB5: // axis a-c, other @
            std::rotate(refs + 2, refs + 3, refs + 5); // abcde -> abdec
            break;
          case Smiley::TB6: // axis a-c, other @@
            std::rotate(refs + 2, refs + 3, refs + 5); // abcde -> abdec
            std::swap(refs[0], refs[4]); // abdec -> cbdea
            break;
          case Smiley::TB7: // axis a-b, other @
            std::rotate(refs + 1, refs + 2, refs + 5); // abcde -> acdeb
            break;
          case Smiley::TB8: // axis a-b, other @@
            std::rotate(refs + 1, refs + 2, refs + 5); // abcde -> acdeb
            std::swap(refs[0], refs[4]); // abdec -> cbdea
            break;
          case Smiley::TB9: // axis b-e, other @
            std::swap(refs[0], refs[1]); // abcde -> bacde
            break;
          case Smiley::TB11: // axis b-e, other @@
            std::swap(refs[0], refs[1]); // abcde -> bacde
            std::swap(refs[0], refs[4]); // bacde -> eacdb
            break;
          case Smiley::TB10: // axis b-d, other @
            std::swap(refs[0], refs[1]); // abcde -> bacde
            std::swap(refs[3], refs[4]); // bacde -> baced
            break;
          case Smiley::TB12: // axis b-d, other @@
            std::swap(refs[0], refs[1]); // abcde -> bacde
            std::swap(refs[3], refs[4]); // bacde -> baced
            std::swap(refs[0], refs[4]); // bacde -> eacdb
            break;
          case Smiley::TB13: // axis b-c, other @
            std::swap(refs[0], refs[1]); // abcde -> bacde
            std::rotate(refs + 2, refs + 3, refs + 5); // bacde -> badec
            break;
          case Smiley::TB14: // axis b-c, other @@
            std::swap(refs[0], refs[1]); // abcde -> bacde
            std::rotate(refs + 2, refs + 3, refs + 5); // bacde -> badec
            std::swap(refs[0], refs[4]); // badec -> cadeb
            break;
          case Smiley::TB15: // axis c-e, other @
            std::swap(refs[1], refs[2]); // abcde -> acbde
            std::swap(refs[0], refs[1]); // acbde -> cabde
            break;
          case Smiley::TB20: // axis c-e, other @@
            std::swap(refs[1], refs[2]); // abcde -> acbde
            std::swap(refs[0], refs[1]); // acbde -> cabde
            std::swap(refs[0], refs[4]); // cabde -> eabdc
            break;
          case Smiley::TB16: // axis c-d, other @
            std::swap(refs[1], refs[2]); // abcde -> acbde
            std::swap(refs[0], refs[1]); // acbde -> cabde
            std::swap(refs[3], refs[4]); // cabde -> cabed
            break;
          case Smiley::TB19: // axis c-d, other @@
            std::swap(refs[1], refs[2]); // abcde -> acbde
            std::swap(refs[0], refs[1]); // acbde -> cabde
            std::swap(refs[3], refs[4]); // cabde -> cabed
            std::swap(refs[0], refs[4]); // cabde -> eabdc
            break;
          case Smiley::TB17: // axis d-e, other @
            std::swap(refs[2], refs[3]); // abcde -> abdce
            std::swap(refs[1], refs[2]); // abdce -> adbce
            std::swap(refs[0], refs[1]); // adbce -> dabce
            break;
          case Smiley::TB18: // axis d-e, other @@
            std::swap(refs[2], refs[3]); // abcde -> abdce
            std::swap(refs[1], refs[2]); // abdce -> adbce
            std::swap(refs[0], refs[1]); // adbce -> dabce
            std::swap(refs[0], refs[4]); // cabde -> eabdc
            break;
          case Smiley::OH1: // axis a-f, other U-shape @
            break;
          case Smiley::OH2: // axis a-f, other U-shape @@
            std::reverse(refs + 1, refs + 5);
            break;
          case Smiley::OH3: // axis a-e, other U-shape @
            std::swap(refs[4], refs[5]); // abcdef -> abcdfe (a-f -> a-e)
            break;
          case Smiley::OH16: // axis a-e, other U-shape @@
            std::swap(refs[4], refs[5]); // abcdef -> abcdfe (a-f -> a-e)
            std::swap(refs[0], refs[5]); // abcdfe -> ebcdfa
            break;
          case Smiley::OH6: // axis a-d, other U-shape @
            std::rotate(refs + 3, refs + 4, refs + 6); // abcdef -> abcefd (a-f -> a-d)
            break;
          case Smiley::OH18: // axis a-d, other U-shape @@
            std::rotate(refs + 3, refs + 4, refs + 6); // abcdef -> abcefd (a-f -> a-d)
            std::swap(refs[0], refs[5]); // abcefd -> dbcefa
            break;
          case Smiley::OH19: // axis a-c, other U-shape @
            std::rotate(refs + 2, refs + 3, refs + 6); // abcdef -> abdefc (a-f -> a-c)
            break;
          case Smiley::OH24: // axis a-c, other U-shape @@
            std::rotate(refs + 2, refs + 3, refs + 6); // abcdef -> abdefc (a-f -> a-c)
            std::swap(refs[0], refs[5]); // abdefc -> cbdefa
            break;
          case Smiley::OH25: // axis a-b, other U-shape @
            std::rotate(refs + 1, refs + 2, refs + 6); // abcdef -> acdefb (a-f -> a-b)
            break;
          case Smiley::OH30: // axis a-b, other U-shape @@
            std::rotate(refs + 1, refs + 2, refs + 6); // abcdef -> acdefb (a-f -> a-b)
            std::swap(refs[0], refs[5]); // acdefb -> bcdefa
            break;
          case Smiley::OH4: // axis a-f, other Z-shape @
            std::swap(refs[3], refs[4]); // Z -> U
            break;
          case Smiley::OH14: // axis a-f, other Z-shape @@
            std::swap(refs[3], refs[4]); // Z -> U
            std::swap(refs[0], refs[5]); // a-f -> f-a
            break;
          case Smiley::OH5: // axis a-e, other Z-shape @
            std::swap(refs[4], refs[5]); // abcdef -> abcdfe (a-f -> a-e)
            std::swap(refs[3], refs[4]); // Z -> U
            break;
          case Smiley::OH15: // axis a-e, other Z-shape @@
            std::swap(refs[4], refs[5]); // abcdef -> abcdfe (a-f -> a-e)
            std::swap(refs[3], refs[4]); // Z -> U
            std::swap(refs[0], refs[5]); // a-f -> f-a
            break;
          case Smiley::OH7: // axis a-d, other Z-shape @
            std::rotate(refs + 3, refs + 4, refs + 6); // abcdef -> abcefd (a-f -> a-d)
            std::swap(refs[3], refs[4]); // Z -> U
            break;
          case Smiley::OH17: // axis a-d, other Z-shape @@
            std::rotate(refs + 3, refs + 4, refs + 6); // abcdef -> abcefd (a-f -> a-d)
            std::swap(refs[3], refs[4]); // Z -> U
            std::swap(refs[0], refs[5]); // a-f -> f-a
            break;
          case Smiley::OH20: // axis a-c, other Z-shape @
            std::rotate(refs + 2, refs + 3, refs + 6); // abcdef -> abdefc (a-f -> a-c)
            std::swap(refs[3], refs[4]); // Z -> U
            break;
          case Smiley::OH23: // axis a-c, other Z-shape @@
            std::rotate(refs + 2, refs + 3, refs + 6); // abcdef -> abdefc (a-f -> a-c)
            std::swap(refs[3], refs[4]); // Z -> U
            std::swap(refs[0], refs[5]); // a-f -> f-a
            break;
          case Smiley::OH26: // axis a-b, other Z-shape @
            std::rotate(refs + 1, refs + 2, refs + 6); // abcdef -> acdefb (a-f -> a-b)
            std::swap(refs[3], refs[4]); // Z -> U
            break;
          case Smiley::OH29: // axis a-b, other Z-shape @@
            std::rotate(refs + 1, refs + 2, refs + 6); // abcdef -> acdefb (a-f -> a-b)
            std::swap(refs[3], refs[4]); // Z -> U
            std::swap(refs[0], refs[5]); // a-f -> f-a
            break;
          case Smiley::OH10: // axis a-f, other 4-shape @
            std::swap(refs[1], refs[4]); // 4 -> U
            break;
          case Smiley::OH8: // axis a-f, other 4-shape @@
            std::swap(refs[1], refs[4]); // 4 -> U
            std::swap(refs[0], refs[5]); // a-f -> f-a
            break;
          case Smiley::OH11: // axis a-e, other 4-shape @
            std::swap(refs[4], refs[5]); // abcdef -> abcdfe (a-f -> a-e)
            std::swap(refs[1], refs[4]); // 4 -> U
            break;
          case Smiley::OH9: // axis a-e, other 4-shape @@
            std::swap(refs[4], refs[5]); // abcdef -> abcdfe (a-f -> a-e)
            std::swap(refs[1], refs[4]); // 4 -> U
            std::swap(refs[0], refs[5]); // a-f -> f-a
            break;
          case Smiley::OH13: // axis a-d, other 4-shape @
            std::rotate(refs + 3, refs + 4, refs + 6); // abcdef -> abcefd (a-f -> a-d)
            std::swap(refs[1], refs[4]); // 4 -> U
            break;
          case Smiley::OH12: // axis a-d, other 4-shape @@
            std::rotate(refs + 3, refs + 4, refs + 6); // abcdef -> abcefd (a-f -> a-d)
            std::swap(refs[1], refs[4]); // 4 -> U
            std::swap(refs[0], refs[5]); // a-f -> f-a
            break;
          case Smiley::OH22: // axis a-c, other 4-shape @
            std::rotate(refs + 2, refs + 3, refs + 6); // abcdef -> abdefc (a-f -> a-c)
            std::swap(refs[1], refs[4]); // 4 -> U
            break;
          case Smiley::OH21: // axis a-c, other 4-shape @@
            std::rotate(refs + 2, refs + 3, refs + 6); // abcdef -> abdefc (a-f -> a-c)
            std::swap(refs[1], refs[4]); // 4 -> U
            std::swap(refs[0], refs[5]); // a-f -> f-a
            break;
          case Smiley::OH28: // axis a-b, other 4-shape @
            std::rotate(refs + 1, refs + 2, refs + 6); // abcdef -> acdefb (a-f -> a-b)
            std::swap(refs[1], refs[4]); // 4 -> U
            break;
          case Smiley::OH27: // axis a-b, other 4-shape @@
            std::rotate(refs + 1, refs + 2, refs + 6); // abcdef -> acdefb (a-f -> a-b)
            std::swap(refs[1], refs[4]); // 4 -> U
            std::swap(refs[0], refs[5]); // a-f -> f-a
            break;
        }

        stereo.add(StereoStorage(type, index, refs, refs + numRefs));
      }

      std::set<Index> upBonds;
      std::set<Index> downBonds;
      EditableMoleculeType &mol;
      Stereochemistry &stereo;
    };

    struct RingNumber
    {
      int number; // the ring number determined by the preprocessor visitor (may be >99)
      int atom1; // index of atom that is written first
      int atom2; // index of atom that is written second
      int order; // the bond order
    };

    // determines ring closure numbers and final atom output order (needed for stereochemistry)
    template<typename MoleculeType>
    struct WriteSmilesPreprocessor : public DFSVisitor<MoleculeType>
    {
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      void initialize(const MoleculeType &mol)
      {
        outputIndex = 0;
        ringNumber = 0;
        ringNumbers.clear();
        outputIndices.resize(num_atoms(mol));
      }

      void atom(const MoleculeType &mol, atom_type prev, atom_type atom)
      {
        outputIndices[get_index(mol, atom)] = outputIndex++;
      }

      void back_bond(const MoleculeType &mol, bond_type bond)
      {
        ++ringNumber;

        int source = get_index(mol, get_source(mol, bond));
        int target = get_index(mol, get_target(mol, bond));

        // both source and target are already visited, use outputIndices to
        // swap them if needed
        if (outputIndices[target] < outputIndices[source])
          std::swap(source, target);

        int order = 0;
        if (get_order(mol, bond) == 1 && !is_aromatic(mol, bond) &&
            is_aromatic(mol, get_source(mol, bond)) && is_aromatic(mol, get_target(mol, bond)))
          order = 1;
        else if (get_order(mol, bond) == 2 && !is_aromatic(mol, bond))
          order = 2;
        else if (get_order(mol, bond) > 2)
          order = get_order(mol, bond);

        RingNumber rn{ringNumber, source, target, order};

        ringNumbers[get_source(mol, bond)].push_back(rn);
        ringNumbers[get_target(mol, bond)].push_back(rn);
      }

      int outputIndex;
      int ringNumber;
      std::map<atom_type, std::vector<RingNumber>> ringNumbers;
      std::vector<int> outputIndices; // maps atom index to output index
    };


    template<typename MoleculeType>
    struct WriteSmilesVisitor : public DFSVisitor<MoleculeType>
    {
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      WriteSmilesVisitor(const Stereochemistry &stereo_,
          const std::map<atom_type, std::vector<RingNumber>> &ringNumbers_,
          const std::vector<int> &outputIndices_, const std::map<Index, int> &atomClasses_,
          int flags_) : stereo(stereo_), ringNumbers(ringNumbers_),
          outputIndices(outputIndices_), atomClasses(atomClasses_),
          explicitBond(0), flags(flags_)
      {
      }

      void initialize(const MoleculeType &mol)
      {
        degrees.resize(num_atoms(mol));
        explicitBond = 0;
        determineUpDownBonds(mol);
      }

      bool isOrganicSubset(int element) const
      {
        switch (element) {
          case 5: // B
          case 6: // C
          case 7: // N
          case 8: // O
          case 9: // F
          case 15: // P
          case 16: // S
          case 17: // Cl
          case 35: // Br
          case 53: // I
            return true;
          default:
            return false;
        }
      }

      void writeStereo(const MoleculeType &mol, atom_type prev, atom_type atom,
          const std::vector<RingNumber> &rings)
      {
        Stereo::Type type = Stereo::None;
        int numRefs = 0;
        if (stereo.isTetrahedral(get_index(mol, atom))) {
          numRefs = 4;
          type = Stereo::Tetrahedral;
        } else if (stereo.isAllene(get_index(mol, atom))) {
          numRefs = 4;
          type = Stereo::Allene;
        } else if (stereo.isSquarePlanar(get_index(mol, atom))) {
          numRefs = 4;
          type = Stereo::SquarePlanar;
        } else if (stereo.isTrigonalBipyramidal(get_index(mol, atom))) {
          numRefs = 5;
          type = Stereo::TrigonalBipyramidal;
        } else if (stereo.isOctahedral(get_index(mol, atom))) {
          numRefs = 6;
          type = Stereo::Octahedral;
        }

        if (type == Stereo::None)
          return;

        std::map<Stereo::Ref, Stereo::Ref> indexMap; // output index -> molecule index
        // check if there is a neighbor before the chiral atom and sort the other neighbors
        // that are not ring closures by output index
        for (auto &nbr : get_nbrs(mol, atom)) {
          if (nbr == prev)
            continue;
          if (ringNumbers.find(nbr) != ringNumbers.end())
            continue;
          indexMap[outputIndices[get_index(mol, nbr)]] = get_index(mol, nbr);
        }

        std::vector<Stereo::Ref> refs;
        // add prev atom (if needed)
        if (prev != molecule_traits<MoleculeType>::null_atom())
          refs.push_back(get_index(mol, prev));
        // add the ring closure neighbors
        for (auto &ring : rings)
          refs.push_back((get_index(mol, atom) == ring.atom1) ? ring.atom2 : ring.atom1);
        // add sorted neighbors to refs
        for (auto i : indexMap)
          refs.push_back(i.second);

        while (refs.size() < numRefs) {
          if (prev != molecule_traits<MoleculeType>::null_atom())
            refs.insert(refs.begin() + 1, Stereo::implRef());
          else
            refs.insert(refs.begin(), Stereo::implRef());
        }

        assert(refs.size() == numRefs);

        switch (type) {
          case Stereo::Tetrahedral:
            if (tetrahedral_class(stereo.tetrahedral(get_index(mol, atom)), refs[0], refs[1], refs[2], refs[3]) == Stereo::TH1)
              smiles << "@";
            else
              smiles << "@@";
            break;
          case Stereo::Allene:
            if (tetrahedral_class(stereo.allene(get_index(mol, atom)), refs[0], refs[1], refs[2], refs[3]) == Stereo::TH1)
              smiles << "@";
            else
              smiles << "@@";
            break;
          case Stereo::SquarePlanar:
            switch (squareplanar_class(stereo.squarePlanar(get_index(mol, atom)), refs[0], refs[1], refs[2], refs[3])) {
              case Stereo::SP1:
                smiles << "@SP1";
                break;
              case Stereo::SP2:
                smiles << "@SP2";
                break;
              case Stereo::SP3:
                smiles << "@SP3";
                break;
              default:
                break;
            }
            break;
          case Stereo::TrigonalBipyramidal:
            smiles << "@TB" << trigonalbipyramidal_class(stereo.trigonalBipyramidal(get_index(mol, atom)),
                refs[0], refs[1], refs[2], refs[3], refs[4]) - Stereo::TB1 + 1;
            break;
          case Stereo::Octahedral:
            smiles << "@OH" << octahedral_class(stereo.octahedral(get_index(mol, atom)),
                refs[0], refs[1], refs[2], refs[3], refs[4], refs[5]) - Stereo::OH1 + 1;
            break;
          default:
            break;
        }
      }

      void atom(const MoleculeType &mol, atom_type prev, atom_type atom)
      {
        if (prev != molecule_traits<MoleculeType>::null_atom()) {
          auto rings = ringNumbers.find(prev);
          int numRings = (rings == ringNumbers.end()) ? 0 : rings->second.size();
          degrees[get_index(mol, prev)]++;

          // the -1 comes from the previous atom of prev
          // (this is compensated below if prev is a first atom in it's component)
          if (degrees[get_index(mol, prev)] < get_degree(mol, prev) - 1 - numRings) {
            smiles << "(";
            branches.push_back(get_index(mol, atom));
          }
        } else {
          degrees[get_index(mol, atom)] = -1; // compensate for -1 above

          if (smiles.str().size())
            smiles << ".";
        }

        // write explicit bond if needed
        if (explicitBond)
          smiles << explicitBond;
        explicitBond = 0;

        // convert symbol to lower case for aromatic atoms
        std::string element = Element::symbol(get_element(mol, atom));
        if (is_aromatic(mol, atom))
          std::transform(element.begin(), element.end(), element.begin(), ::tolower);

        // non organic subset atoms always requre brackets
        bool needBrackets = !isOrganicSubset(get_element(mol, atom));
        // charge requires brackets
        if (get_charge(mol, atom) && (flags & Smiles::Charge))
          needBrackets = true;
        // ignore mass 0
        if (get_mass(mol, atom) && get_mass(mol, atom) !=
            Element::averageMass(get_element(mol, atom)) && (flags & Smiles::Mass))
          needBrackets = true;
        if (flags & Smiles::Hydrogens) {
          int charge = (flags & Smiles::Charge) ? get_charge(mol, atom) : 0;
          //std::cout << get_valence(mol, atom) << " != " << Element::valence(get_element(mol, atom), charge, get_degree(mol, atom)) << std::endl;
          //std::cout << "#H: " << get_hydrogens(mol, atom) << std::endl;
          if (get_valence(mol, atom) < Element::valence(get_element(mol, atom), charge, get_valence(mol, atom)))
            needBrackets = true;
          // special case: halogens with implicit hydrogens and valence >= 3
          if (get_hydrogens(mol, atom) && get_valence(mol, atom) >= 3) {
            switch (get_element(mol, atom)) {
              case 9:
              case 17:
              case 35:
              case 53:
                needBrackets = true;
              default:
                break;
            }
          }

          // handle pyrrole [nH]
          if (is_aromatic(mol, atom) && is_nitrogen(mol, atom) && get_total_hydrogens(mol, atom) == 1)
            needBrackets = true;
        }
        std::map<Index, int>::const_iterator atomClass = atomClasses.find(get_index(mol, atom));
        if (atomClass != atomClasses.end())
          needBrackets = true;
        // stereo requires brackets
        if (stereo.isTetrahedral(get_index(mol, atom)) ||
            stereo.isAllene(get_index(mol, atom)) ||
            stereo.isSquarePlanar(get_index(mol, atom)) ||
            stereo.isTrigonalBipyramidal(get_index(mol, atom)) ||
            stereo.isOctahedral(get_index(mol, atom)))
          needBrackets = true;

        if (needBrackets)
          smiles << "[";

        if (get_mass(mol, atom) && get_mass(mol, atom) != Element::averageMass(get_element(mol, atom)) && (flags & Smiles::Mass))
          smiles << get_mass(mol, atom);

        smiles << element;

        auto rings = ringNumbers.find(atom);

        // write chiral
        if (rings != ringNumbers.end())
          writeStereo(mol, prev, atom, rings->second);
        else
          writeStereo(mol, prev, atom, std::vector<RingNumber>());

        // do not write implicit H for hydrogen
        // place it as an explicit hydrogen next...
        if (get_element(mol, atom) != 1) {
          int numH = get_hydrogens(mol, atom);
          if (needBrackets && (flags & Smiles::Hydrogens) && numH) {
            smiles << "H";
            if (numH > 1)
              smiles << numH;
          }
        }

        if ((flags & Smiles::Charge) && get_charge(mol, atom)) {
          int charge = get_charge(mol, atom);
          if (charge == 1)
            smiles << "+";
          else if (charge == -1)
            smiles << "-";
          else if (charge > 0)
            smiles << "+" << charge;
          else
            smiles << charge;
        }

        if (atomClass != atomClasses.end())
          smiles << ":" << atomClass->second;

        if (needBrackets)
          smiles << "]";

        // handle molecular hydrogen
        if (get_element(mol, atom) == 1 && get_hydrogens(mol, atom) == 1)
          smiles << "[H]";

        if (rings != ringNumbers.end()) {
          for (auto &ring : rings->second) {
            bool firstAtom = (get_index(mol, atom) == ring.atom1);
            if (flags & Smiles::Order && firstAtom)
              switch (ring.order) {
                case 1:
                  smiles << "-";
                  break;
                case 2:
                  smiles << "=";
                  break;
                case 3:
                  smiles << "#";
                  break;
                case 4:
                  smiles << "$";
                  break;
                default:
                  break;
              }

            if (ring.number > 9)
              smiles << "%" << ring.number;
            else
              smiles << ring.number;
          }
        }
      }

      void bond(const MoleculeType &mol, atom_type prev, bond_type bond)
      {
        if (!(flags & Smiles::Order))
          return;

        explicitBond = 0;

        switch (get_order(mol, bond)) {
          case 1:
            if (!is_aromatic(mol, bond) && is_aromatic(mol, get_source(mol, bond)) && is_aromatic(mol, get_target(mol, bond)))
              explicitBond = '-';
            if (upBonds[get_index(mol, bond)])
              explicitBond = '/';
            else if (downBonds[get_index(mol, bond)])
              explicitBond = '\\';
            break;
          case 2:
            if (!is_aromatic(mol, bond))
              explicitBond = '=';
            break;
          case 3:
            explicitBond = '#';
            break;
          case 4:
            explicitBond = '$';
            break;
          default:
            break;
        }
      }

      void backtrack(const MoleculeType &mol, atom_type atom)
      {
        if (branches.size() && branches.back() == get_index(mol, atom)) {
          smiles << ")";
          branches.pop_back();
        }
      }

      void back_bond(const MoleculeType &mol, bond_type bond)
      {
      }


      bool determineUpDownBonds(const MoleculeType &mol, const std::vector<const StereoStorage*> &cistrans)
      {
        std::fill(upBonds.begin(), upBonds.end(), false);
        std::fill(downBonds.begin(), downBonds.end(), false);

        for (auto bs : cistrans) {
          //  a           d
          //   \         /
          //    c1 === c2
          //   /         \
          //  b           c
          CisTransHelper<MoleculeType> ct(mol, *bs);

          //std::cout << *bs << std::endl;

          auto c1 = ct.source();
          auto c2 = ct.target();
          auto a = ct.atomAboveSource();
          auto b = ct.atomBelowSource();
          auto c = ct.atomAboveTarget();
          auto d = ct.atomBelowTarget();

          /*
          std::cout << "c1: " << c1 << std::endl;
          std::cout << "c2: " << c2 << std::endl;
          std::cout << "a: " << a << std::endl;
          std::cout << "b: " << b << std::endl;
          std::cout << "c: " << c << std::endl;
          std::cout << "d: " << d << std::endl;
          */

          bool aIsImpl = (a == molecule_traits<MoleculeType>::null_atom());
          bool bIsImpl = (b == molecule_traits<MoleculeType>::null_atom());
          bool cIsImpl = (c == molecule_traits<MoleculeType>::null_atom());
          bool dIsImpl = (d == molecule_traits<MoleculeType>::null_atom());

          /*
          std::cout << "aIsImpl: " << aIsImpl << std::endl;
          std::cout << "bIsImpl: " << bIsImpl << std::endl;
          std::cout << "cIsImpl: " << cIsImpl << std::endl;
          std::cout << "dIsImpl: " << dIsImpl << std::endl;
          */

          auto c_a = aIsImpl ? molecule_traits<MoleculeType>::null_bond() : get_bond(mol, c1, a);
          auto c_b = bIsImpl ? molecule_traits<MoleculeType>::null_bond() : get_bond(mol, c1, b);
          auto c_c = cIsImpl ? molecule_traits<MoleculeType>::null_bond() : get_bond(mol, c2, c);
          auto c_d = dIsImpl ? molecule_traits<MoleculeType>::null_bond() : get_bond(mol, c2, d);

          /*
          std::cout << "c_a: " << c_a << std::endl;
          std::cout << "c_b: " << c_b << std::endl;
          std::cout << "c_c: " << c_c << std::endl;
          std::cout << "c_d: " << c_d << std::endl;
          */


          bool aIsUp = aIsImpl ? false : upBonds[get_index(mol, c_a)];
          bool bIsUp = bIsImpl ? false : upBonds[get_index(mol, c_b)];
          bool cIsUp = cIsImpl ? false : upBonds[get_index(mol, c_c)];
          bool dIsUp = dIsImpl ? false : upBonds[get_index(mol, c_d)];

          bool aIsDown = aIsImpl ? false : downBonds[get_index(mol, c_a)];
          bool bIsDown = bIsImpl ? false : downBonds[get_index(mol, c_b)];
          bool cIsDown = cIsImpl ? false : downBonds[get_index(mol, c_c)];
          bool dIsDown = dIsImpl ? false : downBonds[get_index(mol, c_d)];



          // check conflicts:
          // a-b
          if ((aIsUp & bIsUp) | (aIsDown & bIsDown))
            return false;
          // a-c
          if ((aIsUp & cIsUp) | (aIsDown & cIsDown))
            return false;
          // a-d
          if ((aIsUp & dIsDown) | (dIsUp & aIsDown))
            return false;
          // b-c
          if ((bIsUp & cIsDown) | (cIsUp & bIsDown))
            return false;
          // b-d
          if ((bIsUp & dIsUp) | (bIsDown & dIsDown))
            return false;
          // c-d
          if ((cIsUp & dIsUp) | (cIsDown & dIsDown))
            return false;

          bool invert = aIsDown | bIsUp | cIsUp | dIsDown;

          if (!(aIsUp | aIsDown | bIsUp | cIsUp)) {
            if (bIsImpl) {
              assert(!aIsImpl);
              assert(c_a != molecule_traits<MoleculeType>::null_bond());
              if (invert)
                downBonds[get_index(mol, c_a)] = true;
              else
                upBonds[get_index(mol, c_a)] = true;
            } else {
              assert(c_b != molecule_traits<MoleculeType>::null_bond());
              if (invert)
                upBonds[get_index(mol, c_b)] = true;
              else
                downBonds[get_index(mol, c_b)] = true;
            }
          }
          if (!(cIsUp | cIsDown | dIsUp | dIsUp)) {
            if (dIsImpl) {
              if (invert)
                upBonds[get_index(mol, c_c)] = true;
              else
                downBonds[get_index(mol, c_c)] = true;
            } else {
              if (invert)
                downBonds[get_index(mol, c_d)] = true;
              else
                upBonds[get_index(mol, c_d)] = true;
            }
          }

        }

        return true;
      }

      void determineUpDownBonds(const MoleculeType &mol)
      {
        upBonds.resize(num_bonds(mol));
        downBonds.resize(num_bonds(mol));

        std::vector<const StereoStorage*> cistrans;
        for (std::size_t i = 0; i < stereo.allStereo().size(); ++i) {
          auto &ct = stereo.allStereo()[i];
          if (ct.type() == Stereo::CisTrans)
            cistrans.push_back(&ct);
        }
        std::sort(cistrans.begin(), cistrans.end());

        bool success = false;
        do {
          success = determineUpDownBonds(mol, cistrans);
          if (success)
            break;
        } while (std::next_permutation(cistrans.begin(), cistrans.end()));

        assert(success);
      }


      const Stereochemistry &stereo;
      std::vector<bool> upBonds;
      std::vector<bool> downBonds;
      const std::map<atom_type, std::vector<RingNumber> > &ringNumbers; // atom -> ring numbers
      const std::vector<int> &outputIndices;
      std::vector<int> degrees;
      std::vector<Index> branches;
      std::stringstream smiles;
      const std::map<Index, int> &atomClasses;
      char explicitBond;
      int flags;
    };

    inline int smiles_valence(int element, int valence)
    {
      switch (element) {
        case 5: // B
          return 3;
        case 6: // C
          return 4;
        case 7: // N
          if (valence <= 3)
            return 3;
          if (valence <= 5)
            return 5;
          return 0;
        case 8: // O
          return 2;
        case 15: // P
          if (valence <= 3)
            return 3;
          if (valence <= 5)
            return 5;
          return 0;
        case 16: // S
          if (valence <= 2)
            return 2;
          if (valence <= 4)
            return 4;
          if (valence <= 6)
            return 6;
          return 0;
        case 9: // F
        case 17: // Cl
        case 35: // Br
        case 53: // I
          if (valence <= 1)
            return 1;
          return 0;
        default:
          return 0;
      }
    }

    // return false if there is no cis/trans stereochemistry
    template<typename EditableMoleculeType, typename AtomType, typename BondType>
    bool processUpDownBonds(const EditableMoleculeType &mol, const AtomType &doubleBondAtom,
        const std::set<Index> &upBonds, const std::set<Index> &downBonds,
        const BondType &bond1, const BondType &bond2, Stereo::Ref &upRef, Stereo::Ref &downRef)
    {
      auto other1 = get_other(mol, bond1, doubleBondAtom);

      //std::cout << "processUpDownBonds(atom: " << get_index(mol, doubleBondAtom) << ")" << std::endl;

      // check if bond is up/down
      bool isRingBond1 = get_source(mol, bond1) > get_target(mol, bond1);
      bool isUp1 = upBonds.find(get_index(mol, bond1)) != upBonds.end();
      bool isDown1 = downBonds.find(get_index(mol, bond1)) != downBonds.end();
      // if other atom is before atom1, meaning of / and \ changes
      if (get_index(mol, other1) < get_index(mol, doubleBondAtom) && !isRingBond1)
        std::swap(isUp1, isDown1);

      // if there is only 1 bond that does not have up/down set, the double bond is not stereogenic
      if (!isUp1 && !isDown1 && bond2 == molecule_traits<EditableMoleculeType>::null_bond())
        return false;
      if (isUp1)
        upRef = get_index(mol, other1);
      if (isDown1)
        downRef = get_index(mol, other1);

      if (bond2 != molecule_traits<EditableMoleculeType>::null_bond()) {
        auto other2 = get_other(mol, bond2, doubleBondAtom);

        bool isRingBond2 = get_source(mol, bond2) > get_target(mol, bond2);
        bool isUp2 = upBonds.find(get_index(mol, bond2)) != upBonds.end();
        bool isDown2 = downBonds.find(get_index(mol, bond2)) != downBonds.end();
        if (get_index(mol, other2) < get_index(mol, doubleBondAtom) && !isRingBond2)
          std::swap(isUp2, isDown2);

        // if bond1 and bond2 are not up/down, there is no stereochemistry
        if (!isUp1 && !isDown1 && !isUp2 && !isDown2)
          return false;

        if (!isUp1 && !isDown1) {
          // case where only bond 2 specifies up/down
          if (isUp2) {
            upRef = get_index(mol, other2);
            downRef = get_index(mol, other1);
          }
          if (isDown2) {
            downRef = get_index(mol, other2);
            upRef = get_index(mol, other1);
          }
        } else if (!isUp2 && !isDown2) {
          // case where only bond 2 specifies up/down
          if (isUp1) {
            upRef = get_index(mol, other1);
            downRef = get_index(mol, other2);
          }
          if (isDown1) {
            downRef = get_index(mol, other1);
            upRef = get_index(mol, other2);
          }
        } else {
          // check for conflict
          if ((isUp1 && isUp2) || (isDown1 && isDown2))
            return false;

          if (isUp2)
            upRef = get_index(mol, other2);
          if (isDown2)
            downRef = get_index(mol, other2);
        }

        if (!isUp1 && !isDown1 && !isUp2 && !isDown2)
          return false;

        assert(upRef != Stereo::implRef() && downRef != Stereo::implRef());
      }

      assert(upRef != Stereo::implRef() || downRef != Stereo::implRef());
      return true;
    }

  } // namespace impl

  template<typename EditableMoleculeType>
  bool Smiles::read(const std::string &smiles, EditableMoleculeType &mol, Stereochemistry &stereo)
  {
    // reset error
    m_error = Error();

    impl::SmileyCallback<EditableMoleculeType> callback(mol, stereo);
    Smiley::Parser<impl::SmileyCallback<EditableMoleculeType> > parser(callback);

    parser.disableExceptions(Smiley::InvalidChiralValence);

    try {
      parser.parse(smiles);
    } catch (Smiley::Exception &e) {
      std::ostringstream errorStream;
      if (e.type() == Smiley::Exception::SyntaxError)
        errorStream << "Syntax";
      else
        errorStream << "Semantics";
      errorStream << "Error: " << e.what() << "." << std::endl;
      errorStream << smiles << std::endl;
      for (std::size_t i = 0; i < e.pos(); ++i)
        errorStream << " ";
      for (std::size_t i = 0; i < e.length(); ++i)
        errorStream << "^";
      errorStream << std::endl;

      // Throw new exception with added details ...
      m_error = Error(errorStream.str());
      return false;
    }

    // add hydrogens
    for (auto &atom : get_atoms(mol)) {
      if (get_hydrogens(mol, atom) != 99)
        continue;
      set_hydrogens(mol, atom, 0);
      if (!Element::addHydrogens(get_element(mol, atom)))
        continue;

      int explicitH = 0;
      for (auto &nbr : get_nbrs(mol, atom))
        if (get_element(mol, nbr) == 1)
          ++explicitH;

      int valence = get_valence(mol, atom);
      assert(get_charge(mol, atom) == 0);
      //int expValence = Element::valence(get_element(mol, atom), get_charge(mol, atom), valence);
      int expValence = impl::smiles_valence(get_element(mol, atom), valence);

      //std::cout << "val: " << valence << "  deg: " << get_degree(mol, atom) << " explH: " << explicitH << " H: " << get_hydrogens(mol, atom) << std::endl;
      if (expValence > valence - explicitH)
        set_hydrogens(mol, atom, expValence - valence);
    }

    // add double bond stereochemistry
    //
    // ref = [ 0 1 2 3 ]
    //
    // 0        3
    //  \      /
    //   C == C
    //  /      \
    // 1        2
    for (auto &bond : get_bonds(mol)) {
      if (get_order(mol, bond) != 2)
        continue;

      // atom1 is the atom on the left side of the double bond
      auto atom1 = get_source(mol, bond);
      auto atom2 = get_target(mol, bond);
      if (get_index(mol, atom2) < get_index(mol, atom1))
        std::swap(atom1, atom2);

      // find the up/down bonds (i.e. all bonds not equal to the double bond)
      std::vector<Index> atom1bonds;
      for (auto &b : get_bonds(mol, atom1))
        if (b != bond)
          atom1bonds.push_back(get_index(mol, b));
      if (atom1bonds.size() == 0 || atom1bonds.size() > 2)
        continue;

      std::vector<Index> atom2bonds;
      for (auto &b : get_bonds(mol, atom2))
        if (b != bond)
          atom2bonds.push_back(get_index(mol, b));
      if (atom2bonds.size() == 0 || atom2bonds.size() > 2)
        continue;

      Stereo::Ref refs[4] = {Stereo::implRef(), Stereo::implRef(), Stereo::implRef(), Stereo::implRef()};

      if (!processUpDownBonds(mol, atom1, callback.upBonds, callback.downBonds, get_bond(mol, atom1bonds[0]),
            (atom1bonds.size() == 2) ? get_bond(mol, atom1bonds[1]) : molecule_traits<EditableMoleculeType>::null_bond(),
            refs[0], refs[1]))
        continue;

      if (!processUpDownBonds(mol, atom2, callback.upBonds, callback.downBonds, get_bond(mol, atom2bonds[0]),
            (atom2bonds.size() == 2) ? get_bond(mol, atom2bonds[1]) : molecule_traits<EditableMoleculeType>::null_bond(),
            refs[3], refs[2]))
        continue;

      stereo.add(StereoStorage(Stereo::CisTrans, get_index(mol, bond), refs, refs + 4));
    }

    return true;
  }

  template<typename MoleculeType>
  std::string Smiles::write(const MoleculeType &mol, const Stereochemistry &stereo,
      const std::map<Index, int> &atomClasses, int flags)
  {
    // determine ring closures and final output order
    impl::WriteSmilesPreprocessor<MoleculeType> preprocessor;
    depth_first_search(mol, preprocessor);

    impl::WriteSmilesVisitor<MoleculeType> smilesWriter(stereo, preprocessor.ringNumbers,
        preprocessor.outputIndices, atomClasses, flags);
    depth_first_search(mol, smilesWriter);

    return smilesWriter.smiles.str();
  }

  template<typename MoleculeType>
  std::string Smiles::write(const MoleculeType &mol, const Stereochemistry &stereo,
      const std::vector<Index> &order, const std::map<Index, int> &atomClasses,
      int flags)
  {
    // determine ring closures and final output order
    impl::WriteSmilesPreprocessor<MoleculeType> preprocessor;
    ordered_depth_first_search(mol, order, preprocessor);

    impl::WriteSmilesVisitor<MoleculeType> smilesWriter(stereo, preprocessor.ringNumbers,
        preprocessor.outputIndices, atomClasses, flags);
    ordered_depth_first_search(mol, order, smilesWriter);

    return smilesWriter.smiles.str();
  }

}

#endif
