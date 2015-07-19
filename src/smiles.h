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
       * @param smiles The SMILES string.
       * @param mol The molecule.
       */
      template<typename EditableMoleculeType>
      bool read(const std::string &smiles, EditableMoleculeType &mol,
          Stereochemistry &stereo);

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
       * @param flags The SMILES features to write (see Flags).
       *
       * @return The SMILES string.
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, int flags = All);

      /**
       * @brief Write a SMILES string for the molecule.
       *
       * @param mol The molecule.
       * @param atomClasses The atom classes.
       * @param flags The SMILES features to write (see Flags).
       *
       * @return The SMILES string.
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol,
          const std::map<Index, int> &atomClasses, int flags = All);

      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, const Stereochemistry &stereo,
          int flags = All);
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, const Stereochemistry &stereo,
          const std::map<Index, int> &atomClasses, int flags = All);

      /**
       * @brief Write a SMILES string for the molecule using a specified order.
       *
       * This version of write() can be used to write canonical SMILES if a
       * canonical atom order is used.
       *
       * @param mol The molecule.
       * @param order The atom order.
       * @param flags The SMILES features to write (see Flags).
       *
       * @return The SMILES string.
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, const std::vector<Index> &order,
          int flags = All);

      /**
       * @brief Write a SMILES string for the molecule using a specified order.
       *
       * This version of write() can be used to write canonical SMILES if a
       * canonical atom order is used.
       *
       * @param mol The molecule.
       * @param order The atom order.
       * @param atomClasses The atom classes.
       * @param flags The SMILES features to write (see Flags).
       *
       * @return The SMILES string.
       */
      template<typename MoleculeType>
      std::string write(const MoleculeType &mol, const std::vector<Index> &order,
          const std::map<Index, int> &atomClasses, int flags = All);

      /**
       * @brief Write a canonical SMILES string for the molecule.
       *
       * @param mol The molecule.
       * @param atomClasses The atom classes.
       * @param flags The SMILES features to write (see Flags).
       *
       * @return The SMILES string.
       */
      template<typename MoleculeType>
      std::string writeCanonical(const MoleculeType &mol, const std::map<Index, int> &atomClasses, int flags = All)
      {
        std::pair<std::vector<Index>, std::vector<unsigned long> > canon = canonicalize(mol,
            extended_connectivities(mol, DefaultAtomInvariant(DefaultAtomInvariant::Element)),
            DefaultAtomInvariant(DefaultAtomInvariant::All), DefaultBondInvariant(DefaultBondInvariant::All),
            connected_atom_components(mol), connected_bond_components(mol));
        return write(mol, canon.first, atomClasses, flags);
      }

      /**
       * @brief Write a canonical SMILES string for the molecule.
       *
       * @param mol The molecule.
       * @param flags The SMILES features to write (see Flags).
       *
       * @return The SMILES string.
       */
      template<typename MoleculeType>
      std::string writeCanonical(const MoleculeType &mol, int flags = All)
      {
        std::map<Index, int> atomClasses;
        return writeCanonical(mol, atomClasses, flags);
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
        typename molecule_traits<EditableMoleculeType>::bond_type bond = add_bond(mol, get_atom(mol, source), get_atom(mol, target));
        if (order == 5)
          set_aromatic(mol, bond, true);
        set_order(mol, bond, order);
      }

      void setChiral(int index, Smiley::Chirality chirality, const std::vector<int> &chiralNbrs)
      {
        Stereo::Ref refs[6];

        switch (chirality) {
          case Smiley::AntiClockwise:
          case Smiley::TH1:
            std::copy(chiralNbrs.begin(), chiralNbrs.end(), refs);
            for (std::size_t i = 0; i < 4; ++i)
              if (refs[i] == Smiley::implicitHydrogen())
                refs[i] = Stereo::implRef();
            stereo.add(StereoStorage::create(Stereo::Tetrahedral, index, refs, refs + 4));
            break;
          case Smiley::Clockwise:
          case Smiley::TH2:
            refs[0] = chiralNbrs[0];
            std::copy(chiralNbrs.rbegin(), chiralNbrs.rbegin() + 3, refs + 1);
            for (std::size_t i = 0; i < 4; ++i)
              if (refs[i] == Smiley::implicitHydrogen())
                refs[i] = Stereo::implRef();
            stereo.add(StereoStorage::create(Stereo::Tetrahedral, index, refs, refs + 4));
            break;
        }
      }

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
        if (stereo.isTetrahedral(get_index(mol, atom))) {
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

          if (refs.size() == 3) {
            if (prev != molecule_traits<MoleculeType>::null_atom())
              refs.insert(refs.begin() + 1, Stereo::implRef());
            else
              refs.insert(refs.begin(), Stereo::implRef());
          }

          assert(refs.size() == 4);

          //std::cout << stereo.tetrahedral(get_index(mol, atom)) << std::endl;
          //std::cout << "refs: " << refs << std::endl;

          if (tetrahedral_class(stereo.tetrahedral(get_index(mol, atom)), refs[0], refs[1], refs[2], refs[3]) == Stereo::TH1)
            smiles << "@";
          else
            smiles << "@@";
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
        if (stereo.isTetrahedral(get_index(mol, atom)))
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

      const Stereochemistry &stereo;
      const std::map<atom_type, std::vector<RingNumber> > &ringNumbers; // atom -> ring numbers
      std::map<int, std::pair<int, int> > orderedRingNumbers; // ring number -> (final ring number, bond order)
      std::map<int, int> closureBondOrders;
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
  std::string Smiles::write(const MoleculeType &mol, const std::map<Index, int> &atomClasses,
      int flags)
  {
    Stereochemistry stereo;
    return write(mol, stereo, atomClasses, flags);
  }

  template<typename MoleculeType>
  std::string Smiles::write(const MoleculeType &mol, const Stereochemistry &stereo, int flags)
  {
    std::map<Index, int> atomClasses;
    return write(mol, stereo, atomClasses, flags);
  }

  template<typename MoleculeType>
  std::string Smiles::write(const MoleculeType &mol, int flags)
  {
    Stereochemistry stereo;
    std::map<Index, int> atomClasses;
    return write(mol, stereo, atomClasses, flags);
  }

  template<typename MoleculeType>
  std::string Smiles::write(const MoleculeType &mol, const std::vector<Index> &order,
      const std::map<Index, int> &atomClasses, int flags)
  {
    impl::WriteSmilesPreprocessor<MoleculeType> preprocessor;
    ordered_depth_first_search(mol, order, preprocessor);

    Stereochemistry stereo;
    impl::WriteSmilesVisitor<MoleculeType> smilesWriter(stereo, preprocessor.ringNumbers,
        preprocessor.outputIndices, atomClasses, flags);
    ordered_depth_first_search(mol, order, smilesWriter);

    return smilesWriter.smiles.str();
  }

  template<typename MoleculeType>
  std::string Smiles::write(const MoleculeType &mol, const std::vector<Index> &order, int flags)
  {
    std::map<Index, int> atomClasses;
    return write(mol, order, atomClasses, flags);
  }

}

#endif
