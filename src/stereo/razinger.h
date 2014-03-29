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
#ifndef HELIUM_RAZINGER_H
#define HELIUM_RAZINGER_H

#include <Helium/stereo/stereounit.h>
#include <Helium/molecule.h>
#include <Helium/tie.h>
#include <Helium/contract.h>

#include <iostream>

#define DEBUG_RAZINGER 1

namespace Helium {

  /**
   * @brief Struct representing a single stereogenic unit.
   */
  struct StereoUnit
  {
    /**
     * Default constructor creating an invalid StereoUnit. The unit can be
     * made valid by setting the type, id and para data members.
     */
    StereoUnit() : type(static_cast<Stereo::Type>(0)), index(Stereo::nullRef()), para(false)
    {
    }

    /**
     * Constructor specifying all data to create a valid StereoUnit.
     */
    StereoUnit(Stereo::Type _type, unsigned long _index, bool _para = false) :
        type(_type), index(_index), para(_para)
    {
    }

    Stereo::Type type; //!< the type for this stereogenic unit
    unsigned long index; //! the atom/bond (depends on type) index
    bool para; //! para- (=ressemble) or true-stereocenter
  };

  /**
   * @brief A single set of StereoUnit objects.
   *
   * This type can be used to represent all stereogenic units in a molecule and
   * is used as return type of find_stereogenic_units(). This set is also the input
   * for many functions requiring this information (e.g. StereoFrom2D, ...).
   */
  typedef std::vector<StereoUnit> StereoUnitSet;

  /**
   * @brief A set of sets of StereoUnit objects.
   *
   * This type is used for cases where there is some relationship between
   * individual StereoUnit objects.
   */
  typedef std::vector<StereoUnitSet> StereoUnitSetOfSets;



  namespace impl {

    /**
     * These are the different types of stereocenter classification used throughout
     * this file.
     */
    enum NeighborSymmetryClasses
    {
      // Tetrahedral
      T1234 = 1234, // 4 different symmetry classes
      T1123 = 1123, // 3 different symmetry classes, 1 class duplicated (2 times)
      T1122 = 1122, // 2 different symmetry classes, 1 class duplicated (3 times)
      T1112 = 1112, // 2 different symmetry classes, each class duplicated (2 times)
      T1111 = 1111, // 1 symmetry class, duplictaed 4 times
      // CisTrans
      C12 = 12, // 2 different symmetry classes
      C11 = 11 // the same symmetry class
    };

    /**
     * Perform a quick check for tetrahedral stereo centers. This function is
     * used by find_stereogenic_units to return quickly if there are no
     * potential stereogenic units.
     */
    template<typename MoleculeType>
    bool may_have_tetrahedral_stereo(const MoleculeType &mol)
    {
      FOREACH_ATOM_T (atom, mol, MoleculeType)
        if (get_heavy_degree(mol, *atom) >= 3)
          return true;
      return false;
    }

    template<typename MoleculeType>
    bool may_have_cistrans_stereo(const MoleculeType &mol)
    {
      FOREACH_BOND_T (bond, mol, MoleculeType)
        if (get_order(mol, *bond) == 2)
          return true;
      return false;
    }

    /**
     * Check if the specified atom is a potential stereogenic atom.
     *
     * Criteria:
     * - not connected to more than 4 atoms
     * - at least 3 "heavy" neighbors
     *
     * Nitrogen (neutral) is treated as a special case since the barrier of inversion is
     * low in many cases making the atom non-stereogenic. Only bridge-head
     * nitrogen atoms (i.e. nitrogen has 3 neighbors in rings) will be
     * considered stereogenic.
     */
    template<typename MoleculeType, typename AtomType>
    bool is_potential_tetrahedral(const MoleculeType &mol, const RingSet<MoleculeType> &rings, AtomType atom)
    {
      // consider only potential steroecenters
      /* FIXME
      if ((atom->GetHyb() != 3 && !(atom->GetHyb() == 5 && atom->IsPhosphorus()))
          || atom->GetImplicitValence() > 4 || atom->GetHvyValence() < 3 || atom->GetHvyValence() > 4)
      */
      if (get_valence(mol, atom) > 4 || get_heavy_degree(mol, atom) < 3 || get_heavy_degree(mol, atom) > 4)
        return false;
      // skip non-chiral N
      if (is_nitrogen(mol, atom) && get_charge(mol, atom) == 0 && rings.numRingBonds(atom) < 3)
        return false;

      if (is_carbon(mol, atom)) {
        if (get_charge(mol, atom))
          return false;
        if (atom_has_nbr(mol, atom, atom_and_predicates(mol, element_eq_predicate(mol, 26), degree_gt_predicate(mol, 7))))
          return false;
      }

      return true;
    }

    template<typename MoleculeType>
    struct StereoRing;

    template<typename MoleculeType>
    struct ParaAtom
    {
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

      ParaAtom(Index idx) : inIdx(idx)
      {
      }

      atom_type center(const MoleculeType &mol) const
      {
        return get_atom(mol, inIdx);
      }

      bool isInRing(const StereoRing<MoleculeType> &ring) const
      {
        for (std::size_t i = 0; i < ring.paraAtoms.size(); ++i)
          if (ring.paraAtoms[i].inIdx == inIdx)
            return true;
        return false;
      }

      union {
        Index index, inIdx, outIdx;
      };
      std::vector<atom_type> insideNbrs, outsideNbrs;
    };

    template<typename MoleculeType>
    struct ParaBond
    {
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      ParaBond(Index _bondIndex, Index _inIdx, Index _outIdx)
        : index(_bondIndex), inIdx(_inIdx), outIdx(_outIdx)
      {
      }

      bond_type center(const MoleculeType &mol) const
      {
        return get_bond(mol, index);
      }

      bool isInRing(const StereoRing<MoleculeType> &ring) const
      {
        for (std::size_t i = 0; i < ring.paraBonds.size(); ++i)
          if (ring.paraBonds[i].inIdx == inIdx)
            return true;
        return false;
      }

      Index index;
      Index inIdx, outIdx;
      std::vector<atom_type> insideNbrs, outsideNbrs;
    };

    template<typename MoleculeType>
    struct StereoRing
    {
      StereoRing() : trueCount(0)
      {
      }

      std::vector<ParaAtom<MoleculeType> > paraAtoms;
      std::vector<ParaBond<MoleculeType> > paraBonds;
      unsigned int trueCount;
    };

    /**
     * Merge the rings in a molecule and return the result as std::vector<bool> objects.
     * Rings are merged if they share at least one atom (e.g. bridged, spiro,
     * adjacent, ...).
     */
    template<typename MoleculeType>
    std::vector<std::vector<bool> > merge_rings(const MoleculeType &mol, const RingSet<MoleculeType> &rings, const std::vector<unsigned long> &ec)
    {
      std::vector<std::vector<bool> > result;

      for (std::size_t i = 0; i < rings.size(); ++i) {
        // check if ring shares atom with previously found ring
        bool found = false;
        for (std::size_t j = 0; j < result.size(); ++j) {
          std::vector<Index> shared;
          // foreach ring atom
          for (std::size_t k = 0; k < rings.ring(i).size(); ++k) {
            // check if the ring atom is in the current result bitvec
            if (result[j][get_index(mol, rings.ring(i).atom(k))]) {
              shared.push_back(get_index(mol, rings.ring(i).atom(k)));
            }
          }

          if (shared.size() > 1) {
            found = true;
          } else if (shared.size() == 1) {
            int classification = classify_tetrahedral_nbr_symmetry_classes(mol, get_atom(mol, shared[0]), ec);
            if (classification == T1122 || classification == T1111)
              found = true;
          }

          if (found) {
            // add bits for the atoms in the ring
            for (std::size_t l = 0; l < rings.ring(i).size(); ++l)
              result[j][get_index(mol, rings.ring(i).atom(l))] = true;
            break;
          }
        }

        // add the ring as a new bitvec if it shares no atom with a previous ring
        if (!found) {
          std::vector<bool> r;
          for (std::size_t l = 0; l < rings.ring(i).size(); ++l)
            r[get_index(mol, rings.ring(i).atom(l))] = true;
          result.push_back(r);
        }
      }

      return result;
    }

    /**
     * Classify the tetrahedral atom using the NeighborSymmetryClasses types.
     */
    template<typename MoleculeType, typename AtomType>
    int classify_tetrahedral_nbr_symmetry_classes(const MoleculeType &mol, AtomType atom, const std::vector<unsigned long> &ec)
    {
      std::vector<unsigned long> nbrClasses, nbrClassesCopy, uniqueClasses;
      FOREACH_NBR_T (nbr, atom, mol, MoleculeType)
        nbrClasses.push_back(ec.at(get_index(mol, *nbr)));
      // add an implicit ref if there are only 3 explicit
      if (nbrClasses.size() == 3)
        nbrClasses.push_back(Stereo::implRef());

      // use some STL to work out the number of unique classes
      nbrClassesCopy = nbrClasses; // keep copy for count below
      std::sort(nbrClasses.begin(), nbrClasses.end());
      std::vector<unsigned long>::iterator endLoc = std::unique(nbrClasses.begin(), nbrClasses.end());
      std::copy(nbrClasses.begin(), endLoc, std::back_inserter(uniqueClasses));

      switch (uniqueClasses.size()) {
        case 4:
          return T1234; // e.g. 1 2 3 4
        case 3:
          return T1123; // e.g. 1 1 2 3
        case 2:
          // differentiate between T1122 and T1112
          if (std::count(nbrClassesCopy.begin(), nbrClassesCopy.end(), uniqueClasses.at(0)) == 2)
            return T1122; // e.g. 1 1 2 2
          else
            return T1112; // e.g. 1 1 1 2
        case 1:
        default:
          return T1111; // e.g. 1 1 1 1
      }
    }

    /**
     * Classify the cis/trans bond using the NeighborSymmetryClasses types.
     */
    template<typename MoleculeType, typename BondType, typename AtomType>
    int classify_cis_trans_nbr_symmetry_classes(const MoleculeType &mol, const std::vector<unsigned long> &ec, BondType bond, AtomType atom)
    {
      std::vector<unsigned long> nbrClasses, uniqueClasses;
      FOREACH_NBR_T (nbr, atom, mol, MoleculeType) {
        if (get_index(mol, *nbr) != get_index(mol, get_other(mol, bond, atom)))
          nbrClasses.push_back(ec.at(get_index(mol, *nbr)));
      }

      if (nbrClasses.size() == 1)
        nbrClasses.push_back(Stereo::implRef());

      if (nbrClasses.at(0) == nbrClasses.at(1))
        return C11; // e.g. 1 1
      else
        return C12; // e.g. 1 2
    }

    /**
     * Helper function for FindStereogenicUnits using automorphisms.
     *
     * Find the duplicated symmetry class for neighbors of atom. This method only works if there is
     * only one duplicated symmetry class (i.e. T1123, T1112, T1111).
     */
    template<typename MoleculeType, typename AtomType>
    unsigned long find_duplicated_symmetry_class(const MoleculeType &mol, AtomType atom, const std::vector<unsigned long> &ec)
    {
      // find the duplicated symmetry class
      unsigned long duplicatedSymClass = -1; // OBGraphSym::NoSymmetryClass; // FIXME
      std::vector<unsigned long> nbrSymClasses;
      FOREACH_NBR_T (nbr, atom, mol, MoleculeType) {
        nbrSymClasses.push_back(ec.at(get_index(mol, *nbr)));
      }

      for (std::size_t i = 0; i < nbrSymClasses.size(); ++i) {
        if (std::count(nbrSymClasses.begin(), nbrSymClasses.end(), nbrSymClasses.at(i)) >= 2) {
          duplicatedSymClass = nbrSymClasses.at(i);
          break;
        }
      }

      return duplicatedSymClass;
    }

    /**
     * Helper function for FindStereogenicUnits using automorphisms.
     *
     * Find the duplicated symmetry classes for neighbors of atom. This method only works for the
     * T1122 case.
     */
    template<typename MoleculeType, typename AtomType>
    void find_duplicated_symmetry_classes(const MoleculeType &mol, AtomType atom, const std::vector<unsigned long> &ec,
        unsigned long &duplicated1, unsigned long &duplicated2)
    {
      std::vector<unsigned long> nbrSymClasses;
      FOREACH_NBR_T (nbr, atom, mol, MoleculeType)
        nbrSymClasses.push_back(ec.at(get_index(mol, *nbr)));
      std::sort(nbrSymClasses.begin(), nbrSymClasses.end());
      duplicated1 = nbrSymClasses[0];
      duplicated2 = nbrSymClasses[2];
    }

    /**
     * Find an atom with the specified symmetry class. The first found atom is returned
     * or 0 when there is no such atom. This function is intended to be used in cases
     * where any atom with the specified symmetry class can be used. For example, when
     * checking a fragments for stereocenters, the result will be the same for any atom
     * with a specified (duplicated) symmetry class.
     */
    template<typename MoleculeType, typename AtomType>
    AtomType find_atom_with_symmetry_class(const MoleculeType &mol, AtomType atom, unsigned long symClass, const std::vector<unsigned long> &ec)
    {
      AtomType ligandAtom = molecule_traits<MoleculeType>::null_atom();
      FOREACH_NBR_T (nbr, atom, mol, MoleculeType)
        if (ec.at(get_index(mol, *nbr)) == symClass)
          ligandAtom = *nbr;
      return ligandAtom;
    }

    /**
     * Helper function for getFragment below.
     */
    template<typename MoleculeType, typename AtomType>
    void add_nbrs(const MoleculeType &mol, std::vector<bool> &fragment, AtomType atom, AtomType skip)
    {
      FOREACH_NBR_T (nbr, atom, mol, MoleculeType) {
        // don't pass through skip
        if (get_index(mol, *nbr) == get_index(mol, skip))
          continue;
        // skip visited atoms
        if (fragment[get_index(mol, *nbr)])
          continue;
        // add the neighbor atom to the fragment
        fragment[get_index(mol, *nbr)] = true;
        // recurse...
        add_nbrs(mol, fragment, *nbr, skip);
      }
    }

    /**
     * Create a bit vector object with bets set for the fragment consisting of all
     * atoms for which there is a path to atom without going through skip.
     */
    template<typename MoleculeType, typename AtomType>
    std::vector<bool> get_fragment(const MoleculeType &mol, AtomType atom, AtomType skip)
    {
      std::vector<bool> fragment(num_atoms(mol));
      fragment[get_index(mol, atom)] = true;
      // start the recursion
      add_nbrs(mol, fragment, atom, skip);
      return fragment;
    }

    /**
     * Check if the specified stereogenic unit is in a fragment.
     */
    template<typename MoleculeType>
    bool is_unit_in_fragment(const MoleculeType &mol, const StereoUnit &unit, const std::vector<bool> &fragment)
    {
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      if (unit.type == Stereo::Tetrahedral) {
        if (fragment[unit.index])
          return true;
      } else if(unit.type == Stereo::CisTrans) {
        bond_type bond = get_bond(mol, unit.index);
        atom_type source = get_source(mol, bond);
        atom_type target = get_target(mol, bond);

        if (fragment[get_index(mol, source)] || fragment[get_index(mol, target)])
          return true;
      }

      return false;
    }

    /**
     * Helper function to determine if a stereogenic center with duplicated symmetry classes
     * really is a stereogenic center.
     *
     * Check if the ligandAtom's fragment (see getFragment()) contains at least one
     * true- or 2 para-stereocenter. This is rule 2a and rule 3 in the Razinger
     * paper on stereoisomer generation.
     */
    template<typename MoleculeType, typename AtomType>
    bool contains_at_least_1true_2para(const MoleculeType &mol, AtomType ligandAtom, AtomType atom,
        const StereoUnitSet &units, const RingSet<MoleculeType> &rings)
    {
      // check if ligand contains at least:
      // - 1 true-stereocenter
      // - 2 para-stereocenters
      std::vector<bool> ligand = get_fragment(mol, ligandAtom, atom);
      bool foundTrueStereoCenter = false;
      int paraStereoCenterCount = 0;
      for (StereoUnitSet::const_iterator unit = units.begin(); unit != units.end(); ++unit) {
        if (is_unit_in_fragment(mol, *unit, ligand)) {
          if (unit->para) {
            paraStereoCenterCount++;
          } else {
            foundTrueStereoCenter = true;
          }
        }
      }

      if (foundTrueStereoCenter || paraStereoCenterCount >= 2)
        return true;
      if (rings.isAtomInRing(ligandAtom) && rings.isAtomInRing(atom) && paraStereoCenterCount)
        return true;
      return false;
    }

    /**
     * Helper function to determine if a stereogenic center with duplicated symmetry classes
     * really is a stereogenic center.
     *
     * Check if the ligandAtom's fragment (see getFragment()) contains at least one
     * true- or 2 separate assemblies of at least 2 para-stereocenter. This is rule
     * 2b in the Razinger paper on stereoisomer generation.
     */
    template<typename MoleculeType, typename AtomType>
    bool contains_at_least_2true_2paraAssemblies(const MoleculeType &mol, AtomType ligandAtom, AtomType atom,
        const StereoUnitSet &units, const std::vector<std::vector<bool> > &mergedRings)
    {
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

      // check if ligand contains at least:
      // - 2 true-stereocenter
      // - 2 separate para-stereocenters assemblies
      std::vector<bool> ligand = get_fragment(mol, ligandAtom, atom);
      int trueStereoCenterCount = 0;
      std::vector<std::size_t> ringIndices;
      for (StereoUnitSet::const_iterator unit = units.begin(); unit != units.end(); ++unit) {
        if (unit->type == Stereo::Tetrahedral) {
          if (ligand[unit->index]) {
            if (unit->para) {
              AtomType paraAtom = get_atom(mol, unit->index);
              for (std::size_t ringIdx = 0; ringIdx < mergedRings.size(); ++ringIdx) {
                if (mergedRings.at(ringIdx)[get_index(mol, paraAtom)])
                  if (std::find(ringIndices.begin(), ringIndices.end(), ringIdx) == ringIndices.end())
                    ringIndices.push_back(ringIdx);
              }
            } else
              trueStereoCenterCount++;
          }
        } else if(unit->type == Stereo::CisTrans) {
          bond_type bond = get_bond(mol, unit->index);
          AtomType begin = get_source(mol, bond);
          AtomType end = get_target(mol, bond);
          if (ligand[get_index(mol, begin)] || ligand[get_index(mol, end)]) {
            if (unit->para) {
              for (std::size_t ringIdx = 0; ringIdx < mergedRings.size(); ++ringIdx) {
                if (mergedRings.at(ringIdx)[get_index(mol, begin)] || mergedRings.at(ringIdx)[get_index(mol, end)]) {
                  if (std::find(ringIndices.begin(), ringIndices.end(), ringIdx) == ringIndices.end()) {
                    ringIndices.push_back(ringIdx);
                  }
                }
              }
            } else {
              trueStereoCenterCount++;
            }
          }
        }
      }

      if (trueStereoCenterCount >= 2 || ringIndices.size() >= 2)
        return true;
      return false;
    }

    template<typename MoleculeType, typename Type>
    bool check_ligands(const MoleculeType &mol, const Type &currentPara, const StereoUnitSet &units)
    {
      if (currentPara.outsideNbrs.size() == 1) {
        //cout << "OK: " << __LINE__ << endl;
        return true;
      }
      //assert(mol->GetAtom(currentPara.outIdx));
      std::vector<bool> ligand = get_fragment(mol, currentPara.outsideNbrs[0], get_atom(mol, currentPara.outIdx));
      for (StereoUnitSet::const_iterator unit = units.begin(); unit != units.end(); ++unit) {
        if (is_unit_in_fragment(mol, *unit, ligand)) {
          //cout << "OK: " << __LINE__ << endl;
          return true;
        }
      }
      //cout << "NOT OK: " << __LINE__ << endl;
      return false;
    }

    template<typename MoleculeType, typename Type>
    bool apply_rule1(const MoleculeType &mol, const Type &currentPara, const std::vector<unsigned long> &symmetry_classes,
        const std::vector<StereoRing<MoleculeType> > &rings, std::vector<bool> &visitedRings, const StereoUnitSet &units,
        std::vector<Index> stereoAtoms)
    {
      bool foundRing = false;
      unsigned int idx = currentPara.inIdx;

      /*
         for (std::size_t i = 0; i < visitedRings.size(); ++i)
         if (visitedRings[i])
         cout << "  ";
         cout << "ApplyRule1(" << currentPara.inIdx << ", " << currentPara.outIdx << ", outside = " << currentPara.outsideNbrs.size() << ")" << endl;
         */

      for (std::size_t i = 0; i < rings.size(); ++i) {
        // skip visited rings
        if (visitedRings[i])
          continue;

        // Check if currentPara is in this ring
        if (!currentPara.isInRing(rings[i]))
          continue;

        //
        // A new ring containing currentPara is found
        //
        foundRing = true;

        // if there are one or more true stereo centers, currentPara is a stereo center
        if (rings[i].trueCount) {
          //cout << "OK: " << __LINE__ << endl;
          return true;
        }

        // check if there is at least one other potential atom
        for (std::size_t j = 0; j < rings[i].paraAtoms.size(); ++j) {
          const ParaAtom<MoleculeType> &paraAtom = rings[i].paraAtoms[j];
          // skip idx
          if (paraAtom.inIdx == idx)
            continue;
          // there is another atom already identified as stereo atom
          if (std::find(stereoAtoms.begin(), stereoAtoms.end(), paraAtom.inIdx) != stereoAtoms.end()) {
            //cout << "OK: " << __LINE__ << endl;
            return true;
          }

          if (paraAtom.outsideNbrs.size() == 1) {
            // only 1 ring substituent, the other is implicit H -> topologically different
            //cout << "OK: " << __LINE__ << endl;
            return true;
          } else {
            if (paraAtom.outsideNbrs.size() != 2)
              return false;
            // two ring substituents, need to check for topological difference
            if (symmetry_classes[get_index(mol, paraAtom.outsideNbrs[0])] != symmetry_classes[get_index(mol, paraAtom.outsideNbrs[1])]) {
              // they are different
              //cout << "OK: " << __LINE__ << endl;
              return true;
            } else {
              // they are the same and they might also be in a ring -> apply rule 1 recursive
              visitedRings[i] = true;
              if (apply_rule1(mol, paraAtom, symmetry_classes, rings, visitedRings, units, stereoAtoms)) {
                //cout << "OK: " << __LINE__ << endl;
                return true;
              }
            }
          }
        }
        // check if there is at least one other potential bond
        for (std::size_t j = 0; j < rings[i].paraBonds.size(); ++j) {
          const ParaBond<MoleculeType> &paraBond = rings[i].paraBonds[j];
          // skip idx
          if (paraBond.inIdx == idx)
            continue;
          // there is another atom already identified as stereo atom
          if (std::find(stereoAtoms.begin(), stereoAtoms.end(), paraBond.inIdx) != stereoAtoms.end()) {
            //cout << "OK: " << __LINE__ << endl;
            return true;
          }

          if (paraBond.outsideNbrs.size() == 1) {
            // only 1 ring substituent, the other is implicit H -> topologically different
            //cout << "OK: " << __LINE__ << endl;
            return true;
          } else {
            if (paraBond.outsideNbrs.size() != 2)
              continue;
            // two ring substituents, need to check for topological difference
            if (symmetry_classes[get_index(mol, paraBond.outsideNbrs[0])] != symmetry_classes[get_index(mol, paraBond.outsideNbrs[1])]) {
              // they are different
              //cout << "OK: " << __LINE__ << endl;
              return true;
            } else {
              // they are the same and they might also be in a ring -> apply rule 1 recursive
              visitedRings[i] = true;
              if (apply_rule1(mol, paraBond, symmetry_classes, rings, visitedRings, units, stereoAtoms)) {
                //cout << "OK: " << __LINE__ << endl;
                return true;
              }
            }
          }
        }

      }

      // if a non-visited ring was found and true was not returned -> it does not
      // contain any stereocenters other than idx
      if (foundRing) {
        //cout << "NOT OK: " << __LINE__ << endl;
        return false;
      }

      //cout << "NOT OK: " << __LINE__ << endl;
      return false;
    }

    template<typename MoleculeType>
    void start_rule1(const MoleculeType &mol, const std::vector<unsigned long> &symmetry_classes,
        const std::vector<StereoRing<MoleculeType> > &rings, StereoUnitSet &units,
        std::vector<Index> &stereoAtoms)
    {
      for (std::size_t i = 0; i < rings.size(); ++i) {
        //cout << "Checking ring: " << i << endl;

        // tetrahedral atoms
        for (std::size_t j = 0; j < rings[i].paraAtoms.size(); ++j) {
          const ParaAtom<MoleculeType> &paraAtom = rings[i].paraAtoms[j];
          // skip the atom if it is already in stereoAtoms
          if (std::find(stereoAtoms.begin(), stereoAtoms.end(), paraAtom.inIdx) != stereoAtoms.end())
            continue;

          std::vector<bool> visitedRings(rings.size(), false);
          //visitedRings[i] = true;
          if (apply_rule1(mol, paraAtom, symmetry_classes, rings, visitedRings, units, stereoAtoms)) {
            bool isStereoUnit = false;
            if (paraAtom.outsideNbrs.size() == 1)
              isStereoUnit = true;
            if (paraAtom.outsideNbrs.size() == 2) {
              if (symmetry_classes[get_index(mol, paraAtom.outsideNbrs[0])] == symmetry_classes[get_index(mol, paraAtom.outsideNbrs[1])]) {
                // check for spiro atom
                bool isSpiro = false;
                for (std::size_t k = 0; k < rings[i].paraAtoms.size(); ++k) {
                  const ParaAtom<MoleculeType> &paraAtom2 = rings[i].paraAtoms[k];
                  if (paraAtom.inIdx == paraAtom2.outIdx && paraAtom.insideNbrs == paraAtom2.outsideNbrs) {
                    isSpiro = true;
                    if (apply_rule1(mol, paraAtom2, symmetry_classes, rings, visitedRings, units, stereoAtoms))
                      isStereoUnit = true;
                  }
                }
                if (!isSpiro)
                  isStereoUnit = check_ligands(mol, paraAtom, units);
                //cout << "isStereoUnit = " << isStereoUnit << endl;
              } else {
                isStereoUnit = true;
              }
            }

            if (isStereoUnit) {
              stereoAtoms.push_back(paraAtom.inIdx);
              units.push_back(StereoUnit(Stereo::Tetrahedral, paraAtom.index, true));
            }
          }

        }

        // cistrans bonds
        for (std::size_t j = 0; j < rings[i].paraBonds.size(); ++j) {
          const ParaBond<MoleculeType> &paraBond = rings[i].paraBonds[j];
          // skip the atom if it is already in stereoAtoms
          if (std::find(stereoAtoms.begin(), stereoAtoms.end(), paraBond.inIdx) != stereoAtoms.end())
            continue;

          std::vector<bool> visitedRings(rings.size(), false);
          //visitedRings[i] = true;
          if (apply_rule1(mol, paraBond, symmetry_classes, rings, visitedRings, units, stereoAtoms)) {
            bool isStereoUnit = false;
            if (paraBond.outsideNbrs.size() == 1)
              isStereoUnit = true;
            if (paraBond.outsideNbrs.size() == 2) {
              if (symmetry_classes[get_index(mol, paraBond.outsideNbrs[0])] == symmetry_classes[get_index(mol, paraBond.outsideNbrs[1])]) {
                // check for spiro bond
                bool isSpiro = false;
                for (std::size_t k = 0; k < rings[i].paraBonds.size(); ++k) {
                  const ParaBond<MoleculeType> &paraBond2 = rings[i].paraBonds[k];
                  if (paraBond.inIdx == paraBond2.outIdx && paraBond.insideNbrs == paraBond2.outsideNbrs) {
                    isSpiro = true;
                    if (apply_rule1(mol, paraBond2, symmetry_classes, rings, visitedRings, units, stereoAtoms))
                      isStereoUnit = true;
                  }
                }
                if (!isSpiro)
                  isStereoUnit = check_ligands(mol, paraBond, units);
                //cout << "isStereoUnit = " << isStereoUnit << endl;
              } else {
                isStereoUnit = true;
              }
            }

            if (isStereoUnit) {
              stereoAtoms.push_back(paraBond.inIdx);
              stereoAtoms.push_back(paraBond.outIdx);
              units.push_back(StereoUnit(Stereo::CisTrans, paraBond.index, true));
            }
          }

        }


      }

    }

  }

  template<typename MoleculeType>
  std::vector<StereoUnit> find_stereogenic_units(const MoleculeType &mol,
      const RingSet<MoleculeType> &rings, const std::vector<bool> &aromaticBonds,
      const std::vector<unsigned long> &ec)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;

    PRE(ec.size() == num_atoms(mol));

    // quick check to see if there are any potential stereogenic units
    if (!impl::may_have_tetrahedral_stereo(mol) && !impl::may_have_cistrans_stereo(mol))
      return std::vector<StereoUnit>();

    std::vector<StereoUnit> units;

    // para-stereocenters candidates
    std::vector<Index> stereoAtoms; // Tetrahedral = idx, CisTrans = begin & end idx
    std::vector<Index> paraAtoms;
    std::vector<Index> paraBonds;

    /**
     * true Tetrahedral stereocenters:
     * - have four different symmetry classes for the ligands to the central atom
     */
    bool ischiral;
    FOREACH_ATOM_T (atom, mol, MoleculeType) {
      if (!impl::is_potential_tetrahedral(mol, rings, *atom))
        continue;

      // list containing neighbor symmetry classes
      std::vector<unsigned long> tlist;
      ischiral = true;

      // check neighbors to see if this atom is stereogenic
      FOREACH_NBR_T (nbr, *atom, mol, MoleculeType) {
        // check if we already have a neighbor with this symmetry class
        std::vector<unsigned long>::iterator k;
        for (k = tlist.begin(); k != tlist.end(); ++k)
          if (ec[get_index(mol, *nbr)] == *k) {
            ischiral = false;
            // if so, might still be a para-stereocenter
            paraAtoms.push_back(get_index(mol, *atom));
          }

        if (ischiral)
          // keep track of all neighbors, so we can detect duplicates
          tlist.push_back(ec[get_index(mol, *nbr)]);
        else
          break;
      }

      if (ischiral) {
        // true-stereocenter found
        stereoAtoms.push_back(get_index(mol, *atom));
        units.push_back(StereoUnit(Stereo::Tetrahedral, get_index(mol, *atom)));
      }
    }

    /**
     * true CisTrans stereocenters:
     * - each terminal has two different symmetry classes for it's ligands
     */
    bool isCisTransBond;
    FOREACH_BOND_T (bond, mol, MoleculeType) {
      if (rings.isBondInRing(*bond) && aromaticBonds[get_index(mol, *bond)])
        continue; // Exclude C=C in phenyl rings for example

      if (get_order(mol, *bond) == 2) {
        atom_type begin = get_source(mol, *bond);
        atom_type end = get_target(mol, *bond);
        if (begin == molecule_traits<MoleculeType>::null_atom() || end == molecule_traits<MoleculeType>::null_atom())
          continue;

        if (get_valence(mol, begin) > 3 || get_valence(mol, end) > 3)
          continue; // e.g. C=Ru where the Ru has four substituents

        // Needs to have at least one explicit single bond at either end
        // FIXME: timvdm: what about C=C=C=C
        if (!atom_has_bond(mol, begin, order_eq_predicate(mol, 1)) || !atom_has_bond(mol, end, order_eq_predicate(mol, 1)))
          continue;

        isCisTransBond = true;

        if (get_degree(mol, begin) == 2) {
          // Begin atom has two explicit neighbors. One is the end atom. The other should
          // be a heavy atom - this is what we test here.
          // (There is a third, implicit, neighbor which is either a hydrogen
          // or a lone pair.)
          if (atom_has_nbr(mol, begin, element_eq_predicate(mol, 1)))
            isCisTransBond = false;
        } else if (get_degree(mol, begin) == 3) {
          std::vector<Index> tlist;

          FOREACH_NBR_T (nbr, begin, mol, MoleculeType) {
            // skip end atom
            if (get_index(mol, *nbr) == get_index(mol, end))
              continue;
            // do we already have an atom with this symmetry class?
            if (tlist.size()) {
              // compare second with first
              if (ec[get_index(mol, *nbr)] == tlist.at(0)) {
                isCisTransBond = false;
                // if same, might still be a para-stereocenter
                paraBonds.push_back(get_index(mol, *bond));
              }
              break;
            }

            // save first symmetry class
            tlist.push_back(ec[get_index(mol, *nbr)]);
          }
        } else {
          // Valence is not 2 or 3, for example SR3=NR
          isCisTransBond = false;
        }

        if (!isCisTransBond)
          continue;

        if (get_degree(mol, end) == 2) {
          // see comment above for begin atom
          if (atom_has_nbr(mol, end, element_eq_predicate(mol, 1)))
            isCisTransBond = false;
        } else if (get_degree(mol, end) == 3) {
          std::vector<Index> tlist;

          FOREACH_NBR_T (nbr, end, mol, MoleculeType) {
            // skip end atom
            if (get_index(mol, *nbr) == get_index(mol, begin))
              continue;
            // do we already have an atom with this symmetry class?
            if (tlist.size()) {
              // compare second with first
              if (ec[get_index(mol, *nbr)] == tlist.at(0)) {
                // if same, might still be a para-stereocenter
                paraBonds.push_back(get_index(mol, *bond));
                isCisTransBond = false;
              }
              break;
            }

            // save first symmetry class
            tlist.push_back(ec[get_index(mol, *nbr)]);
          }
        } else {
          // Valence is not 2 or 3, for example SR3=NR
          isCisTransBond = false;
        }

        if (isCisTransBond)
          // true-stereocenter found
          units.push_back(StereoUnit(Stereo::CisTrans, get_index(mol, *bond)));
      }
    }

    /**
     * Apply rule 1 from the Razinger paper recusively:
     *
     * All rings are merged "mergedRings". A merged ring is simply a fragment consisting
     * of all atoms of a ring system (bridged, spiro, adjacent, ...). If two rings in the
     * SSSR set share an atom, they are merged.
     *
     * Each merged must at least have two para-stereocenters (or 1 true + 1 para) in order
     * for the para-stereocenter to be valid. This is repeated until no new stereocenters
     * are identified.
     *
     * rule 1a for double bonds:
     * - bond atom in ring has two identical symmetry classes for it's neighbor atoms (-> para)
     * - other bond atom:
     *   - has two different symmetry classes for it's neighbours -> new stereocenter
     *   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
     *
     * rule 1b for tetracoord atoms:
     * - at least two neighbour symmetry classes are the same (-> para)
     * - other pair:
     *   - has two different symmetry classes for it's neighbours -> new stereocenter
     *   - has two identical symmetry classes, but the ligand contains at least 1 true or para stereocenter -> new stereocenter
     *
     * NOTE: there must always be at least 2 new stereocenters (or one existing + 1 newly found) in order for them to be valid
     */
    std::vector<impl::StereoRing<MoleculeType> > stereoRings;

    //cout << "=====================================================" << endl;
    for (std::size_t i = 0; i < rings.size(); ++i) {
      stereoRings.push_back(impl::StereoRing<MoleculeType>());
      impl::StereoRing<MoleculeType> &ring = stereoRings.back();


      for (std::size_t j = 0; j < stereoAtoms.size(); ++j)
        if (rings.ring(i).containsAtom(get_atom(mol, stereoAtoms[j])))
          ring.trueCount++;

      //cout << "StereoRing: trueCount = " << ring.trueCount << endl;
      for (std::size_t j = 0; j < paraAtoms.size(); ++j) {
        if (rings.ring(i).containsAtom(get_atom(mol, paraAtoms[j]))) {
          atom_type atom = get_atom(mol, paraAtoms[j]);
          ring.paraAtoms.push_back(impl::ParaAtom<MoleculeType>(paraAtoms[j]));

          FOREACH_NBR_T (nbr, atom, mol, MoleculeType) {
            if (rings.ring(i).containsAtom(*nbr))
              ring.paraAtoms.back().insideNbrs.push_back(*nbr);
            else
              ring.paraAtoms.back().outsideNbrs.push_back(*nbr);
          }

          //cout << "  ParaAtom(idx = " << ring.paraAtoms.back().inIdx << ", outside = " << ring.paraAtoms.back().outsideNbrs.size() << ")" << endl;
          if (ring.paraAtoms.back().insideNbrs.size() != 2)
            ring.paraAtoms.pop_back();
        }
      }

      for (std::size_t j = 0; j < paraBonds.size(); ++j) {
        bond_type bond = get_bond(mol, paraBonds[j]);
        Index beginIdx = get_index(mol, get_source(mol, bond));
        Index endIdx = get_index(mol, get_target(mol, bond));

        if (rings.ring(i).containsAtom(get_atom(mol, beginIdx))) {
          ring.paraBonds.push_back(impl::ParaBond<MoleculeType>(get_index(mol, bond), beginIdx, endIdx));

          FOREACH_NBR_T (nbr, get_source(mol, bond), mol, MoleculeType) {
            if (get_index(mol, *nbr) == endIdx)
              continue;
            ring.paraBonds.back().insideNbrs.push_back(*nbr);
          }
          FOREACH_NBR_T (nbr, get_target(mol, bond), mol, MoleculeType) {
            if (get_index(mol, *nbr) == beginIdx)
              continue;
            ring.paraBonds.back().outsideNbrs.push_back(*nbr);
          }

          //cout << "  ParaBond(inIdx = " << beginIdx << ", outIdx = " << endIdx << ", outside = " << ring.paraBonds.back().outsideNbrs.size() << ")" << endl;
          if (ring.paraBonds.back().insideNbrs.size() != 2)
            ring.paraBonds.pop_back();
        }

        if (rings.ring(i).containsAtom(get_atom(mol, endIdx))) {
          ring.paraBonds.push_back(impl::ParaBond<MoleculeType>(get_index(mol, bond), endIdx, beginIdx));

          FOREACH_NBR_T (nbr, get_target(mol, bond), mol, MoleculeType) {
            if (get_index(mol, *nbr) == beginIdx)
              continue;
            ring.paraBonds.back().insideNbrs.push_back(*nbr);
          }
          FOREACH_NBR_T (nbr, get_source(mol, bond), mol, MoleculeType) {
            if (get_index(mol, *nbr) == endIdx)
              continue;
            ring.paraBonds.back().outsideNbrs.push_back(*nbr);
          }

          //cout << "  ParaBond(inIdx = " << endIdx << ", outIdx = " << beginIdx << ", outside = " << ring.paraBonds.back().outsideNbrs.size() << ")" << endl;
          if (ring.paraBonds.back().insideNbrs.size() != 2)
            ring.paraBonds.pop_back();
        }

      }

      if (ring.paraAtoms.size() + ring.paraBonds.size() == 1) {
        ring.paraAtoms.clear();
        ring.paraBonds.clear();
      }

    }
    //cout << "=====================================================" << endl;

    unsigned int numStereoUnits;
    do {
      numStereoUnits = units.size();
      impl::start_rule1(mol, ec, stereoRings, units, stereoAtoms);
    } while (units.size() > numStereoUnits);


    // XXXXXXXXXXXXXXXXXXXX
    std::vector<std::vector<bool> > mergedRings = impl::merge_rings(mol, rings, ec);
    /**
     * Apply rule 2a for tetracoordinate carbon:
     * - 1 or 2 pair identical ligands
     * - each pair contains at least 1 true-stereocenter or 2 para-stereocenters
     *
     * Apply rule 2b for tetracoordinate carbon:
     * - 3 or 4 identical ligands with at least
     *   - 2 true-stereocenters
     *   - 2 separate assemblies of para-stereocenters
     */
    for (std::vector<Index>::iterator idx = paraAtoms.begin(); idx != paraAtoms.end(); ++idx) {
      atom_type atom = get_atom(mol, *idx);
      // make sure we didn't add this atom already from rule 1
      bool alreadyAdded = false;
      for (StereoUnitSet::iterator u2 = units.begin(); u2 != units.end(); ++u2) {
        if ((*u2).type == Stereo::Tetrahedral)
          if (get_index(mol, atom) == (*u2).index) {
            alreadyAdded = true;
          }
      }
      if (alreadyAdded)
        continue;

      int classification = impl::classify_tetrahedral_nbr_symmetry_classes(mol, atom, ec);
      switch (classification) {
        case impl::T1123:
          // rule 2a with 1 pair
          {
            unsigned long duplicatedSymClass = impl::find_duplicated_symmetry_class(mol, atom, ec);
            atom_type ligandAtom = impl::find_atom_with_symmetry_class(mol, atom, duplicatedSymClass, ec);
            if (impl::contains_at_least_1true_2para(mol, ligandAtom, atom, units, rings))
              units.push_back(StereoUnit(Stereo::Tetrahedral, get_index(mol, atom), true));
          }
          break;
        case impl::T1122:
          // rule 2a with 2 pairs
          {
            unsigned long duplicatedSymClass1, duplicatedSymClass2;
            impl::find_duplicated_symmetry_classes(mol, atom, ec, duplicatedSymClass1, duplicatedSymClass2);
            atom_type ligandAtom1 = impl::find_atom_with_symmetry_class(mol, atom, duplicatedSymClass1, ec);
            atom_type ligandAtom2 = impl::find_atom_with_symmetry_class(mol, atom, duplicatedSymClass2, ec);
            if (impl::contains_at_least_1true_2para(mol, ligandAtom1, atom, units, rings) &&
                impl::contains_at_least_1true_2para(mol, ligandAtom2, atom, units, rings))
              units.push_back(StereoUnit(Stereo::Tetrahedral, get_index(mol, atom), true));
          }
          break;
        case impl::T1112:
          // rule 2b with 3 identical
          {
            unsigned long duplicatedSymClass = impl::find_duplicated_symmetry_class(mol, atom, ec);
            atom_type ligandAtom = impl::find_atom_with_symmetry_class(mol, atom, duplicatedSymClass, ec);
            if (impl::contains_at_least_2true_2paraAssemblies(mol, ligandAtom, atom, units, mergedRings))
              units.push_back(StereoUnit(Stereo::Tetrahedral, get_index(mol, atom), true));
          }
          break;
        case impl::T1111:
          // rule 2b with 4 identical
          {
            unsigned long duplicatedSymClass = impl::find_duplicated_symmetry_class(mol, atom, ec);
            atom_type ligandAtom = impl::find_atom_with_symmetry_class(mol, atom, duplicatedSymClass, ec);
            if (impl::contains_at_least_2true_2paraAssemblies(mol, ligandAtom, atom, units, mergedRings))
              units.push_back(StereoUnit(Stereo::Tetrahedral, get_index(mol, atom), true));
          }
          break;

      }

    }

    /**
     * Apply rule 3 for double bonds.
     * - 1 or 2 pair identical ligands (on begin and end atom)
     * - each pair contains at least 1 true-stereocenter or 2 para-stereocenters (from rule1)
     */
    for (std::vector<Index>::iterator idx = paraBonds.begin(); idx != paraBonds.end(); ++idx) {
      bond_type bond = get_bond(mol, *idx);

      // make sure we didn't add this atom already from rule 1
      bool alreadyAdded = false;
      for (StereoUnitSet::iterator unit = units.begin(); unit != units.end(); ++unit) {
        if (unit->type == Stereo::CisTrans)
          if (get_index(mol, bond) == unit->index) {
            alreadyAdded = true;
          }
      }
      if (alreadyAdded)
        continue;

      atom_type begin = get_source(mol, bond);
      atom_type end = get_target(mol, bond);

      int beginClassification = impl::classify_cis_trans_nbr_symmetry_classes(mol, ec, bond, begin);
      bool beginValid = false;
      switch (beginClassification) {
        case impl::C12:
          beginValid = true;
          break;
        case impl::C11:
          {
            // find the ligand
            atom_type ligandAtom = molecule_traits<MoleculeType>::null_atom();
            FOREACH_NBR_T (nbr, begin, mol, MoleculeType) {
              if ((get_index(mol, *nbr) != get_index(mol, begin)) && (get_index(mol, *nbr) != get_index(mol, end))) {
                ligandAtom = *nbr;
                break;
              }
            }

            std::vector<bool> ligand = impl::get_fragment(mol, ligandAtom, begin);
            for (StereoUnitSet::iterator unit = units.begin(); unit != units.end(); ++unit) {
              if (unit->type == Stereo::Tetrahedral) {
                if (ligand[unit->index])
                  beginValid = true;
              } else if(unit->type == Stereo::CisTrans) {
                bond_type bond = get_bond(mol, unit->index);
                atom_type begin = get_source(mol, bond);
                atom_type end = get_target(mol, bond);
                if (ligand[get_index(mol, begin)] || ligand[get_index(mol, end)])
                  beginValid = true;
              }
            }
          }
          break;
      }

      if (!beginValid)
        continue;

      int endClassification = impl::classify_cis_trans_nbr_symmetry_classes(mol, ec, bond, end);
      bool endValid = false;
      switch (endClassification) {
        case impl::C12:
          endValid = true;
          break;
        case impl::C11:
          {
            // find the ligand
            atom_type ligandAtom = molecule_traits<MoleculeType>::null_atom();
            FOREACH_NBR_T (nbr, end, mol, MoleculeType) {
              if ((get_index(mol, *nbr) != get_index(mol, begin)) && (get_index(mol, *nbr) != get_index(mol, end))) {
                ligandAtom = *nbr;
                break;
              }
            }

            std::vector<bool> ligand = impl::get_fragment(mol, ligandAtom, end);
            for (StereoUnitSet::iterator unit = units.begin(); unit != units.end(); ++unit) {
              if (unit->type == Stereo::Tetrahedral) {
                if (ligand[unit->index])
                  endValid = true;
              } else if(unit->type == Stereo::CisTrans) {
                bond_type bond = get_bond(mol, unit->index);
                atom_type begin = get_source(mol, bond);
                atom_type end = get_target(mol, bond);
                if (ligand[get_index(mol, begin)] || ligand[get_index(mol, end)])
                  endValid = true;
              }
            }
          }
          break;
      }

      if (endValid)
        units.push_back(StereoUnit(Stereo::CisTrans, get_index(mol, bond), true));
    }

    if (DEBUG_RAZINGER) {
      for (StereoUnitSet::iterator unit = units.begin(); unit != units.end(); ++unit) {
        if (unit->type == Stereo::Tetrahedral)
          std::cout << "Tetrahedral(center = " << unit->index << ", para = " << unit->para << ")" << std::endl;
        if (unit->type == Stereo::CisTrans)
          std::cout << "CisTrans(bond = " << unit->index << ", para = " << unit->para << ")" << std::endl;
        if (unit->type == Stereo::SquarePlanar)
          std::cout << "SquarePlanar(bond = " << unit->index << ", para = " << unit->para << ")" << std::endl;
      }
    }

    return units;
  }

}

#endif
