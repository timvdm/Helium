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
#ifndef HELIUM_SMIRKS_H
#define HELIUM_SMIRKS_H

#include <Helium/smarts.h>

#include <iostream>

namespace Helium {

  /**
   * @brief Class representing possible SMIRKS errors.
   */
  class SmirksError
  {
    public:
      /**
       * @brief The error type.
       */
      enum Type {
        /**
         * @brief No error.
         */
        None,
        /**
         * @brief The SMIRKS does not contain '>>'.
         */
        NoReaction,
        /**
         * @brief Parse error in reactant SMARTS.
         */
        ReactantSmarts,
        /**
         * @brief Parse error in product SMARTS.
         */
        ProductSmarts,
        /**
         * @brief The atom classes are not pair-wise.
         */
        AtomClassPairWise,
        /**
         * @brief The product SMARTS contains an OR expression.
         */
        ProductContainsOr,
        /**
         * @brief The product SMARTS contains a NOT expression.
         */
        ProductContainsNot,
        /**
         * @brief The product contains a complex bond expression.
         */
        InvalidProductBond
      };

      /**
       * @brief Default constructor.
       *
       * This constructor initializes the error to None.
       */
      SmirksError() : m_type(None)
      {
      }

      /**
       * @brief Constructor.
       *
       * @param type The error type.
       * @param what The error message.
       */
      SmirksError(Type type, const std::string &what) : m_type(type), m_what(what)
      {
      }

      /**
       * @brief Get the error type.
       *
       * @return The error type.
       */
      Type type() const
      {
        return m_type;
      }

      /**
       * @brief Get the error message.
       *
       * @return The error message.
       */
      const std::string& what() const
      {
        return m_what;
      }

    private:
      Type m_type; //!< The error type.
      std::string m_what; //!< The error message.
  };

  /**
   * @brief Class for applying SMIRKS transformations.
   */
  class Smirks
  {
      /**
       * @brief Class representing a change to a bond.
       */
      struct BondChange
      {
        /**
         * @brief The type of change.
         */
        enum Type {
          /**
           * @brief Bond is changed.
           */
          Changed,
          /**
           * @brief Bond is added.
           */
          Added,
          /**
           * @brief Bond is removed.
           */
          Removed
        };

        /**
         * @brief Constructor.
         *
         * @param type_ The error type.
         * @param source_ The reactant source index.
         * @param target_ The reactant target index.
         * @param expr_ The product SMARTS bond expression.
         */
        BondChange(int type_, int source_, int target_, impl::SmartsBondExpr *expr_)
          : type(type_), source(source_), target(target_), expr(expr_)
        {
        }

        int type; //!< The type of change.
        int source; //!< Reactant source index.
        int target; //!< Reactant target index.
        impl::SmartsBondExpr *expr; //!< Product SMARTS bond expression.
      };

    public:
      /**
       * @brief Initialize the SMIRKS.
       *
       * @param smirks The SMIRKS string.
       *
       * @return True if successful.
       */
      bool init(const std::string &smirks);

      /**
       * @brief Initialize the SMIRKS.
       *
       * @param reactant The reactant SMARTS string.
       * @param product The product SMARTS string.
       *
       * @return True if successful.
       */
      bool init(const std::string &reactant, const std::string &product);

      /**
       * @brief Get the error from the last call to init().
       *
       * @return The SMIRKS erorr.
       */
      const SmirksError& error() const
      {
        return m_error;
      }

      /**
       * @brief Apply the SMIRKS transformation to a molecule.
       *
       * @param mol The molecule.
       * @param rings The molecule's rings (needed for cyclic queries).
       *
       * @return True if changes were made to the molecule.
       */
      template<typename EditableMoleculeType>
      bool apply(EditableMoleculeType &mol, const RingSet<EditableMoleculeType> &rings)
      {
        MappingList mapping;
        if (!m_reactant.search(mol, mapping, rings, true))
          return false;

        // apply atom changes
        for (std::size_t i = 0; i < mapping.maps.size(); ++i) {
          const IsomorphismMapping &map = mapping.maps[i];
          for (std::size_t j = 0; j < map.size(); ++j) {
            int atomClass = m_reactant.atomClass(j);
            if (atomClass == -1)
              continue;

            impl::SmartsAtomExpr *productExpr = m_productExpr[atomClass];

            apply(mol, get_atom(mol, map[j]), productExpr);
          }
        }

        // apply bond changes
        std::vector<Index> removeBonds;
        for (std::size_t i = 0; i < mapping.maps.size(); ++i) {
          const IsomorphismMapping &map = mapping.maps[i];

          for (std::size_t j = 0; j < m_bondChanges.size(); ++j) {
            const BondChange &change = m_bondChanges[j];

            switch (change.type) {
              case BondChange::Changed:
                {
                  molecule_traits<HeMol>::atom_type source = get_atom(mol, map[change.source]);
                  molecule_traits<HeMol>::atom_type target = get_atom(mol, map[change.target]);
                  molecule_traits<HeMol>::bond_type bond = get_bond(mol, source, target);
                  apply(mol, bond, change.expr);
                }
                break;
              case BondChange::Removed:
                {
                  molecule_traits<HeMol>::atom_type source = get_atom(mol, map[change.source]);
                  molecule_traits<HeMol>::atom_type target = get_atom(mol, map[change.target]);
                  molecule_traits<HeMol>::bond_type bond = get_bond(mol, source, target);
                  removeBonds.push_back(get_index(mol, bond));
                }
                break;
              case BondChange::Added:
                {
                  molecule_traits<HeMol>::atom_type source = get_atom(mol, map[change.source]);
                  molecule_traits<HeMol>::atom_type target = get_atom(mol, map[change.target]);
                  molecule_traits<HeMol>::bond_type bond = add_bond(mol, source, target);
                  apply(mol, bond, change.expr);
                }
                break;
            }
          }
        }

        // remove the planned bonds
        std::sort(removeBonds.begin(), removeBonds.end(), std::greater<Index>());
        for (std::size_t i = 0; i < removeBonds.size(); ++i)
          remove_bond(mol, get_bond(mol, removeBonds[i]));

        return true;
      }

    private:
      /**
       * @brief Apply changes to an atom.
       *
       * @param mol The molecule.
       * @param atom The atom to change.
       * @param expr The product SMARTS atom expression.
       */
      template<typename EditableMoleculeType, typename AtomType>
      void apply(EditableMoleculeType &mol, AtomType atom, impl::SmartsAtomExpr *expr)
      {
        switch (expr->type) {
          case Smiley::OP_AndHi:
          case Smiley::OP_AndLo:
          case Smiley::OP_And:
          case Smiley::OP_Or:
            apply(mol, atom, expr->left);
            apply(mol, atom, expr->right);
            break;
          case Smiley::AE_Isotope:
            set_mass(mol, atom, expr->value);
            break;
          case Smiley::AE_AtomicNumber:
            set_element(mol, atom, expr->value);
            break;
          case Smiley::AE_AromaticElement:
            set_element(mol, atom, expr->value);
            set_aromatic(mol, atom, true);
            break;
          case Smiley::AE_AliphaticElement:
            set_element(mol, atom, expr->value);
            set_aromatic(mol, atom, false);
            break;
          case Smiley::AE_TotalH:
            {
              int explicitH = 0;
              FOREACH_NBR (nbr, atom, mol, EditableMoleculeType)
                if (is_hydrogen(mol, *nbr))
                  ++explicitH;
              set_hydrogens(mol, atom, expr->value - explicitH);
            }
            break;
          case Smiley::AE_ImplicitH:
            set_hydrogens(mol, atom, expr->value);
            break;
          case Smiley::AE_Charge:
            set_charge(mol, atom, expr->value);
            break;
          case Smiley::AE_Chirality:
          case Smiley::OP_Not:
          default:
            break;
        }
      }

      /**
       * @brief Apply changes to a bond.
       *
       * @param mol The molecule.
       * @param bond The bond to change.
       * @param expr The product SMARTS bond expression.
       */
      template<typename EditableMoleculeType, typename BondType>
      void apply(EditableMoleculeType &mol, BondType bond, impl::SmartsBondExpr *expr)
      {
        switch (expr->type) {
          case Smiley::BE_Single:
            set_order(mol, bond, 1);
            set_aromatic(mol, bond, false);
            break;
          case Smiley::BE_Double:
            set_order(mol, bond, 2);
            set_aromatic(mol, bond, false);
            break;
          case Smiley::BE_Triple:
            set_order(mol, bond, 3);
            set_aromatic(mol, bond, false);
            break;
          case Smiley::BE_Quadriple:
            set_order(mol, bond, 4);
            set_aromatic(mol, bond, false);
            break;
          case Smiley::BE_Aromatic:
            set_order(mol, bond, 5);
            set_aromatic(mol, bond, true);
            break;
          case Smiley::BE_Up:
          case Smiley::BE_Down:
            break;
        }
      }

      /**
       * @brief Extract atom class to SMARTS expression mapping from SMARTS object.
       *
       * @param smarts The SMARTS object.
       *
       * @return The atom class to SMARTS expression mapping.
       */
      std::map<int, impl::SmartsAtomExpr*> atomClassToExpr(const Smarts &smarts) const;

      /**
       * @brief Find a bond in a SMARTS between the specified atom classes.
       *
       * @param smarts The SMARTS object.
       * @param sourceAtomClass The source atom class.
       * @param targetAtomClass The target atom class.
       *
       * @return The found bond or molecule_traits<HeMol>::null_bond().
       */
      molecule_traits<HeMol>::bond_type getBond(const Smarts &smarts, int sourceAtomClass, int targetAtomClass);

      Smarts m_reactant; //!< The reactant SMARTS.
      Smarts m_product; //!< The product SMARTS.
      std::map<int, impl::SmartsAtomExpr*> m_reactantExpr; //!< Reactant atom class to atom expr.
      std::map<int, impl::SmartsAtomExpr*> m_productExpr; //!< Product atom class to atom expr.
      std::vector<BondChange> m_bondChanges; //!< The bond changes.
      SmirksError m_error; //!< Error from last call to init().
  };

}

#endif
