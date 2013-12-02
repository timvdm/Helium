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

#include <Helium/algorithms/smarts.h>

#include <iostream>

namespace Helium {

  class SmirksError
  {
    public:
      enum Type {
        None,
        NoReaction,
        ReactantSmarts,
        ProductSmarts,
        AtomClassPairWise,
        ProductContainsOr,
        ProductContainsNot,
        InvalidProductBond
      };

      SmirksError() : m_type(None)
      {
      }

      SmirksError(Type type, const std::string &what) : m_type(type), m_what(what)
      {
      }

      Type type() const
      {
        return m_type;
      }

      const std::string& what() const
      {
        return m_what;
      }

    private:
      Type m_type;
      std::string m_what;
  };

  class Smirks
  {
      // bond exists in both reactant and product
      struct BondChange
      {
        enum Type {
          Changed,
          Added,
          Removed
        };

        BondChange(int type_, int source_, int target_, impl::SmartsBondExpr *expr_)
          : type(type_), source(source_), target(target_), expr(expr_)
        {
        }

        int type;
        int source; // reactant source index
        int target; // reactant target index
        impl::SmartsBondExpr *expr; // product bond expr
      };

    public:
      bool init(const std::string &smirks);

      bool init(const std::string &reactant, const std::string &product);

      const SmirksError& error() const
      {
        return m_error;
      }

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

      template<typename ExprType>
      bool exprContains(ExprType *expr, int type)
      {
        if (expr->type == type)
          return true;

        switch (expr->type) {
          case Smiley::OP_AndHi:
          case Smiley::OP_AndLo:
          case Smiley::OP_And:
          case Smiley::OP_Or:
            if (exprContains(expr->left, type))
              return true;
            if (exprContains(expr->right, type))
              return true;
            break;
          case Smiley::OP_Not:
            if (exprContains(expr->arg, type))
              return true;
            break;
          default:
            break;
        }

        return false;
      }

      std::map<int, impl::SmartsAtomExpr*> atomClassToExpr(const Smarts &smarts) const;

      molecule_traits<HeMol>::bond_type getBond(const Smarts &smarts, int sourceAtomClass, int targetAtomClass);

      Smarts m_reactant;
      Smarts m_product;
      std::map<int, impl::SmartsAtomExpr*> m_reactantExpr;
      std::map<int, impl::SmartsAtomExpr*> m_productExpr;
      std::vector<BondChange> m_bondChanges;
      SmirksError m_error;
  };

}

#endif
