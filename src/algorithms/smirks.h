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
        ProductContainsNot
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
      struct AtomClass
      {
        AtomClass(int component_, int index_, int atomClass_)
          : component(component_), index(index_), atomClass(atomClass_)
        {
        }

        int component;
        int index;
        int atomClass;
      };

      struct AtomExprMap
      {
        AtomExprMap(impl::SmartsAtomExpr *reactant_ = 0, impl::SmartsAtomExpr *product_ = 0)
          : reactant(reactant_), product(product_)
        {
        }

        impl::SmartsAtomExpr *reactant;
        impl::SmartsAtomExpr *product;
      };

    public:
      bool init(const std::string &smirks)
      {
        std::size_t pos = smirks.find(">>");
        if (pos == std::string::npos) {
          m_error = SmirksError(SmirksError::NoReaction, "SMIRKS does not contain '>>'");
          return false;
        }

        return init(smirks.substr(0, pos), smirks.substr(pos + 2));
      }

      bool init(const std::string &reactant, const std::string &product)
      {
        //std::cout << "Smirks::init(" << reactant << ">>" << product << ")" << std::endl;

        m_error = SmirksError();

        if (!m_reactant.init(reactant)) {
          m_error = SmirksError(SmirksError::ReactantSmarts, m_reactant.error().what());
          return false;
        }
        if (!m_product.init(product)) {
          m_error = SmirksError(SmirksError::ReactantSmarts, m_product.error().what());
          return false;
        }

        m_reactantExpr = atomClassToExpr(m_reactant);
        m_productExpr = atomClassToExpr(m_product);

        // atom mapping must be pair-wise
        if (m_reactantExpr.size() != m_productExpr.size()) {
          m_error = SmirksError(SmirksError::AtomClassPairWise, "Atom class mapping must be pairwise");
          return false;
        }

        for (std::map<int, impl::SmartsAtomExpr*>::const_iterator i = m_reactantExpr.begin(); i != m_reactantExpr.end(); ++i)
          if (m_productExpr.find(i->first) == m_productExpr.end()) {
            m_error = SmirksError(SmirksError::AtomClassPairWise, "Atom class mapping must be pairwise");
            return false;
          }

        // product may not contain OR
        for (std::size_t i = 0; i < num_atoms(m_product.query()); ++i)
          if (exprContains(m_product.trees().atom(i), Smiley::OP_Or)) {
            m_error = SmirksError(SmirksError::ProductContainsOr, "Product may not contain OR expression");
            return false;
          }

        // product may not contain NOT
        for (std::size_t i = 0; i < num_atoms(m_product.query()); ++i)
          if (exprContains(m_product.trees().atom(i), Smiley::OP_Not)) {
            m_error = SmirksError(SmirksError::ProductContainsOr, "Product may not contain NOT expression");
            return false;
          }

        return true;
      }

      const SmirksError& error() const
      {
        return m_error;
      }

      template<typename EditableMoleculeType>
      bool apply(EditableMoleculeType &mol, const RingSet<EditableMoleculeType> &rings)
      {
        MappingList mapping;
        if (!m_reactant.search(mol, mapping, rings))
          return false;

        for (std::size_t i = 0; i < mapping.maps.size(); ++i) {
          const IsomorphismMapping &map = mapping.maps[i];
          for (std::size_t j = 0; j < map.size(); ++j) {
            int atomClass = m_reactant.atomClass(j);
            if (atomClass == -1)
              continue;

            //impl::SmartsAtomExpr *reactantExpr = m_reactantExpr[atomClass];
            impl::SmartsAtomExpr *productExpr = m_productExpr[atomClass];

            apply(mol, get_atom(mol, map[j]), productExpr);
          }
        }

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

      std::map<int, impl::SmartsAtomExpr*> atomClassToExpr(const Smarts &smarts) const
      {
        std::map<int, impl::SmartsAtomExpr*> result;

        for (std::size_t i = 0; i < num_atoms(smarts.query()); ++i) {
          int atomClass = smarts.atomClass(i);
          if (atomClass != -1)
            result[atomClass] = smarts.trees().atom(i);
        }

        return result;
      }


      Smarts m_reactant;
      Smarts m_product;
      std::map<int, impl::SmartsAtomExpr*> m_reactantExpr;
      std::map<int, impl::SmartsAtomExpr*> m_productExpr;
      SmirksError m_error;
  };

}

#endif
