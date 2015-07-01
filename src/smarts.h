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
#ifndef HELIUM_SMARTS_H
#define HELIUM_SMARTS_H

#include <Helium/config.h>
#include <Helium/hemol.h>
#include <Helium/ring.h>
#include <Helium/algorithms/isomorphism.h>
#include <Helium/algorithms/components.h>
#include <Helium/algorithms/cycles.h>
#include <Helium/error.h>
#include <Helium/smiley.h>

#include <iostream>

namespace Helium {

  /**
   * @file smarts.h
   * @brief SMARTS substructure search.
   */

  namespace impl {

    struct SmartsAtomExpr
    {
      SmartsAtomExpr(int type_)
        : type(type_)
      {
      }
      SmartsAtomExpr(int type_, int value_)
        : type(type_), value(value_)
      {
      }
      SmartsAtomExpr(int type_, SmartsAtomExpr *arg_)
        : type(type_), arg(arg_)
      {
      }
      SmartsAtomExpr(int type_, SmartsAtomExpr *left_, SmartsAtomExpr *right_)
        : type(type_), left(left_), right(right_)
      {
      }

      int type;
      union {
        int value;
        SmartsAtomExpr *arg;
        SmartsAtomExpr *left;
      };
      SmartsAtomExpr *right;
    };

    struct SmartsBondExpr
    {
      SmartsBondExpr(int type_)
        : type(type_)
      {
      }
      SmartsBondExpr(int type_, SmartsBondExpr *arg_)
        : type(type_), arg(arg_)
      {
      }
      SmartsBondExpr(int type_, SmartsBondExpr *left_, SmartsBondExpr *right_)
        : type(type_), left(left_), right(right_)
      {
      }


      int type;
      union {
        int value;
        SmartsBondExpr *arg;
        SmartsBondExpr *left;
      };
      SmartsBondExpr *right;
    };

    class SmartsTrees
    {
      public:
        SmartsTrees();

        SmartsTrees(const SmartsTrees &other);

        ~SmartsTrees();

        SmartsTrees& operator=(const SmartsTrees &other);

        void clear();

        void addAtom(SmartsAtomExpr *expr)
        {
          m_atoms.push_back(expr);
        }

        SmartsAtomExpr* atom(std::size_t index) const
        {
          return m_atoms[index];
        }

        const std::vector<SmartsAtomExpr*>& atoms() const
        {
          return m_atoms;
        }

        void addBond(SmartsBondExpr *expr)
        {
          m_bonds.push_back(expr);
        }

        SmartsBondExpr* bond(std::size_t index) const
        {
          return m_bonds[index];
        }

        const std::vector<SmartsBondExpr*>& bonds() const
        {
          return m_bonds;
        }

        SmartsAtomExpr* copy(SmartsAtomExpr *expr);

        SmartsBondExpr* copy(SmartsBondExpr *expr);

      private:
        template<typename ExprType>
        void cleanup(ExprType *expr)
        {
          switch (expr->type) {
            case Smiley::OP_AndHi:
            case Smiley::OP_AndLo:
            case Smiley::OP_And:
            case Smiley::OP_Or:
              cleanup(expr->left);
              cleanup(expr->right);
              break;
            case Smiley::OP_Not:
              cleanup(expr->arg);
              break;
            default:
              break;
          }
          delete expr;
        }

        std::vector<SmartsAtomExpr*> m_atoms;
        std::vector<SmartsBondExpr*> m_bonds;
    };

    /**
     * @brief Check if an expression tree contains a specific type.
     *
     * @param expr The SMARTS expression tree.
     * @param type The type to search for.
     *
     * @return True if the @p type is found in the expression tree.
     */
    template<typename ExprType>
    bool smarts_expr_contains(ExprType *expr, int type)
    {
      if (expr->type == type)
        return true;

      switch (expr->type) {
        case Smiley::OP_AndHi:
        case Smiley::OP_AndLo:
        case Smiley::OP_And:
        case Smiley::OP_Or:
          if (smarts_expr_contains(expr->left, type))
            return true;
          if (smarts_expr_contains(expr->right, type))
            return true;
          break;
        case Smiley::OP_Not:
          if (smarts_expr_contains(expr->arg, type))
            return true;
          break;
        default:
          break;
      }

      return false;
    }

    /**
     * @brief Check if an expression tree contains a specific type and value.
     *
     * @param expr The SMARTS expression tree.
     * @param type The type to search for.
     * @param value The value to search for.
     *
     * @return True if the @p type and @p value is found in the expression tree.
     */
    template<typename ExprType>
    bool smarts_expr_contains(ExprType *expr, int type, int value)
    {
      if (expr->type == type && expr->value == value)
        return true;

      switch (expr->type) {
        case Smiley::OP_AndHi:
        case Smiley::OP_AndLo:
        case Smiley::OP_And:
        case Smiley::OP_Or:
          if (smarts_expr_contains(expr->left, type, value))
            return true;
          if (smarts_expr_contains(expr->right, type, value))
            return true;
          break;
        case Smiley::OP_Not:
          if (smarts_expr_contains(expr->arg, type, value))
            return false;
          break;
        default:
          break;
      }

      return false;
    }
    /**
     * @brief Change all occurences of type and value with the new type and value.
     *
     * @param expr The SMARTS expression tree.
     * @param type The type to search for.
     * @param value The value to search for.
     * @param newType The new type.
     * @param newValue The new value.
     */
    template<typename ExprType>
    void smarts_expr_change(ExprType *expr, int type, int value, int newType, int newValue)
    {
      if (expr->type == type && expr->value == value) {
        expr->type = newType;
        expr->value = newValue;
        return;
      }

      switch (expr->type) {
        case Smiley::OP_AndHi:
        case Smiley::OP_AndLo:
        case Smiley::OP_And:
        case Smiley::OP_Or:
          smarts_expr_change(expr->left, type, value, newType, newValue);
          smarts_expr_change(expr->right, type, value, newType, newValue);
          break;
        case Smiley::OP_Not:
          smarts_expr_change(expr->arg, type, value, newType, newValue);
          break;
        default:
          break;
      }
    }

    template<typename QueryType, typename MoleculeType>
    class SmartsBondMatcher;

    template<typename QueryType, typename MoleculeType>
    class SmartsAtomMatcher
    {
      public:
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<QueryType>::atom_type query_atom_type;

        SmartsAtomMatcher(const std::vector<SmartsAtomExpr*> &atoms, const RingSet<MoleculeType> &rings,
            const std::vector<HeMol> &recursiveMols, const std::vector<SmartsTrees> &recursiveTrees, int mode)
          : m_atoms(atoms), m_rings(rings), m_recursiveMols(recursiveMols),
            m_recursiveTrees(recursiveTrees), m_mode(mode)
        {
        }

        bool operator()(const QueryType &query, query_atom_type queryAtom,
            const MoleculeType &mol, atom_type atom) const
        {
          SmartsAtomExpr *expr = m_atoms[get_index(query, queryAtom)];
          return match(mol, atom, expr);
        }

      private:
        bool match(const MoleculeType &mol, atom_type atom, SmartsAtomExpr *expr) const
        {
          switch (expr->type) {
            case Smiley::AE_True:
              return true;
            case Smiley::AE_False:
              return false;
            case Smiley::AE_Aromatic:
              return is_aromatic(mol, atom);
            case Smiley::AE_Aliphatic:
              return !is_aromatic(mol, atom);
            case Smiley::AE_Cyclic:
              return m_rings.isAtomInRing(atom);
            case Smiley::AE_Acyclic:
              return !m_rings.isAtomInRing(atom);
            case Smiley::AE_Isotope:
              return get_mass(mol, atom) == expr->value;
            case Smiley::AE_AtomicNumber:
              return get_element(mol, atom) == expr->value;
            case Smiley::AE_AromaticElement:
              return get_element(mol, atom) == expr->value && is_aromatic(mol, atom);
            case Smiley::AE_AliphaticElement:
              return get_element(mol, atom) == expr->value && !is_aromatic(mol, atom);
            case Smiley::AE_Degree:
              return get_degree(mol, atom) == expr->value;
            case Smiley::AE_Valence:
              return get_valence(mol, atom) == expr->value;
            case Smiley::AE_Connectivity:
              return get_connectivity(mol, atom) == expr->value;
            case Smiley::AE_TotalH:
              {
                int h = 0;
                for (auto &nbr : get_nbrs(mol, atom))
                  if (get_element(mol, nbr) == 1)
                    ++h;
                return (h + get_hydrogens(mol, atom)) == expr->value;
              }
            case Smiley::AE_ImplicitH:
              if (expr->value == -1) // default: at least 1
                return get_hydrogens(mol, atom) >= 1;
              else
                return get_hydrogens(mol, atom) == expr->value;
            case Smiley::AE_RingMembership:
              if (m_mode == 3 /*Smarts::OpenEye*/) {
                // R<n> same as x<n>
                if (expr->value == -1) // default: at least 1
                  return m_rings.numRingNbrs(atom) >= 1;
                else
                  return m_rings.numRingNbrs(atom) == expr->value;
              } else {
                return m_rings.numRings(atom) == expr->value;
              }
            case Smiley::AE_RingSize:
              return m_rings.isAtomInRingSize(atom, expr->value);
            case Smiley::AE_RingConnectivity:
              if (expr->value == -1) // default: at least 1
                return m_rings.numRingNbrs(atom) >= 1;
              else
                return m_rings.numRingNbrs(atom) == expr->value;
            case Smiley::AE_Charge:
              return get_charge(mol, atom) == expr->value;
            case Smiley::AE_Chirality:
              return true;
            case Smiley::AE_AtomClass:
              return true;
            case Smiley::AE_Recursive:
              {
                assert(expr->value < m_recursiveMols.size());
                assert(expr->value < m_recursiveTrees.size());
                impl::SmartsAtomMatcher<HeMol, MoleculeType> atomMatcher(m_recursiveTrees[expr->value].atoms(),
                    m_rings, m_recursiveMols, m_recursiveTrees, m_mode);
                impl::SmartsBondMatcher<HeMol, MoleculeType> bondMatcher(m_recursiveTrees[expr->value].bonds(), m_rings);
                NoMapping mapping;
                return isomorphism_search(mol, atom, m_recursiveMols[expr->value], mapping, atomMatcher, bondMatcher);
              }
            case Smiley::OP_Not:
              return !match(mol, atom, expr->arg);
            case Smiley::OP_AndHi:
            case Smiley::OP_AndLo:
            case Smiley::OP_And:
              return match(mol, atom, expr->left) && match(mol, atom, expr->right);
            case Smiley::OP_Or:
              return match(mol, atom, expr->left) || match(mol, atom, expr->right);
            default:
              return true;
          }
        }

        const std::vector<SmartsAtomExpr*> m_atoms;
        const RingSet<MoleculeType> &m_rings;
        const std::vector<HeMol> &m_recursiveMols;
        const std::vector<SmartsTrees> &m_recursiveTrees;
        int m_mode;
    };

    template<typename QueryType, typename MoleculeType>
    class SmartsBondMatcher
    {
      public:
        typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
        typedef typename molecule_traits<QueryType>::bond_type query_bond_type;

        SmartsBondMatcher(const std::vector<SmartsBondExpr*> &bonds, const RingSet<MoleculeType> &rings)
          : m_bonds(bonds), m_rings(rings)
        {
        }

        bool operator()(const QueryType &query, query_bond_type queryBond,
            const MoleculeType &mol, bond_type bond) const
        {
          SmartsBondExpr *expr = m_bonds[get_index(query, queryBond)];
          return match(mol, bond, expr);
        }

      private:
        bool match(const MoleculeType &mol, bond_type bond, SmartsBondExpr *expr) const
        {
          switch (expr->type) {
            case Smiley::BE_True:
              return true;
            case Smiley::BE_False:
              return false;
            case Smiley::BE_Single:
              return get_order(mol, bond) == 1;
            case Smiley::BE_Double:
              return get_order(mol, bond) == 2;
            case Smiley::BE_Triple:
              return get_order(mol, bond) == 3;
            case Smiley::BE_Quadriple:
              return get_order(mol, bond) == 4;
            case Smiley::BE_Aromatic:
              return is_aromatic(mol, bond);
            case Smiley::BE_Up:
            case Smiley::BE_Down:
              return true;
            case Smiley::BE_Ring:
              return m_rings.isBondInRing(bond);
            case Smiley::OP_Not:
              return !match(mol, bond, expr->arg);
            case Smiley::OP_AndHi:
            case Smiley::OP_AndLo:
            case Smiley::OP_And:
              return match(mol, bond, expr->left) && match(mol, bond, expr->right);
            case Smiley::OP_Or:
              return match(mol, bond, expr->left) || match(mol, bond, expr->right);
            default:
              return true;
          }
        }

        const std::vector<SmartsBondExpr*> m_bonds;
        const RingSet<MoleculeType> &m_rings;
    };

  }

  /**
   * @class Smarts smarts.h <Helium/smarts.h>
   * @brief Class for matching SMARTS.
   */
  class Smarts
  {
    public:
      /**
       * @brief SMARTS cycle matching modes.
       */
      enum Modes {
        /**
         * @brief Use the SSSR ring set for matching the R<n> atom primitive.
         */
        SSSR = 1,
        /**
         * @brief Use the relevant cycles ring set for matching the R<n> atom primitive.
         */
        RelevantCycles = 2,
        /**
         * @brief Use OpenEye R<n> semmantics (i.e. R<n> is the same as x<n>).
         */
        OpenEye = 3
      };

      /**
       * @brief Initialize with the specified SMARTS.
       *
       * The default cycle match mode is OpenEye (i.e. R<n> is the same as x<n>).
       * Setting the mode to SSSR or RelevantCycles will only have effect when the
       * overloaded find functions are used without explicit ring set or cyclic atoms
       * parameters in which case the correct ring set is computed when necessary.
       *
       * @param smarts The SMARTS string.
       * @param mode The cycle matching mode (i.e. Smarts::Modes).
       *
       * @return True when succesfull. Use error() when false is returned.
       */
      bool init(const std::string &smarts, int mode = OpenEye);

      /**
       * @brief Get the query molecule.
       *
       * @return The query molecule.
       */
      const HeMol& query() const
      {
        return m_query;
      }

      /**
       * @brief Get the SMARTS expression trees.
       *
       * @return The SMARTS expression trees.
       */
      const impl::SmartsTrees& trees() const
      {
        return m_trees;
      }

      /**
       * @brief Get the list of recursive molecules.
       *
       * These are the molecules for the recursive SMARTS (i.e. $(...))
       * encountered in the SMARTS.
       *
       * @return The list of recursive molecules.
       */
      const std::vector<HeMol>& recursiveMols() const
      {
        return m_recursiveMols;
      }

      /**
       * @brief Get the list of recursive SMARTS expression trees.
       *
       * These are the SMARTS expression trees for the recursive SMARTS
       * (i.e. $(...)) encountered in the SMARTS.
       *
       * @return The list of recursive SMARTS expression trees.
       */
      const std::vector<impl::SmartsTrees>& recursiveTrees() const
      {
        return m_recursiveTrees;
      }

      /**
       * @brief Get the error resulting from calling init().
       *
       * @return The parse error.
       */
      const Error& error() const
      {
        return m_error;
      }

      /**
       * @brief Get the atom class for an atom in the SMARTS query.
       *
       * @param index The index of the atom in the original SMARTS query.
       *
       * @return The atom's atom class or -1 if none was specified.
       */
      int atomClass(Index index) const;

      /**
       * @brief Check if the SMARTS has atom primitives that require a ring set.
       *
       * In OpenEye mode, the only atom primitive that requires a ring set is
       * the r<n> primitive (e.g. [r5], [r6], ...). In other modes the R<n>
       * atom primitive also requires a ring set.
       *
       * @return True if the SMARTS has cylce atom/bond primitives.
       */
      bool requiresRingSet() const
      {
        // computed and cached in init()
        return m_requiresRingSet;
      }

      /**
       * @brief Check if the SMARTS has cycle atom/bond primitives elements that
       *        require cyclicity information.
       *
       * Examples are [R3], [x2], *@*.
       *
       * @return True if the SMARTS has cylce atom/bond primitives.
       */
      bool requiresCyclicity() const
      {
        // computed and cached in init()
        return m_requiresCyclicity;
      }

      /**
       * @brief Check if the SMARTS has explicit hydrogens.
       *
       * @return True if the SMARTS has explicit hydrogens.
       */
      bool requiresExplicitHydrogens() const
      {
        // computed and cached in init()
        return m_requiresExplicitHydrogens;
      }

      /**
       * @brief Get the cycle match mode.
       */
      int mode() const
      {
        return m_mode;
      }

      /**
       * @brief Perform a SMARTS search on the specified molecule.
       *
       * This function requires a ring set as argument which is needed to match
       * the r<n> atom primitive. This ring set can be any ring set but the used
       * ring set determines the matching behavior of the R<n> atom primitive
       * (unless OpenEyeMode is used). When OpenEyeMode is used, OpenEye
       * semantics will be used for R<n> (i.e. R<n> is the same as x<n>).
       * Other cycle information will also be derived from the ring set.
       *
       * This function should be used when requiresRingSet() returns true.
       *
       * @param mol The queried molecule.
       * @param rings The ring set (needed for R<n>, r<n>, ...).
       * @param mapping The mapping to store the result.
       * @param uniqueComponents If true, an additional check will be performed
       *        to ensure the resulting MappingList(s) contain only unique maps.
       *        For example, SMARTS 'C.C' matched against 'C.C' will give 2
       *        mappings with uniqueComponents set to false and only 1 with
       *        uniqueComponents set to true. Similarly, SMARTS 'CO.CO.CC' will
       *        result in 4 maps for SMILES 'OCCCCCO' with uniqueComponents set
       *        to true and only 2 when set to false.
       *
       * @return True if the SMARTS matches the molecule.
       */
      template<typename MoleculeType, typename MappingType>
      bool findMapping(const MoleculeType &mol, const RingSet<MoleculeType> &rings,
          MappingType &mapping, bool uniqueComponents = false);

      /**
       * @brief Perform a SMARTS search on the specified molecule.
       *
       * When OpenEyeMode is set, OpenEye semantics will be used for R<n>
       * (i.e. R<n> is the same as x<n>).
       *
       * This function should be used when requiresRingSet() returns false but
       * requiresCyclicity() return true.
       *
       * @param mol The queried molecule.
       * @param cyclicAtoms The cyclic atoms.
       * @param cyclicBonds The cyclic bonds.
       * @param mapping The mapping to store the result.
       * @param uniqueComponents If true, an additional check will be performed
       *        to ensure the resulting MappingList(s) contain only unique maps.
       *        For example, SMARTS 'C.C' matched against 'C.C' will give 2
       *        mappings with uniqueComponents set to false and only 1 with
       *        uniqueComponents set to true. Similarly, SMARTS 'CO.CO.CC' will
       *        result in 4 maps for SMILES 'OCCCCCO' with uniqueComponents set
       *        to true and only 2 when set to false.
       *
       * @return True if the SMARTS matches the molecule.
       */
      template<typename MoleculeType, typename MappingType>
      bool findMapping(const MoleculeType &mol, const std::vector<bool> &cyclicAtoms,
          const std::vector<bool> &cyclicBonds, MappingType &mapping,
          bool uniqueComponents = false)
      {
        RingSet<MoleculeType> rings(mol, cyclicAtoms, cyclicBonds);
        return findMapping(mol, rings, mapping, uniqueComponents);
      }

      /**
       * @overload
       */
      template<typename MoleculeType, typename MappingType>
      bool findMapping(const MoleculeType &mol, MappingType &mapping,
          bool uniqueComponents = false)
      {
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

        switch (m_mode) {
          case SSSR:
            return findMapping(mol, relevant_cycles(mol), mapping, uniqueComponents); // FIXME
          case RelevantCycles:
            return findMapping(mol, relevant_cycles(mol), mapping, uniqueComponents);
          default:
          case OpenEye:
            {
              // get cycle membership
              std::vector<bool> cyclicAtoms(num_atoms(mol));
              std::vector<bool> cyclicBonds(num_bonds(mol));
              cycle_membership(mol, cyclicAtoms, cyclicBonds);

              return findMapping(mol, cyclicAtoms, cyclicBonds, mapping, uniqueComponents);
            }
        };
      }

      /**
       * @overload
       */
      template<typename MoleculeType>
      bool find(const MoleculeType &mol, const RingSet<MoleculeType> &rings,
          bool uniqueComponents = false)
      {
        NoMapping mapping;
        return findMapping(mol, rings, mapping, uniqueComponents);
      }

      /**
       * @overload
       */
      template<typename MoleculeType>
      bool find(const MoleculeType &mol, const std::vector<bool> &cyclicAtoms,
          const std::vector<bool> &cyclicBonds, bool uniqueComponents = false)
      {
        NoMapping mapping;
        return findMapping(mol, cyclicBonds, cyclicAtoms, mapping, uniqueComponents);
      }

      /**
       * @overload
       */
      template<typename MoleculeType>
      bool find(const MoleculeType &mol, bool uniqueComponents = false)
      {
        NoMapping mapping;
        return findMapping(mol, mapping, uniqueComponents);
      }

    private:
      HeMol m_query; //!< The query
      impl::SmartsTrees m_trees; //!< The SMARTS expression trees
      std::vector<HeMol> m_components; //!< The query components
      std::vector<impl::SmartsTrees> m_componentTrees; //!< The SMARTS component expression trees
      std::vector<HeMol> m_recursiveMols; //!< The recursive molecules
      std::vector<impl::SmartsTrees> m_recursiveTrees; //!< The recursive expression trees
      std::vector<std::vector<Index> > m_atomMaps; //!< m_components atom index to original smarts index
      //std::vector<std::vector<Index> > m_bondMaps; //!< m_components bond index to original smarts index
      Error m_error;
      int m_mode;
      bool m_requiresRingSet;
      bool m_requiresCyclicity;
      bool m_requiresExplicitHydrogens;
  };

  namespace impl {

    inline bool mappings_overlap(const IsomorphismMapping &map1, const IsomorphismMapping &map2)
    {
      for (std::size_t i = 0; i < map1.size(); ++i)
        if (std::find(map2.begin(), map2.end(), map1[i]) != map2.end())
          return true;
      return false;
    }

    template<typename MappingType>
    bool is_mapping_unique(const MappingType &mapping, const IsomorphismMapping &map)
    {
      return true;
    }

    template<>
    bool is_mapping_unique<MappingList>(const MappingList &mappings, const IsomorphismMapping &map);

    template<typename MappingType>
    void enumerate_mappings(std::size_t fragment, const std::vector<MappingList> &mappings,
        const IsomorphismMapping &current, const std::vector<std::vector<Index> > atomMaps,
        MappingType &output, bool &match, bool uniqueComponents)
    {
      for (std::size_t i = 0; i < mappings[fragment].maps.size(); ++i) {
        // check for overlap
        if (mappings_overlap(current, mappings[fragment].maps[i]))
          continue;
        // add mapping of considered fragment to current by translating indices
        IsomorphismMapping map = current;
        const IsomorphismMapping &next = mappings[fragment].maps[i];
        for (std::size_t j = 0; j < next.size(); ++j)
          map[atomMaps[fragment][j]] = next[j];

        if (fragment + 1 < mappings.size()) {
          // recursive call to consider next fragment
          enumerate_mappings(fragment + 1, mappings, map, atomMaps, output, match, uniqueComponents);
        } else {
          // last fragment done, add mapping to output
          assert(std::find(map.begin(), map.end(), -1) == map.end());

          // check if mapping is unique
          if (uniqueComponents && !is_mapping_unique(output, map))
            continue;
          impl::add_mapping(output, map);
          match = true;
        }
      }
    }

  }

  template<typename MoleculeType, typename MappingType>
  bool Smarts::findMapping(const MoleculeType &mol, const RingSet<MoleculeType> &rings,
      MappingType &mapping, bool uniqueComponents)
  {
    if (m_components.empty())
      return false;

    if (m_components.size() == 1) {
      // simple case: single SMARTS fragment
      impl::SmartsAtomMatcher<HeMol, MoleculeType> atomMatcher(m_componentTrees[0].atoms(),
          rings, m_recursiveMols, m_recursiveTrees, m_mode);
      impl::SmartsBondMatcher<HeMol, MoleculeType> bondMatcher(m_componentTrees[0].bonds(), rings);
      return isomorphism_search(mol, m_components[0], mapping, atomMatcher, bondMatcher);
    } else {
      // match each fragment seperatly and store results in mappings
      int numQueryAtoms = 0;
      std::vector<MappingList> mappings(m_components.size());
      for (std::size_t i = 0; i < m_components.size(); ++i) {
        numQueryAtoms += num_atoms(m_components[i]);
        impl::SmartsAtomMatcher<HeMol, MoleculeType> atomMatcher(m_componentTrees[i].atoms(),
            rings, m_recursiveMols, m_recursiveTrees, m_mode);
        impl::SmartsBondMatcher<HeMol, MoleculeType> bondMatcher(m_componentTrees[i].bonds(), rings);
        if (!isomorphism_search(mol, m_components[i], mappings[i], atomMatcher, bondMatcher))
          return false;
      }

      bool match = false;
      IsomorphismMapping map(numQueryAtoms, -1);
      impl::enumerate_mappings(0, mappings, map, m_atomMaps, mapping, match, uniqueComponents);

      return match;
    }
  }

}

#endif
