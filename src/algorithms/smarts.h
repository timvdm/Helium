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

#include <iostream>

#include <Helium/hemol.h>
#include <Helium/ring.h>
#include <Helium/algorithms/isomorphism.h>

#include <smiley.h>

namespace Helium {

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
        void *recursive;
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
        ~SmartsTrees()
        {
          for (std::size_t i = 0; i < m_atoms.size(); ++i)
            cleanup(m_atoms[i]);
          for (std::size_t i = 0; i < m_bonds.size(); ++i)
            cleanup(m_bonds[i]);
        }

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
            //case AE_Recursive:
            //  cleanup(expr->recursive);
            //  break;
            default:
              break;
          }
          delete expr;
        }

        std::vector<SmartsAtomExpr*> m_atoms;
        std::vector<SmartsBondExpr*> m_bonds;
    };

    struct SmartsCallback : public Smiley::CallbackBase
    {
      SmartsCallback(HeMol &mol_, SmartsTrees &trees_)
        : mol(mol_), trees(trees_)
      {
      }

      //@name SMILES/SMARTS
      //@{
      /**
       * Prepare the callback functor for a new SMILES/SMARTS. This method is
       * always invoked at the start of parsing before any of the other methods.
       */
      void clear()
      {
        bondRing = false;
      }
      /**
       * Set the chirality for an atom. This method is invoked when the entire
       * SMILES/SMARTS is parsed.
       */
      void setChiral(int index, Smiley::Chirality chirality, const std::vector<int> &chiralNbrs) {}
      void end() {}
      //@}
      //@name SMILES
      /**
       * Invoked when an atom is completly parsed.
       */
      void addAtom(int element, bool aromatic, int isotope, int hCount, int charge, int atomClass) {}

      void processBond()
      {
        // pop stack
        while (!bondStack.empty()) {
          bondPostfix.push_back(bondStack.back());
          bondStack.pop_back();
        }

        // convert postfix to tree
        typedef std::vector<SmartsBondExpr*>::const_iterator PostfixIter;
        for (PostfixIter expr = bondPostfix.begin(); expr != bondPostfix.end(); ++expr) {
          if (isUnaryOp((*expr)->type)) {
            assert(bondStack.size() >= 1);

            (*expr)->arg = bondStack.back();
            bondStack.pop_back();

            bondStack.push_back(*expr);
          } else if (isBinaryOp((*expr)->type)) {
            assert(bondStack.size() >= 2);

            (*expr)->left = bondStack.back();
            bondStack.pop_back();
            (*expr)->right = bondStack.back();
            bondStack.pop_back();
            bondStack.push_back(*expr);
          } else
            bondStack.push_back(*expr);
        }

        bondPostfix.clear();
      }

      /**
       * Invoked once both bon atom are added using addAtom().
       */
      void addBond(int source, int target, int order, bool isUp, bool isDown)
      {
        mol.addBond(mol.atom(source), mol.atom(target));

        // process parsed bond primitives
        if (!bondRing)
          processBond();
        bondRing = false;

        // default bond
        if (bondStack.empty()) {
          if (order == 5) {
            // note: cc means any pair of attached aromatic carbons (i.e. c-c or c:c)
            SmartsBondExpr *single = new SmartsBondExpr(Smiley::BE_Single);
            SmartsBondExpr *aromatic = new SmartsBondExpr(Smiley::BE_Aromatic);
            bondStack.push_back(new SmartsBondExpr(Smiley::OP_Or, single, aromatic));
          } else
            bondStack.push_back(new SmartsBondExpr(Smiley::BE_Single));
        }

        assert(!bondStack.empty());
        trees.addBond(bondStack.back());
        bondStack.clear();
        bondPostfix.clear();
      }
      //@}

      int precedence(int type)
      {
        switch (type) {
          case Smiley::OP_Not:
            return 0;
          case Smiley::OP_AndHi:
            return 1;
          case Smiley::OP_Or:
            return 2;
          case Smiley::OP_AndLo:
            return 3;
        }
        return 0;
      }

      //@name SMARTS
      //@{
      /**
       * Invoked when a unary or binary logical operator is parsed
       * (i.e. '&', ';' or ','). This method is also invoked for implicit AND.
       */
      void atomOperation(int type)
      {
        // pop until top of stack has lower precedence
        while (atomStack.size() && precedence(atomStack.back()->type) < precedence(type)) {
          atomPostfix.push_back(atomStack.back());
          atomStack.pop_back();
        }

        // push operator to stack
        atomStack.push_back(new SmartsAtomExpr(type));
      }

      void bondOperation(int type)
      {
        // pop until top of stack has lower precedence
        while (bondStack.size() && precedence(bondStack.back()->type) < precedence(type)) {
          bondPostfix.push_back(bondStack.back());
          bondStack.pop_back();
        }

        // push operator to stack
        bondStack.push_back(new SmartsBondExpr(type));
      }

      /**
       * Invoked when an unbracketed atom (i.e. organic subset) is parsed.
       */
      void addOrganicSubsetAtom(int element, bool aromatic)
      {
        // add an atom to the molecule
        mol.addAtom();

        // create the associated smarts expression tree
        if (element == -1) {
          int type = aromatic ? Smiley::AE_Aromatic : Smiley::AE_Aliphatic;
          SmartsAtomExpr *expr = new SmartsAtomExpr(type);
          trees.addAtom(expr);
        } else if (element == 0) {
          SmartsAtomExpr *expr = new SmartsAtomExpr(Smiley::AE_True);
          trees.addAtom(expr);
        } else {
          int type = aromatic ? Smiley::AE_AromaticElement : Smiley::AE_AliphaticElement;
          SmartsAtomExpr *expr = new SmartsAtomExpr(type, element);
          trees.addAtom(expr);
        }
      }

      /**
       * Invoked when an atom primitive is parsed.
       */
      void beginAtom()
      {
        // add an atom to the molecule
        mol.addAtom();

        atomStack.clear();
        atomPostfix.clear();
      }

      void atomPrimitive(int type, int value)
      {
        // add operand to postfix
        atomPostfix.push_back(new SmartsAtomExpr(type, value));
      }

      bool isBinaryOp(int type)
      {
        switch (type) {
          case Smiley::OP_AndHi:
          case Smiley::OP_AndLo:
          case Smiley::OP_Or:
            return true;
          default:
            return false;
        }
      }

      bool isUnaryOp(int type)
      {
        switch (type) {
          case Smiley::OP_Not:
            return true;
          default:
            return false;
        }
      }

      void endAtom()
      {
        // pop stack
        while (!atomStack.empty()) {
          atomPostfix.push_back(atomStack.back());
          atomStack.pop_back();
        }

        // convert postfix to tree
        typedef std::vector<SmartsAtomExpr*>::const_iterator PostfixIter;
        for (PostfixIter expr = atomPostfix.begin(); expr != atomPostfix.end(); ++expr) {
          if (isUnaryOp((*expr)->type)) {
            assert(atomStack.size() >= 1);

            (*expr)->arg = atomStack.back();
            atomStack.pop_back();

            atomStack.push_back(*expr);
          } else if (isBinaryOp((*expr)->type)) {
            assert(atomStack.size() >= 2);

            (*expr)->left = atomStack.back();
            atomStack.pop_back();
            (*expr)->right = atomStack.back();
            atomStack.pop_back();
            atomStack.push_back(*expr);
          } else
            atomStack.push_back(*expr);
        }

        assert(atomStack.size() == 1);

        trees.addAtom(atomStack.back());
        atomStack.clear();
      }

      /**
       * Invoked when a bond primitive is parsed. This method is also invoked for
       * implicit bonds.
       */
      void bondPrimitive(int type)
      {
        // add operand to postfix
        bondPostfix.push_back(new SmartsBondExpr(type));
      }

      /**
       * Invoked when a new ring bond number is parsed.
       */
      void startRingBond(int number)
      {
        processBond();
        if (!bondStack.empty()) {
          bondRings[number] = bondStack.back();
          bondStack.clear();
        }

        assert(bondStack.empty());
        assert(bondPostfix.empty());
      }
      /**
       * Invoked when a prviously found ring bond number is parsed to add the bond.
       */
      void endRingBond(int number)
      {
        std::map<int, SmartsBondExpr*>::iterator expr = bondRings.find(number);
        if (expr != bondRings.end()) {
          bondStack.push_back(expr->second);
          assert(bondStack.size() == 1);
          bondRings.erase(expr);
          bondRing = true;
        }

      }
      //@}


      HeMol &mol;
      SmartsTrees &trees;

      std::vector<SmartsAtomExpr*> atomStack;
      std::vector<SmartsAtomExpr*> atomPostfix;
      std::vector<SmartsBondExpr*> bondStack;
      std::vector<SmartsBondExpr*> bondPostfix;
      std::map<int, SmartsBondExpr*> bondRings;
      bool bondRing;
    };

    template<typename QueryType, typename MoleculeType>
    class SmartsAtomMatcher
    {
      public:
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<QueryType>::atom_type query_atom_type;

        SmartsAtomMatcher(const std::vector<SmartsAtomExpr*> &atoms, const RingSet<MoleculeType> &rings)
          : m_atoms(atoms), m_rings(rings)
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
                FOREACH_NBR (nbr, atom, mol, MoleculeType)
                  if (get_element(mol, *nbr) == 1)
                    ++h;
                return (h + num_hydrogens(mol, atom)) == expr->value;
              }
            case Smiley::AE_ImplicitH:
              if (expr->value == -1) // default: at least 1
                return num_hydrogens(mol, atom) >= 1;
              else
                return num_hydrogens(mol, atom) == expr->value;
            case Smiley::AE_RingMembership:
              return m_rings.numRings(atom) == expr->value;
            case Smiley::AE_RingSize:
              return m_rings.isAtomInRingSize(atom, expr->value);
            case Smiley::AE_RingConnectivity:
              if (expr->value == -1) // default: at least 1
                return m_rings.numRingBonds(atom) >= 1;
              else
                return m_rings.numRingBonds(atom) == expr->value;
            case Smiley::AE_Charge:
              return get_charge(mol, atom) == expr->value;
            case Smiley::AE_Chirality:
              return true;
            case Smiley::AE_AtomClass:
              return true;
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

  class Smarts
  {
    public:
      Smarts(const std::string &smarts)
      {
        impl::SmartsCallback callback(m_query, m_trees);
        Smiley::Parser<impl::SmartsCallback> parser(callback, Smiley::Parser<impl::SmartsCallback>::SmartsMode);

        try {
          parser.parse(smarts);
        } catch (Smiley::Exception &e) {
          std::ostringstream errorStream;
          if (e.type() == Smiley::Exception::SyntaxError)
            errorStream << "Syntax";
          else
            errorStream << "Semantics";
          errorStream << "Error: " << e.what() << "." << std::endl;
          errorStream << smarts << std::endl;
          for (std::size_t i = 0; i < e.pos(); ++i)
            errorStream << " ";
          for (std::size_t i = 0; i < e.length(); ++i)
            errorStream << "^";
          errorStream << std::endl;

          // FIXME... bool init(...) -> getError()
          // Throw new exception with added details ...
          throw Smiley::Exception(e.type(), e.errorCode(),
              errorStream.str(), e.pos(), e.length());
        }
      }

      HeMol& query()
      {
        return m_query;
      }

      impl::SmartsTrees& trees()
      {
        return m_trees;
      }

      template<typename MoleculeType, typename MappingType>
      bool search(MoleculeType &mol, MappingType &mapping, const RingSet<MoleculeType> &rings);

      template<typename MoleculeType>
      bool search(MoleculeType &mol, const RingSet<MoleculeType> &rings)
      {
        NoMapping mapping;
        return search(mol, mapping, rings);
      }

    private:
      HeMol m_query;
      impl::SmartsTrees m_trees;
  };


  template<typename MoleculeType, typename MappingType>
  bool Smarts::search(MoleculeType &mol, MappingType &mapping, const RingSet<MoleculeType> &rings)
  {
    impl::SmartsAtomMatcher<HeMol, MoleculeType> atomMatcher(m_trees.atoms(), rings);
    impl::SmartsBondMatcher<HeMol, MoleculeType> bondMatcher(m_trees.bonds(), rings);
    return isomorphism_search(mol, m_query, mapping, atomMatcher, bondMatcher);
  }

}

#endif
