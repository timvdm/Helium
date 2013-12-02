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
#include <Helium/algorithms/smarts.h>

namespace Helium {

  namespace impl {

    SmartsTrees::SmartsTrees()
    {
    }

    SmartsTrees::SmartsTrees(const SmartsTrees &other)
    {
      for (std::size_t i = 0; i < other.m_atoms.size(); ++i)
        m_atoms.push_back(copy(other.m_atoms[i]));
      for (std::size_t i = 0; i < other.m_bonds.size(); ++i)
        m_bonds.push_back(copy(other.m_bonds[i]));
    }

    SmartsTrees::~SmartsTrees()
    {
      clear();
    }

    SmartsTrees& SmartsTrees::operator=(const SmartsTrees &other)
    {
      for (std::size_t i = 0; i < other.m_atoms.size(); ++i)
        m_atoms.push_back(copy(other.m_atoms[i]));
      for (std::size_t i = 0; i < other.m_bonds.size(); ++i)
        m_bonds.push_back(copy(other.m_bonds[i]));
      return *this;
    }

    void SmartsTrees::clear()
    {
      for (std::size_t i = 0; i < m_atoms.size(); ++i)
        cleanup(m_atoms[i]);
      for (std::size_t i = 0; i < m_bonds.size(); ++i)
        cleanup(m_bonds[i]);
      m_atoms.clear();
      m_bonds.clear();
    }

    SmartsAtomExpr* SmartsTrees::copy(SmartsAtomExpr *expr)
    {
      switch (expr->type) {
        case Smiley::OP_AndHi:
        case Smiley::OP_AndLo:
        case Smiley::OP_And:
        case Smiley::OP_Or:
          return new SmartsAtomExpr(expr->type, copy(expr->left), copy(expr->right));
        case Smiley::OP_Not:
          return new SmartsAtomExpr(expr->type, copy(expr->arg));
        default:
          return new SmartsAtomExpr(expr->type, expr->value);
      }
    }

    SmartsBondExpr* SmartsTrees::copy(SmartsBondExpr *expr)
    {
      switch (expr->type) {
        case Smiley::OP_AndHi:
        case Smiley::OP_AndLo:
        case Smiley::OP_And:
        case Smiley::OP_Or:
          return new SmartsBondExpr(expr->type, copy(expr->left), copy(expr->right));
        case Smiley::OP_Not:
          return new SmartsBondExpr(expr->type, copy(expr->arg));
        default:
          return new SmartsBondExpr(expr->type);
      }
    }

    struct SmartsCallback : public Smiley::CallbackBase
    {
      public:
        struct State
        {
          State() : bondRing(false)
          {
          }

          std::vector<SmartsAtomExpr*> atomStack;
          std::vector<SmartsAtomExpr*> atomPostfix;
          std::vector<SmartsBondExpr*> bondStack;
          std::vector<SmartsBondExpr*> bondPostfix;
          std::map<int, SmartsBondExpr*> bondRings;
          bool bondRing;
        };

        SmartsCallback(HeMol &mol, SmartsTrees &trees, std::vector<HeMol> &recursiveMols,
            std::vector<SmartsTrees> &recursiveTrees)
          : m_mol(mol), m_trees(trees), m_recursiveMols(recursiveMols), m_recursiveTrees(recursiveTrees)
        {
        }

        State& state()
        {
          return m_states.back();
        }

        HeMol& mol()
        {
          if (m_states.size() == 1)
            return m_mol;
          assert(!m_recursiveMols.empty());
          return m_recursiveMols.back();
        }

        SmartsTrees& trees()
        {
          if (m_states.size() == 1)
            return m_trees;
          assert(!m_recursiveTrees.empty());
          return m_recursiveTrees.back();
        }

        std::vector<SmartsAtomExpr*>& atomStack()
        {
          return state().atomStack;
        }

        std::vector<SmartsAtomExpr*>& atomPostfix()
        {
          return state().atomPostfix;
        }

        std::vector<SmartsBondExpr*>& bondStack()
        {
          return state().bondStack;
        }

        std::vector<SmartsBondExpr*>& bondPostfix()
        {
          return state().bondPostfix;
        }

        std::map<int, SmartsBondExpr*>& bondRings()
        {
          return state().bondRings;
        }

        bool& bondRing()
        {
          return state().bondRing;
        }


        // SMILES/SMARTS
        /**
         * Prepare the callback functor for a new SMILES/SMARTS. This method is
         * always invoked at the start of parsing before any of the other methods.
         */
        void clear()
        {
          m_states.clear();
          m_states.resize(1);
          bondRing() = false;
        }
        /**
         * Set the chirality for an atom. This method is invoked when the entire
         * SMILES/SMARTS is parsed.
         */
        void setChiral(int index, Smiley::Chirality chirality, const std::vector<int> &chiralNbrs) {}
        void end() {}
        // SMILES
        /**
         * Invoked when an atom is completly parsed.
         */
        void addAtom(int element, bool aromatic, int isotope, int hCount, int charge, int atomClass) {}

        void processBond()
        {
          // pop stack
          while (!bondStack().empty()) {
            bondPostfix().push_back(bondStack().back());
            bondStack().pop_back();
          }

          // convert postfix to tree
          typedef std::vector<SmartsBondExpr*>::const_iterator PostfixIter;
          for (PostfixIter expr = bondPostfix().begin(); expr != bondPostfix().end(); ++expr) {
            if (isUnaryOp((*expr)->type)) {
              assert(bondStack().size() >= 1);

              (*expr)->arg = bondStack().back();
              bondStack().pop_back();

              bondStack().push_back(*expr);
            } else if (isBinaryOp((*expr)->type)) {
              assert(bondStack().size() >= 2);

              (*expr)->left = bondStack().back();
              bondStack().pop_back();
              (*expr)->right = bondStack().back();
              bondStack().pop_back();
              bondStack().push_back(*expr);
            } else
              bondStack().push_back(*expr);
          }

          bondPostfix().clear();
        }

        /**
         * Invoked once both bon atom are added using addAtom().
         */
        void addBond(int source, int target, int order, bool isUp, bool isDown)
        {
          mol().addBond(mol().atom(source), mol().atom(target));

          // process parsed bond primitives
          if (!bondRing())
            processBond();
          bondRing() = false;

          // default bond
          if (bondStack().empty()) {
            if (order == 5) {
              // note: cc means any pair of attached aromatic carbons (i.e. c-c or c:c)
              SmartsBondExpr *single = new SmartsBondExpr(Smiley::BE_Single);
              SmartsBondExpr *aromatic = new SmartsBondExpr(Smiley::BE_Aromatic);
              bondStack().push_back(new SmartsBondExpr(Smiley::OP_Or, single, aromatic));
            } else
              bondStack().push_back(new SmartsBondExpr(Smiley::BE_Single));
          }

          assert(!bondStack().empty());
          trees().addBond(bondStack().back());
          bondStack().clear();
          bondPostfix().clear();
        }

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

        // SMARTS
        /**
         * Invoked when a unary or binary logical operator is parsed
         * (i.e. '&', ';' or ','). This method is also invoked for implicit AND.
         */
        void atomOperation(int type)
        {
          // pop until top of stack has lower precedence
          while (atomStack().size() && precedence(atomStack().back()->type) < precedence(type)) {
            atomPostfix().push_back(atomStack().back());
            atomStack().pop_back();
          }

          // push operator to stack
          atomStack().push_back(new SmartsAtomExpr(type));
        }

        void bondOperation(int type)
        {
          // pop until top of stack has lower precedence
          while (bondStack().size() && precedence(bondStack().back()->type) < precedence(type)) {
            bondPostfix().push_back(bondStack().back());
            bondStack().pop_back();
          }

          // push operator to stack
          bondStack().push_back(new SmartsBondExpr(type));
        }

        /**
         * Invoked when an unbracketed atom (i.e. organic subset) is parsed.
         */
        void addOrganicSubsetAtom(int element, bool aromatic)
        {
          // add an atom to the molecule
          mol().addAtom();

          // create the associated smarts expression tree
          if (element == -1) {
            int type = aromatic ? Smiley::AE_Aromatic : Smiley::AE_Aliphatic;
            SmartsAtomExpr *expr = new SmartsAtomExpr(type);
            trees().addAtom(expr);
          } else if (element == 0) {
            SmartsAtomExpr *expr = new SmartsAtomExpr(Smiley::AE_True);
            trees().addAtom(expr);
          } else {
            int type = aromatic ? Smiley::AE_AromaticElement : Smiley::AE_AliphaticElement;
            SmartsAtomExpr *expr = new SmartsAtomExpr(type, element);
            trees().addAtom(expr);
          }
        }

        /**
         * Invoked when an atom primitive is parsed.
         */
        void beginAtom()
        {
          // add an atom to the molecule
          mol().addAtom();

          atomStack().clear();
          atomPostfix().clear();
        }

        void atomPrimitive(int type, int value)
        {
          // add operand to postfix
          atomPostfix().push_back(new SmartsAtomExpr(type, value));
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
          while (!atomStack().empty()) {
            atomPostfix().push_back(atomStack().back());
            atomStack().pop_back();
          }

          // convert postfix to tree
          typedef std::vector<SmartsAtomExpr*>::const_iterator PostfixIter;
          for (PostfixIter expr = atomPostfix().begin(); expr != atomPostfix().end(); ++expr) {
            if (isUnaryOp((*expr)->type)) {
              assert(atomStack().size() >= 1);

              (*expr)->arg = atomStack().back();
              atomStack().pop_back();

              atomStack().push_back(*expr);
            } else if (isBinaryOp((*expr)->type)) {
              assert(atomStack().size() >= 2);

              (*expr)->left = atomStack().back();
              atomStack().pop_back();
              (*expr)->right = atomStack().back();
              atomStack().pop_back();
              atomStack().push_back(*expr);
            } else
              atomStack().push_back(*expr);
          }

          assert(atomStack().size() == 1);

          trees().addAtom(atomStack().back());
          atomStack().clear();
        }

        /**
         * Invoked when a bond primitive is parsed. This method is also invoked for
         * implicit bonds.
         */
        void bondPrimitive(int type)
        {
          // add operand to postfix
          bondPostfix().push_back(new SmartsBondExpr(type));
        }

        /**
         * Invoked when a new ring bond number is parsed.
         */
        void startRingBond(int number)
        {
          processBond();
          if (!bondStack().empty()) {
            bondRings()[number] = bondStack().back();
            bondStack().clear();
          }

          assert(bondStack().empty());
          assert(bondPostfix().empty());
        }
        /**
         * Invoked when a prviously found ring bond number is parsed to add the bond.
         */
        void endRingBond(int number)
        {
          std::map<int, SmartsBondExpr*>::iterator expr = bondRings().find(number);
          if (expr != bondRings().end()) {
            bondStack().push_back(expr->second);
            assert(bondStack().size() == 1);
            bondRings().erase(expr);
            bondRing() = true;
          }
        }

        void pushState()
        {
          // add recusrive atom expr to postfix
          int index = m_recursiveMols.size();
          atomPostfix().push_back(new SmartsAtomExpr(Smiley::AE_Recursive, index));
          // push state
          m_states.push_back(State());
          m_recursiveMols.resize(m_recursiveMols.size() + 1);
          m_recursiveTrees.resize(m_recursiveTrees.size() + 1);
        }

        void popState()
        {
          m_states.pop_back();
        }

    private:
        HeMol &m_mol;
        SmartsTrees &m_trees;
        std::vector<HeMol> &m_recursiveMols;
        std::vector<SmartsTrees> &m_recursiveTrees;
        std::vector<State> m_states;
    };

  }

  bool Smarts::init(const std::string &smarts)
  {
    // clear error, ...
    m_error = Smiley::Exception();

    m_query.clear();
    m_trees.clear();
    m_components.clear();
    m_componentTrees.clear();
    m_recursiveMols.clear();
    m_recursiveTrees.clear();
    m_atomMaps.clear();
    //m_bondMaps.clear();

    // create SMARTS parser callback
    impl::SmartsCallback callback(m_query, m_trees, m_recursiveMols, m_recursiveTrees);
    // create the parser
    Smiley::Parser<impl::SmartsCallback> parser(callback, Smiley::Parser<impl::SmartsCallback>::SmartsMode);

    try {
      // try to parse
      parser.parse(smarts);
    } catch (Smiley::Exception &e) {
      // parsing failed: do some error handling
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

      m_error = Smiley::Exception(e.type(), e.errorCode(),
          errorStream.str(), e.pos(), e.length());

      return false;
    }

    int numComponents = num_connected_components(m_query);

    if (numComponents <= 1) {
      m_components.push_back(m_query);
      m_componentTrees.push_back(m_trees);
    } else {
      m_atomMaps.resize(numComponents);
      //m_bondMaps.resize(numComponents);

      std::vector<unsigned int> atomComponents = connected_atom_components(m_query);
      std::vector<unsigned int> bondComponents = connected_bond_components(m_query);

      for (unsigned int c = 0; c < numComponents; ++c) {
        m_atomMaps[c].resize(std::count(atomComponents.begin(), atomComponents.end(), c));
        //m_bondMaps[c].resize(std::count(bondComponents.begin(), bondComponents.end(), c));

        HeMol q;
        impl::SmartsTrees t;

        // original smarts index to m_components index
        std::map<Index, Index> reverseAtomMap;

        Index index = 0;
        for (Index i = 0; i < num_atoms(m_query); ++i) {
          if (atomComponents[i] == c) {
            reverseAtomMap[i] = index;
            m_atomMaps[c][index++] = i;
            q.addAtom();
            t.addAtom(m_trees.copy(m_trees.atom(i)));
          }
        }

        //index = 0;
        for (Index i = 0; i < num_bonds(m_query); ++i) {
          if (bondComponents[i] == c) {
            //m_bondMaps[c][index++] = i;
            molecule_traits<HeMol>::bond_type bond = get_bond(m_query, i);
            molecule_traits<HeMol>::atom_type source = get_atom(q, reverseAtomMap[get_index(m_query, get_source(m_query, bond))]);
            molecule_traits<HeMol>::atom_type target = get_atom(q, reverseAtomMap[get_index(m_query, get_target(m_query, bond))]);
            q.addBond(source, target);
            t.addBond(m_trees.copy(m_trees.bond(i)));
          }
        }

        assert(num_atoms(q) == t.atoms().size());
        assert(num_bonds(q) == t.bonds().size());

        m_components.push_back(q);
        m_componentTrees.push_back(t);
      }
    }

    return true;
  }

  namespace impl {

    int extract_atom_class(SmartsAtomExpr *expr)
    {
      int cls = -1;
      switch (expr->type) {
        case Smiley::AE_AtomClass:
          return expr->value;
        case Smiley::OP_AndHi:
        case Smiley::OP_AndLo:
        case Smiley::OP_Or:
          cls = extract_atom_class(expr->left);
          if (cls != -1)
            return cls;
          cls = extract_atom_class(expr->right);
          if (cls != -1)
            return cls;
          return -1;
        default:
          return -1;
      }
    }

  }

  int Smarts::atomClass(Index index) const
  {
    assert(index < num_atoms(m_query));
    return impl::extract_atom_class(m_trees.atom(index));
  }

  namespace impl {

    template<>
    bool is_mapping_unique<MappingList>(const MappingList &mappings, const IsomorphismMapping &map)
    {
      for (std::size_t i = 0; i < mappings.maps.size(); ++i) {
        const IsomorphismMapping &ref = mappings.maps[i];

        bool unique = false;
        for (std::size_t j = 0; j < map.size(); ++j)
          if (std::find(ref.begin(), ref.end(), map[j]) == ref.end()) {
            unique = true;
            break;
          }

        if (!unique)
          return false;
      }

      return true;
    }

  }

}
