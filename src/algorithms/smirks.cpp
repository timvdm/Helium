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
#include <Helium/algorithms/smirks.h>

#include <iostream>

namespace Helium {

  bool Smirks::init(const std::string &smirks)
  {
    std::size_t pos = smirks.find(">>");
    if (pos == std::string::npos) {
      m_error = SmirksError(SmirksError::NoReaction, "SMIRKS does not contain '>>'");
      return false;
    }

    return init(smirks.substr(0, pos), smirks.substr(pos + 2));
  }

  bool Smirks::init(const std::string &reactant, const std::string &product)
  {
    //std::cout << "Smirks::init(" << reactant << ">>" << product << ")" << std::endl;

    m_error = SmirksError();

    if (!m_reactant.init(reactant)) {
      m_error = SmirksError(SmirksError::ReactantSmarts, m_reactant.error().what());
      return false;
    }
    if (!m_product.init(product)) {
      m_error = SmirksError(SmirksError::ProductSmarts, m_product.error().what());
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
        m_error = SmirksError(SmirksError::ProductContainsNot, "Product may not contain NOT expression");
        return false;
      }

    //
    // Detect bond changes
    //

    // iterate over reactant bonds: find removed and changed bonds
    for (std::size_t i = 0; i < num_bonds(m_reactant.query()); ++i) {
      molecule_traits<HeMol>::bond_type bond = get_bond(m_reactant.query(), i);
      molecule_traits<HeMol>::atom_type source = get_source(m_reactant.query(), bond);
      molecule_traits<HeMol>::atom_type target = get_target(m_reactant.query(), bond);
      int sourceAtomClass = m_reactant.atomClass(get_index(m_reactant.query(), source));
      int targetAtomClass = m_reactant.atomClass(get_index(m_reactant.query(), target));

      // find matching bond in product
      molecule_traits<HeMol>::bond_type productBond = getBond(m_product, sourceAtomClass, targetAtomClass);

      if (productBond == molecule_traits<HeMol>::null_bond()) {
        // bond is removed
        m_bondChanges.push_back(BondChange(BondChange::Removed,
              get_index(m_reactant.query(), source),
              get_index(m_reactant.query(), target), 0));
      } else {
        // bond is (possibly) changed
        m_bondChanges.push_back(BondChange(BondChange::Changed,
              get_index(m_reactant.query(), source),
              get_index(m_reactant.query(), target),
              m_product.trees().bond(get_index(m_product.query(), productBond))));
      }
    }

    // iterate over reactant bonds: find added bonds
    for (std::size_t i = 0; i < num_bonds(m_product.query()); ++i) {
      // check if bond expr is "simple"
      switch (m_product.trees().bond(i)->type) {
        case Smiley::BE_Single:
        case Smiley::BE_Double:
        case Smiley::BE_Triple:
        case Smiley::BE_Quadriple:
        case Smiley::BE_Aromatic:
        case Smiley::BE_Up:
        case Smiley::BE_Down:
          break;
        default:
          m_error = SmirksError(SmirksError::InvalidProductBond, "Bonds in product must be valid SMILES");
          return false;
      }

      molecule_traits<HeMol>::bond_type bond = get_bond(m_product.query(), i);
      molecule_traits<HeMol>::atom_type source = get_source(m_product.query(), bond);
      molecule_traits<HeMol>::atom_type target = get_target(m_product.query(), bond);
      int sourceAtomClass = m_product.atomClass(get_index(m_product.query(), source));
      int targetAtomClass = m_product.atomClass(get_index(m_product.query(), target));

      // find matching bond in reactant
      molecule_traits<HeMol>::bond_type reactantBond = getBond(m_reactant, sourceAtomClass, targetAtomClass);

      if (reactantBond == molecule_traits<HeMol>::null_bond()) {
        int sourceIndex = molecule_traits<HeMol>::null_index();
        for (std::size_t j = 0; j < num_atoms(m_reactant.query()); ++j)
          if (m_reactant.atomClass(j) == sourceAtomClass) {
            sourceIndex = j;
            break;
          }

        int targetIndex = molecule_traits<HeMol>::null_index();
        for (std::size_t j = 0; j < num_atoms(m_reactant.query()); ++j)
          if (m_reactant.atomClass(j) == targetAtomClass) {
            targetIndex = j;
            break;
          }

        assert(sourceIndex != molecule_traits<HeMol>::null_index());
        assert(targetIndex != molecule_traits<HeMol>::null_index());

        // bond is added
        m_bondChanges.push_back(BondChange(BondChange::Added, sourceIndex, targetIndex,
              m_product.trees().bond(i)));
      }
    }

    return true;
  }


  std::map<int, impl::SmartsAtomExpr*> Smirks::atomClassToExpr(const Smarts &smarts) const
  {
    std::map<int, impl::SmartsAtomExpr*> result;

    for (std::size_t i = 0; i < num_atoms(smarts.query()); ++i) {
      int atomClass = smarts.atomClass(i);
      if (atomClass != -1)
        result[atomClass] = smarts.trees().atom(i);
    }

    return result;
  }

  molecule_traits<HeMol>::bond_type Smirks::getBond(const Smarts &smarts, int sourceAtomClass, int targetAtomClass)
  {
    const HeMol &mol = smarts.query();
    for (std::size_t i = 0; i < num_bonds(mol); ++i) {
      molecule_traits<HeMol>::bond_type bond = get_bond(mol, i);
      molecule_traits<HeMol>::atom_type source = get_source(mol, bond);
      molecule_traits<HeMol>::atom_type target = get_target(mol, bond);

      if (smarts.atomClass(get_index(mol, source)) == sourceAtomClass &&
          smarts.atomClass(get_index(mol, target)) == targetAtomClass)
        return bond;

      if (smarts.atomClass(get_index(mol, source)) == targetAtomClass &&
          smarts.atomClass(get_index(mol, target)) == sourceAtomClass)
        return bond;
    }

    return molecule_traits<HeMol>::null_bond();
  }

}
