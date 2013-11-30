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
#include <Helium/hemol.h>
#include <Helium/smiles.h>

namespace Helium {
  
  //@cond dev

  molecule_traits<HeMol>::atom_type HeMol::addAtom()
  {
    Index index = m_element.size();
    m_adjList.resize(m_adjList.size() + 1);
    m_atomAromatic.resize(m_atomAromatic.size() + 1);
    m_element.resize(m_element.size() + 1);
    m_mass.resize(m_mass.size() + 1);
    m_hydrogens.resize(m_hydrogens.size() + 1);
    m_charge.resize(m_charge.size() + 1);

    return atom_type(this, index);
  }

  molecule_traits<HeMol>::bond_type HeMol::addBond(const molecule_traits<HeMol>::atom_type &source,
                                                   const molecule_traits<HeMol>::atom_type &target)
  {
    Index index = m_order.size();

    m_source.push_back(source.index());
    m_target.push_back(target.index());
    m_bondAromatic.resize(m_bondAromatic.size() + 1);
    m_order.resize(m_order.size() + 1);

    bond_type bond(this, index);

    m_adjList[source.index()].push_back(bond);
    m_adjList[target.index()].push_back(bond);

    return bond;
  }

  void HeMol::clear()
  {
    m_adjList.clear();
    m_atomAromatic.clear();
    m_element.clear();
    m_mass.clear();
    m_hydrogens.clear();
    m_charge.clear();

    m_source.clear();
    m_target.clear();
    m_bondAromatic.clear();
    m_order.clear();
  }

  void HeMol::renumberAtoms(const std::vector<Index> &permutation)
  {
    assert(permutation.size() == m_adjList.size());
    impl::apply_permutation(m_adjList, permutation);
    impl::apply_permutation(m_atomAromatic, permutation);
    impl::apply_permutation(m_element, permutation);
    impl::apply_permutation(m_mass, permutation);
    impl::apply_permutation(m_hydrogens, permutation);
    impl::apply_permutation(m_charge, permutation);

    std::vector<Index> copy;
    for (std::size_t i = 0; i < m_source.size(); ++i)
      copy.push_back(index_of(permutation, m_source[i]));
    m_source = copy;

    copy.clear();
    for (std::size_t i = 0; i < m_target.size(); ++i)
      copy.push_back(index_of(permutation, m_target[i]));
    m_target = copy;
  }

  void hemol_from_smiles(const std::string &smiles, HeMol &mol)
  {
    try {
      parse_smiles(smiles, mol);
    } catch(Smiley::Exception &e) {
      std::cerr << e.what();
    }
  }

  HeMol hemol_from_smiles(const std::string &smiles)
  {
    HeMol mol;
    hemol_from_smiles(smiles, mol);
    return mol;
  }
  
  //@endcond dev

}
