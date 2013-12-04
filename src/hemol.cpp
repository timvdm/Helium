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
  HeMol::HeMol(const HeMol &other)
  {
    // copy m_adjList
    m_adjList.resize(other.m_adjList.size());
    for (std::size_t i = 0; i < other.m_adjList.size(); ++i)
      for (std::size_t j = 0; j < other.m_adjList[i].size(); ++j)
        m_adjList[i].push_back(bond_type(this, other.m_adjList[i][j].index()));

    // copy atom properties
    m_atomAromatic = other.m_atomAromatic;
    m_element = other.m_element;
    m_mass = other.m_mass;
    m_hydrogens = other.m_hydrogens;
    m_charge = other.m_charge;

    // copy bond properties
    m_source = other.m_source;
    m_target = other.m_target;
    m_bondAromatic = other.m_bondAromatic;
    m_order = other.m_order;
  }

  HeMol& HeMol::operator=(const HeMol &other)
  {
    // copy m_adjList
    m_adjList.clear();
    m_adjList.resize(other.m_adjList.size());
    for (std::size_t i = 0; i < other.m_adjList.size(); ++i)
      for (std::size_t j = 0; j < other.m_adjList[i].size(); ++j)
        m_adjList[i].push_back(bond_type(this, other.m_adjList[i][j].index()));

    // copy atom properties
    m_atomAromatic = other.m_atomAromatic;
    m_element = other.m_element;
    m_mass = other.m_mass;
    m_hydrogens = other.m_hydrogens;
    m_charge = other.m_charge;

    // copy bond properties
    m_source = other.m_source;
    m_target = other.m_target;
    m_bondAromatic = other.m_bondAromatic;
    m_order = other.m_order;

    return *this;
  }

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

  namespace impl {

    struct SortBondsByDecreasingIndex
    {
      template<typename BondType>
      bool operator()(const BondType &bond1, const BondType &bond2) const
      {
        return bond1.index() > bond2.index();
      }
    };

  }

  void HeMol::removeAtom(const atom_type &atom)
  {
    Index index = atom.index();
    std::vector<bond_type> bonds = m_adjList[index];

    // remove properties
    m_adjList.erase(m_adjList.begin() + index);
    m_atomAromatic.erase(m_atomAromatic.begin() + index);
    m_element.erase(m_element.begin() + index);
    m_mass.erase(m_mass.begin() + index);
    m_hydrogens.erase(m_hydrogens.begin() + index);
    m_charge.erase(m_charge.begin() + index);

    // sort bonds by decreasing bond index so they can be correctly removed
    std::sort(bonds.begin(), bonds.end(), impl::SortBondsByDecreasingIndex());

    // update m_source & m_target
    for (std::size_t i = 0; i < m_source.size(); ++i)
      if (m_source[i] > index)
        --m_source[i];
    for (std::size_t i = 0; i < m_target.size(); ++i)
      if (m_target[i] > index)
        --m_target[i];

    // remove bonds to removed atom
    for (std::size_t i = 0; i < bonds.size(); ++i)
      removeBond(bonds[i]);
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

  void HeMol::removeBond(const bond_type &bond)
  {
    Index index = bond.index();

    // remove properties
    m_source.erase(m_source.begin() + index);
    m_target.erase(m_target.begin() + index);
    m_bondAromatic.erase(m_bondAromatic.begin() + index);
    m_order.erase(m_order.begin() + index);

    // update m_adjList
    for (std::size_t i = 0; i < m_adjList.size(); ++i) {
      int remove = -1;
      for (std::size_t j = 0; j < m_adjList[i].size(); ++j) {
        if (m_adjList[i][j].index() == index)
          remove = j;
        else if (m_adjList[i][j].index() > index)
          m_adjList[i][j] = bond_type(this, m_adjList[i][j].index() - 1);
      }
      if (remove != -1)
        m_adjList[i].erase(m_adjList[i].begin() + remove);
    }
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

  std::ostream& operator<<(std::ostream &os, HeMol &mol)
  {
    os << "Molecule:" << std::endl;
    os << "    Atoms:\tindex\telement" << std::endl;
    molecule_traits<HeMol>::atom_iter atom, end_atoms;
    TIE(atom, end_atoms) = mol.atoms();
    for (; atom != end_atoms; ++atom)
      os << "          \t" << (*atom).index() << "\t" << (*atom).element() << std::endl;

    os << "    Bonds:\tsource\ttarget\torder" << std::endl;
    molecule_traits<HeMol>::bond_iter bond, end_bonds;
    TIE(bond, end_bonds) = mol.bonds();
    for (; bond != end_bonds; ++bond)
      os << "          \t" << (*bond).source().index() << "\t" << (*bond).target().index() << "\t" << (*bond).order() << std::endl;
    return os;
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
