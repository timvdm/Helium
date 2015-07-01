/*
 * Copyright (c) 2014, Tim Vandermeersch
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
#ifndef HELIUM_KEKULIZE_H
#define HELIUM_KEKULIZE_H

#include <Helium/molecule.h>
#include <Helium/ring.h>
#include <Helium/algorithms/cycles.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

namespace Helium {

  /**
   * @file algorithms/kekulize.h
   * @brief Kekulization.
   */

  /**
   * @brief Kekulize a molecule.
   *
   * This algorithm changes aromatic bonds to alternating single and double bonds.
   * It uses Edmonds' maximum matching algorithm on the graph of aromatic bonds
   * with some bonds excluded.
   *
   * Exclude bonds around hetero atoms in 5-membered rings if the atom is:
   * - oxygen (e.g. furan)
   * - sulfur (e.g. thiophene)
   * - nitrogen with: charge = 0 and (1 implicit hydrogen or degree 3)
   *   (e.g. pyrrole)
   *
   * @param mol The molecule.
   * @param rings The ring set.
   *
   * @return True if successful.
   */
  template<typename EditableMoleculeType>
  bool kekulize(EditableMoleculeType &mol, const RingSet<EditableMoleculeType> &rings)
  {
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> graph_type;
    typedef typename molecule_traits<EditableMoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<EditableMoleculeType>::bond_type bond_type;

    // create a list of aromatic bonds
    std::vector<bool> bonds(num_bonds(mol));
    for (auto &bond : get_bonds(mol))
      if (is_aromatic(mol, bond))
        bonds[get_index(mol, bond)] = true;

    // exclude bonds around hetero atoms in 5-membered rings if the atom is:
    // - oxygen (e.g. furan)
    // - sulfur (e.g. thiophene)
    // - nitrogen with: charge = 0 and (1 implicit hydrogen or degree 3)
    //   (e.g. pyrrole)
    for (std::size_t i = 0; i < rings.size(); ++i) {
      const Ring<EditableMoleculeType> &ring = rings.ring(i);
      if (ring.size() != 5)
        continue;

      for (std::size_t j = 0; j < ring.size(); ++j) {
        atom_type atom = ring.atom(j);

        bool exclude = false;
        if (is_oxygen(mol, atom))
          exclude = true;
        else if (is_sulfur(mol, atom))
          exclude = true;
        else if (is_nitrogen(mol, atom) && get_charge(mol, atom) == 0 &&
            (get_degree(mol, atom) == 3 || get_hydrogens(mol, atom) == 1))
          exclude = true;

        if (exclude) {
          for (auto &bond : get_bonds(mol, atom))
            bonds[get_index(mol, bond)] = false;
        }
      }
    }

    // create boost graph
    graph_type g(num_atoms(mol));

    for (auto &bond : get_bonds(mol))
      if (bonds[get_index(mol, bond)])
        boost::add_edge(get_index(mol, get_source(mol, bond)), get_index(mol, get_target(mol, bond)), g);

    // run maximum matching algorithm
    std::vector<boost::graph_traits<graph_type>::vertex_descriptor> mate(num_atoms(mol));
    if (!boost::checked_edmonds_maximum_cardinality_matching(g, &mate[0]))
      return false;

    // mark all aromatic bonds as single
    for (auto &bond : get_bonds(mol))
      if (is_aromatic(mol, bond)) {
        set_order(mol, bond, 1);
        set_aromatic(mol, bond, false);
      }

    // mark all atoms ad non-aromatic
    for (auto &atom : get_atoms(mol))
      set_aromatic(mol, atom, false);

    // mark all bonds in the maximum matching as double
    boost::graph_traits<graph_type>::vertex_iterator vi, vi_end;
    for(boost::tie(vi,vi_end) = boost::vertices(g); vi != vi_end; ++vi)
      if (mate[*vi] != boost::graph_traits<graph_type>::null_vertex() && *vi < mate[*vi]) {
        bond_type bond = get_bond(mol, get_atom(mol, *vi), get_atom(mol, mate[*vi]));
        set_order(mol, bond, 2);
      }

    return true;
  }

  /**
   * @overload
   */
  template<typename EditableMoleculeType>
  bool kekulize(EditableMoleculeType &mol)
  {
    return kekulize(mol, relevant_cycles(mol));
  }

}

#endif
