/**
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
#ifndef HELIUM_CYCLES_H
#define HELIUM_CYCLES_H

#include "components.h"
#include "isomorphism.h"
#include "smiles.h"

namespace Helium {

  template<typename MoleculeType>
  Size cyclomatic_number(MoleculeType &mol, Size numComponents)
  {
    return num_bonds(mol) - num_atoms(mol) + numComponents;
  }

  template<typename MoleculeType>
  Size cyclomatic_number(MoleculeType &mol)
  {
    return cyclomatic_number(mol, num_connected_components(mol));
  }

  namespace impl {

    template<typename MoleculeType, typename QueryType>
    struct CycleAtomMatcher
    {
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<QueryType>::atom_type query_atom_type;

      bool operator()(QueryType &query, query_atom_type queryAtom, MoleculeType &mol, atom_type atom) const
      {
        return is_cyclic(mol, atom);
      }
    };

    template<typename MoleculeType, typename QueryType>
    struct CycleBondMatcher
    {
      typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
      typedef typename molecule_traits<QueryType>::bond_type query_bond_type;

      bool operator()(QueryType &query, query_bond_type queryAtom, MoleculeType &mol, bond_type bond) const
      {
        return is_cyclic(mol, bond);
      }
    };

  }

  template<typename MoleculeType>
  std::vector<std::vector<Index> > relevant_cycles(MoleculeType &mol, Size cyclomaticNumber)
  {
    std::vector<std::vector<Index> > cycles;

    // current cycle size
    unsigned int size = 3;

    while (cycles.size() < cyclomaticNumber) {
      // create query
      std::string smiles = "*1" + std::string(size - 1, '*') + "1";
      HeMol cycleMol;
      parse_smiles(smiles, cycleMol);

      // find all cycles of size
      MappingList mappings;
      if (isomorphism_search<impl::CycleAtomMatcher, impl::CycleBondMatcher, MoleculeType, HeMol, MappingList>(mol, cycleMol, mappings)) {
        for (std::size_t i = 0; i < mappings.maps.size(); ++i)
          cycles.push_back(mappings.maps[i]);
      }

      // increment cycle size
      ++size;
    }

    return cycles;
  }

  template<typename MoleculeType>
  std::vector<std::vector<Index> > relevant_cycles(MoleculeType &mol)
  {
    return relevant_cycles(mol, cyclomatic_number(mol));
  }

}

#endif
