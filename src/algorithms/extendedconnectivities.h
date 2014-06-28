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
#ifndef HELIUM_EXTENDEDCONNECTIVITIES_H
#define HELIUM_EXTENDEDCONNECTIVITIES_H

#include <Helium/algorithms/invariants.h>
#include <Helium/util/vector.h>
#include <set>

namespace Helium {

  /**
   * @file algorithms/extendedconnectivities.h
   * @brief Morgan's extended-connectivities algorithm.
   */

  namespace impl {

    /**
     * Perform a single Extended Connectivities (EC) iteration.
     */
    template<typename MoleculeType>
    void extended_connectivities_iterate(const MoleculeType &mol, std::vector<unsigned long> &ec)
    {
      std::vector<unsigned long> next = ec;

      FOREACH_ATOM (atom, mol)
        FOREACH_NBR (nbr, *atom, mol)
          next[get_index(mol, *atom)] += ec[get_index(mol, *nbr)];

      ec.swap(next);
    }

    /**
     * Renumber the extended connectivities values to be in the range [0,n-1]
     * where n is the number of unique values.
     */
    inline void extended_connectivities_renumber(std::vector<unsigned long> &ec)
    {
      std::set<unsigned long> classes;
      for (std::size_t i = 0; i < ec.size(); ++i)
        classes.insert(ec[i]);
      unsigned long cls = 0;
      for (std::set<unsigned long>::iterator i = classes.begin(); i != classes.end(); ++i) {
        for (std::size_t j = 0; j < ec.size(); ++j)
          if (ec[j] == *i)
            ec[j] = cls;
        ++cls;
      }
    }

    inline unsigned int num_symmetry_classes(std::vector<unsigned long> symmetry)
    {
      std::sort(symmetry.begin(), symmetry.end());
      unsigned int result = 0;
      for (std::size_t i = 1; i < symmetry.size(); ++i)
        if (symmetry[i - 1] != symmetry[i])
          result++;
      return result;
    }

  }

  /**
   * Calculate the Morgan's Extended Connectivities for the specified molecule.
   *
   * @param mol The molecule.
   * @param atomInvariant The atom invariant to use.
   *
   @verbatim
   Morgan, H. L. The Generation of a Unique Machine Description for Chemical
   Structures - A Technique Developed at Chemical Abstracts Service. J. Chem.
   Doc. 1965, 5: 107-112.
   @endverbatim
   *
   * @return The extended connectivities.
   */
  template<typename MoleculeType, typename AtomInvariant>
  std::vector<unsigned long> extended_connectivities(const MoleculeType &mol,
      const AtomInvariant &atomInvariant)
  {
    // initial atom invariants
    std::vector<unsigned long> ec;
    FOREACH_ATOM (atom, mol)
      ec.push_back(atomInvariant(mol, *atom));

    // iterate
    unsigned int numClasses = impl::num_symmetry_classes(ec);
    for (int i = 0; i < 100; ++i) { // should never reach 100...
      impl::extended_connectivities_iterate(mol, ec);
      //unsigned int nextNumClasses = unique_elements(ec);
      unsigned int nextNumClasses = impl::num_symmetry_classes(ec);
      // if the number of unique values didn't change, stop the iteration
      if (numClasses == nextNumClasses)
        break;
      numClasses = nextNumClasses;
    }

    // renumber the EC values
    impl::extended_connectivities_renumber(ec);
    return ec;
  }

}

#endif
