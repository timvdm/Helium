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
#ifndef HELIUM_AROMATIZE_H
#define HELIUM_AROMATIZE_H

#include <Helium/molecule.h>
#include <Helium/ring.h>
#include <Helium/algorithms/cycles.h>

namespace Helium {

  /**
   * @file algorithms/aromatize.h
   * @brief Aromatize.
   */

  /**
   * @brief Aromatize a molecule.
   *
   * This algorithm changes alternating single and double bond paths to aromatic
   * paths using Huckel's rule.
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
  bool aromatize(EditableMoleculeType &mol, const RingSet<EditableMoleculeType> &rings)
  {
    typedef typename molecule_traits<EditableMoleculeType>::atom_type atom_type;

    for (std::size_t i = 0; i < rings.size(); ++i) {
      const Ring<EditableMoleculeType> &ring = rings.ring(i);

      int electrons = 0;
      for (std::size_t j = 0; j < ring.size(); ++j) {
        atom_type atom = ring.atom(j);

        bool hasDoubleBond = false;
        FOREACH_INCIDENT (bond, atom, mol)
          if (get_order(mol, *bond) == 2) {
            hasDoubleBond = true;
            break;
          }

        switch (get_element(mol, atom)) {
          case 6:
            if (hasDoubleBond)
              electrons++;
            break;
          case 7:
            if (ring.size() == 6 && hasDoubleBond)
              electrons++;
            else if (ring.size() == 5 && !hasDoubleBond)
              electrons += 2;
            break;
          case 8:
          case 16:
            if (ring.size() == 5 && !hasDoubleBond)
              electrons += 2;
            break;
        }
      }

      bool aromatic = ((electrons % 4) == 2);

      if (aromatic) {
        for (std::size_t j = 0; j < ring.size(); ++j) {
          set_aromatic(mol, ring.atom(j), true);
          set_aromatic(mol, ring.bond(j), true);
        }
      }
    }

    return true;
  }

  /**
   * @overload
   */
  template<typename EditableMoleculeType>
  bool aromatize(EditableMoleculeType &mol)
  {
    return aromatize(mol, relevant_cycles(mol));
  }

}

#endif
