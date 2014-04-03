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
#ifndef HELIUM_INVARIANTS_H
#define HELIUM_INVARIANTS_H

#include <Helium/tie.h>
#include <Helium/molecule.h>

namespace Helium {

  /**
   * @brief Atom invariant to be used for canonical coding.
   *
   * @note Complexity: @f$O(1)@f$
   * @ingroup Production
   * @note Phase: Production
   */
  class DefaultAtomInvariant
  {
    public:
      /**
       * @brief Possible atom invariants.
       */
      enum Invariants {
        /**
         * @brief Include none.
         */
        None = 0,
        /**
         * @brief Include atom element.
         */
        Element = 1,
        /**
         * @brief Include atom mass.
         */
        Mass = 2,
        /**
         * @brief Include atom charge.
         */
        Charge = 4,
        /**
         * @brief Include atom degree.
         */
        Degree = 8,
        /**
         * @brief Include atom aromaticity.
         */
        Aromatic = 16,
        /**
         * @brief Include all invariants.
         */
        All = Element | Mass | Charge | Degree | Aromatic
      };

      /**
       * @brief Constructor.
       *
       * @param invariants The atom attributes to include in the invariants.
       */
      DefaultAtomInvariant(int invariants = All) : m_invariants(invariants)
      {
      }

      /**
       * @brief Get an invariant label classifying the atom.
       *
       * This invariant label is used for generating and comparing canonical
       * codes.
       *
       * @param mol The molecule.
       * @param atom The atom.
       *
       * @return The invariant label.
       */
      template<typename MoleculeType, typename AtomType>
      unsigned long operator()(const MoleculeType &mol, AtomType atom) const
      {
        unsigned long invariant = 0;
        if (m_invariants & Element)
          invariant += get_element(mol, atom);
        if (m_invariants & Mass)
          invariant += get_mass(mol, atom) * 1000;
        if (m_invariants & Charge)
          invariant += (10 + get_charge(mol, atom)) * 1000000;
        if (m_invariants & Degree)
          invariant += get_degree(mol, atom) * 100000000;
        if (m_invariants & Aromatic)
          invariant += is_aromatic(mol, atom) * 1000000000;
        return invariant;
      }

    private:
      int m_invariants; //!< The invariants.
  };

  /**
   * @brief Bond invariant to be used for canonical coding.
   *
   * @note Complexity: @f$O(1)@f$
   * @ingroup Production
   * @note Phase: Production
   */
  class DefaultBondInvariant
  {
    public:
      /**
       * @brief Possible bond invariants.
       */
      enum Invariants {
        /**
         * @brief Include none.
         */
        None = 0,
        /**
         * @brief Include bond order.
         */
        Order = 1,
        /**
         * @brief Include bond aromaticity.
         */
        Aromatic = 2,
        /**
         * @brief Include all invariants.
         */
        All = Order | Aromatic
      };

      /**
       * @brief Constructor.
       *
       * @param invariants The atom attributes to include in the invariants.
       */
      DefaultBondInvariant(int invariants = All) : m_invariants(invariants)
      {
      }

      /**
       * @brief Get an invariant label classifying the bond.
       *
       * This invariant label is used for generating and comparing canonical
       * codes.
       *
       * @param mol The molecule.
       * @param bond The bond.
       *
       * @return The invariant label.
       */
      template<typename MoleculeType, typename BondType>
      unsigned long operator()(const MoleculeType &mol, BondType bond) const
      {
        unsigned long invariant = 0;
        if (m_invariants & Order)
          invariant += get_order(mol, bond);
        if (m_invariants & Aromatic)
          invariant += is_aromatic(mol, bond) * 100;
        return invariant;
      }

    private:
      int m_invariants; //!< The invariants.
  };

}

#endif
