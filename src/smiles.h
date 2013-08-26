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
#ifndef HELIUM_SMILES_H
#define HELIUM_SMILES_H

#include <Helium/hemol.h>

#include <smiley.h>

#include <sstream>


namespace Helium {

  namespace impl {

    template<typename MoleculeType>
    struct SmileyCallback : public Smiley::CallbackBase
    {
      SmileyCallback(MoleculeType &mol_) : mol(mol_)
      {
      }

      void clear()
      {
        mol.clear();
      }

      void addAtom(int element, bool aromatic, int isotope, int hCount, int charge, int atomClass)
      {
        typename molecule_traits<MoleculeType>::atom_type atom = mol.addAtom();
        atom.setElement(element);
        atom.setAromatic(aromatic);
        atom.setMass(isotope);
        atom.setHydrogens(hCount);
        atom.setCharge(charge);
      }

      void addBond(int source, int target, int order, bool isUp, bool isDown)
      {
        typename molecule_traits<MoleculeType>::bond_type bond = mol.addBond(mol.atom(source), mol.atom(target));
        if (order == 5)
          bond.setAromatic(true);
        bond.setOrder(order);
      }

      MoleculeType &mol;
    };

  }

  /**
   * @brief Parse a SMILES string.
   *
   * @tparam MoleculeType The type of the molecule, this must be a model of the
   *         Molecule concept.
   *
   * @param smiles The SMILES string.
   * @param mol The molecule.
   */
  template<typename MoleculeType>
  void parse_smiles(const std::string &smiles, MoleculeType &mol)
  {
    impl::SmileyCallback<MoleculeType> callback(mol);
    Smiley::Parser<impl::SmileyCallback<MoleculeType> > parser(callback);

    //parser.parse(smiles);
    try {
      parser.parse(smiles);
    } catch (Smiley::Exception &e) {
      std::ostringstream errorStream;
      if (e.type() == Smiley::Exception::SyntaxError)
        errorStream << "Syntax";
      else
        errorStream << "Semantics";
      errorStream << "Error: " << e.what() << "." << std::endl;
      errorStream << smiles << std::endl;
      for (std::size_t i = 0; i < e.pos(); ++i)
        errorStream << " ";
      for (std::size_t i = 0; i < e.length(); ++i)
        errorStream << "^";
      errorStream << std::endl;

      // Throw new exception with added details ...
      throw Smiley::Exception(e.type(), e.errorCode(),
          errorStream.str(), e.pos(), e.length());

    }
  }

}

#endif
