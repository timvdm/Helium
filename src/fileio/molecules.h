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
#ifndef HELIUM_FILEIO_MOLECULES_H
#define HELIUM_FILEIO_MOLECULES_H

#include "../hemol.h"
#include "../util.h"
#include "file.h"

#include <json/json.h>

namespace Helium {

  template<typename MoleculeType>
  bool read_molecule(std::istream &is, MoleculeType &mol)
  {
    mol.clear();

    if (!is)
      return false;

    // read number of atoms and bonds
    unsigned short numAtoms, numBonds;
    read16(is, numAtoms);
    read16(is, numBonds);

    if (!is)
      return false;

    unsigned char element, cyclic, aromatic, mass, hydrogens;
    signed char charge;
    for (int i = 0; i < numAtoms; ++i) {
      // read the element
      read8(is, element);
      read8(is, cyclic);
      read8(is, aromatic);
      read8(is, mass);
      read8(is, hydrogens);
      read8(is, charge);

      // create the atom
      HeAtom atom = mol.addAtom();
      // set the atom properties
      atom.setAromatic(aromatic);
      atom.setCyclic(cyclic);
      atom.setElement(element);
      atom.setMass(mass);
      atom.setHydrogens(hydrogens);
      atom.setCharge(charge);
    }

    unsigned short source, target;
    unsigned char props;
    for (int i = 0; i < numBonds; ++i) {
      // source and target indices
      read16(is, source);
      read16(is, target);
      // read remaining properties: cyclic, aromatic & bond order
      read8(is, props);

      // get the atoms
      HeAtom s = mol.atom(source);
      HeAtom t = mol.atom(target);

      // create the bond
      HeBond bond = mol.addBond(s, t);
      // set the bond properties
      bond.setAromatic(props & 128);
      bond.setCyclic(props & 64);
      bond.setOrder(props & 63);
    }

    return (bool)is;
  }

  class MoleculeFile
  {
    public:
      MoleculeFile() : m_numMolecules(0)
      {
      }

      MoleculeFile(const std::string &filename)
      {
        load(filename);
      }

      void load(const std::string &filename)
      {
        TIMER("MoleculeFile::load():");

        m_file.open(filename);

        // read the josn header
        Json::Reader reader;
        Json::Value data;
        if (!reader.parse(m_file.header(), data))
          throw std::runtime_error(reader.getFormattedErrorMessages());

        // make sure the required attributes are present
        if (!data.isMember("filetype") || data["filetype"].asString() != "molecules")
          throw std::runtime_error(make_string("JSON header for file ", filename, " does not contain 'filetype' attribute or is not 'molecules'"));
        if (!data.isMember("num_molecules"))
          throw std::runtime_error(make_string("JSON header for file ", filename, " does not contain 'num_molecules' attribute"));
        if (!data.isMember("molecule_indexes"))
          throw std::runtime_error(make_string("JSON header for file ", filename, " does not contain 'molecule_indexes' attribute"));

        // extract needed attributes
        m_numMolecules = data["num_molecules"].asUInt();
        Json::UInt64 positionsPos = data["molecule_indexes"].asUInt64();

        // read the molecule indexes
        m_positions.resize(m_numMolecules);
        m_file.stream().seekg(positionsPos);
        m_file.read(&m_positions[0], m_positions.size() * sizeof(std::size_t));

        // reset stream position to read first molecule
        m_file.seek(0);
      }

      unsigned int numMolecules() const
      {
        return m_numMolecules;
      }

      /**
       * Read the next molecule from the file.
       *
       * @return True if successfull.
       */
      template<typename MoleculeType>
      bool read_molecule(MoleculeType &mol)
      {
        if (!(bool)m_file)
          return false;
        Helium::read_molecule(m_file.stream(), mol);
        return true;
      }

      /**
       * Read the molecule with the specified index from the file.
       *
       * @return True if successfull.
       */
      template<typename MoleculeType>
      bool read_molecule(unsigned int index, MoleculeType &mol)
      {
        if (!(bool)m_file)
          return false;
        m_file.stream().seekg(m_positions[index]);
        Helium::read_molecule(m_file.stream(), mol);
        return true;
      }

      std::ifstream& stream()
      {
        return m_file.stream();
      }


    private:
      BinaryInputFile m_file;
      std::vector<std::size_t> m_positions;
      unsigned int m_numMolecules;
  };

  template<typename MoleculeType>
  void write_sdf(std::ostream &os, MoleculeType &mol)
  {
    std::streamsize width = os.width();
    os << std::endl;
    os << " Helium 0.0.1" << std::endl;
    os << std::endl;
    os.width(3);
    os << num_atoms(mol);
    os.width(3);
    os << num_bonds(mol) << "  0  0  0  0  0  0  0  0999 V2000" << std::endl;
  }

}

#endif
