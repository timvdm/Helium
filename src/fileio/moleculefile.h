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
#ifndef HELIUM_FILEIO_MOLECULES_H
#define HELIUM_FILEIO_MOLECULES_H

#include <Helium/molecule.h>
#include <Helium/util.h>
#include <Helium/fileio/file.h>

#include <json/json.h>

#include <boost/iostreams/device/mapped_file.hpp>

namespace Helium {

  /**
   * @brief Read a molecule from a STL input stream.
   *
   * @param is The input stream.
   * @param mol The molecule.
   *
   * @return True if successful.
   */
  template<typename EditableMoleculeType>
  bool read_molecule(std::istream &is, EditableMoleculeType &mol)
  {
    typedef typename molecule_traits<EditableMoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<EditableMoleculeType>::bond_type bond_type;

    clear_molecule(mol);

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
      atom_type atom = add_atom(mol);
      // set the atom properties
      set_aromatic(mol, atom, aromatic);
      set_element(mol, atom, element);
      set_mass(mol, atom, mass);
      set_hydrogens(mol, atom, hydrogens);
      set_charge(mol, atom, charge);
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
      atom_type s = get_atom(mol, source);
      atom_type t = get_atom(mol, target);

      // create the bond
      bond_type bond = add_bond(mol, s, t);
      // set the bond properties
      set_aromatic(mol, bond, props & 128);
      set_order(mol, bond, props & 63);
    }

    return true;
  }

  /**
   * @brief Read a molecule from a memory location.
   *
   * @param data Pointer to the memory.
   * @param mol The molecule.
   *
   * @return True if successful.
   */
  template<typename EditableMoleculeType>
  bool read_molecule(const char *data, EditableMoleculeType &mol)
  {
    typedef typename molecule_traits<EditableMoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<EditableMoleculeType>::bond_type bond_type;

    clear_molecule(mol);

    // read number of atoms and bonds
    unsigned short numAtoms = *reinterpret_cast<const unsigned short*>(data);
    unsigned short numBonds = *reinterpret_cast<const unsigned short*>(data + sizeof(unsigned short));

    unsigned char element, aromatic, mass, hydrogens;
    signed char charge;
    for (int i = 0; i < numAtoms; ++i) {
      unsigned int offset = 2 * sizeof(unsigned short) + 6 * i;
      // read the element
      element = *reinterpret_cast<const unsigned char*>(data + offset);
      aromatic = *reinterpret_cast<const unsigned char*>(data + offset + 2);
      mass = *reinterpret_cast<const unsigned char*>(data + offset + 3);
      hydrogens = *reinterpret_cast<const unsigned char*>(data + offset + 4);
      charge = *reinterpret_cast<const signed char*>(data + offset + 5);

      // create the atom
      atom_type atom = add_atom(mol);
      // set the atom properties
      set_aromatic(mol, atom, aromatic);
      set_element(mol, atom, element);
      set_mass(mol, atom, mass);
      set_hydrogens(mol, atom, hydrogens);
      set_charge(mol, atom, charge);
    }

    unsigned short source, target;
    unsigned char props;
    for (int i = 0; i < numBonds; ++i) {
      unsigned int offset = 2 * sizeof(unsigned short) + 6 * numAtoms + 5 * i;
      // source and target indices
      source = *reinterpret_cast<const unsigned short*>(data + offset);
      target = *reinterpret_cast<const unsigned short*>(data + offset + sizeof(unsigned short));
      // read remaining properties: cyclic, aromatic & bond order
      props = *reinterpret_cast<const unsigned char*>(data + offset + 2 * sizeof(unsigned short));

      // get the atoms
      atom_type s = get_atom(mol, source);
      atom_type t = get_atom(mol, target);

      // create the bond
      bond_type bond = add_bond(mol, s, t);
      // set the bond properties
      set_aromatic(mol, bond, props & 128);
      set_order(mol, bond, props & 63);
    }

    return true;
  }

  /**
   * @brief Class for reading molecules from Helium binary file.
   *
   * When load() is called, this class reads the file header and loads the
   * position index for fast random access to molecules.
   */
  class MoleculeFile
  {
    public:
      /**
       * @brief Constructor.
       */
      MoleculeFile() : m_numMolecules(0), m_positionsPos(0)
      {
      }

      /**
       * @brief Constructor.
       *
       * @param filename The filename.
       */
      MoleculeFile(const std::string &filename)
      {
        load(filename);
      }

      /**
       * @brief Load a file.
       *
       * This function reads the file header and prepares the file for use.
       *
       * @param filename The filename.
       */
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
        m_positionsPos = data["molecule_indexes"].asUInt64();

        // read the molecule indexes
        m_positions.resize(m_numMolecules);
        m_file.stream().seekg(m_positionsPos);
        m_file.read(&m_positions[0], m_positions.size() * sizeof(uint64_t));

        // reset stream position to read first molecule
        m_file.seek(0);
      }

      /**
       * @brief Get the number of molecules.
       *
       * @return The number of molecules.
       */
      unsigned int numMolecules() const
      {
        return m_numMolecules;
      }

      /**
       * Read the next molecule from the file.
       *
       * @return True if successfull.
       */
      template<typename EditableMoleculeType>
      bool readMolecule(EditableMoleculeType &mol)
      {
        if (!(bool)m_file)
          return false;
        return Helium::read_molecule(m_file.stream(), mol);
      }

      /**
       * Read the molecule with the specified index from the file.
       *
       * @return True if successfull.
       */
      template<typename EditableMoleculeType>
      bool readMolecule(unsigned int index, EditableMoleculeType &mol)
      {
        if (!(bool)m_file)
          return false;
        m_file.stream().seekg(m_positions[index]);
        return Helium::read_molecule(m_file.stream(), mol);
      }

      /**
       * @brief Get the STL input stream.
       *
       * @return The STL input stream.
       */
      std::ifstream& stream()
      {
        return m_file.stream();
      }

      /**
       * @brief Get the stream positions of the molecules.
       *
       * @return The stream positions of the molecules.
       */
      std::vector<uint64_t> positions() const
      {
        return m_positions;
      }

      /**
       * @brief Get the stream position of the molecule position index.
       *
       * @return The stream position of the molecule position index.
       */
      Json::UInt64 positionsPos() const
      {
        return m_positionsPos;
      }

      /**
       * @brief Close the file.
       */
      void close()
      {
        m_file.close();
        m_positions.clear();
        m_numMolecules = 0;
        m_positionsPos = 0;
      }

    private:
      BinaryInputFile m_file; //!< The binary file.
      std::vector<uint64_t> m_positions; //!< The molecule positions.
      unsigned int m_numMolecules; //!< The number of molecules.
      Json::UInt64 m_positionsPos; //!< Position of molecule positions index.
  };

  /**
   * @brief Class for reading molecules from Helium binary file.
   */
  class MemoryMappedMoleculeFile
  {
    public:
      /**
       * @brief Constructor.
       */
      MemoryMappedMoleculeFile() : m_numMolecules(0)
      {
      }

      /**
       * @brief Constructor.
       *
       * @param filename The filename.
       */
      MemoryMappedMoleculeFile(const std::string &filename)
      {
        load(filename);
      }

      /**
       * @brief Load a file.
       *
       * This function reads the file header and prepares the file for use.
       *
       * @param filename The filename.
       */
      void load(const std::string &filename)
      {
        TIMER("MemoryMappedMoleculeFile::load():");

        // use a temporary BinaryInputFile to easily read the header
        BinaryInputFile file(filename);

        // read the josn header
        Json::Reader reader;
        Json::Value data;
        if (!reader.parse(file.header(), data))
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

        // open the memory mapped file
        m_mappedFile.open(filename);
        assert(m_mappedFile.is_open());
        assert(m_mappedFile.data());

        // read the molecule indexes
        m_positions.resize(m_numMolecules);
        for (unsigned int i = 0; i < m_positions.size(); ++i)
          m_positions[i] = *reinterpret_cast<const uint64_t*>(m_mappedFile.data() + positionsPos + i * sizeof(uint64_t));
      }

      /**
       * @brief Get the number of molecules.
       *
       * @return The number of molecules.
       */
      unsigned int numMolecules() const
      {
        return m_numMolecules;
      }

      /**
       * Read the molecule with the specified index from the file.
       *
       * @return True if successfull.
       */
      template<typename MoleculeType>
      bool readMolecule(unsigned int index, MoleculeType &mol)
      {
        Helium::read_molecule(m_mappedFile.data() + m_positions[index], mol);
        return true;
      }

    private:
      boost::iostreams::mapped_file_source m_mappedFile; //!< The memory mapped file.
      std::vector<uint64_t> m_positions; //!< The positions of the molecules in the file.
      unsigned int m_numMolecules; //!< The number of molecules in the file.
  };

}

#endif
