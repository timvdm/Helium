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
#ifndef HELIUM_FILEIO_MOLECULEFILE_H
#define HELIUM_FILEIO_MOLECULEFILE_H

#include <Helium/molecule.h>
#include <Helium/util.h>
#include <Helium/fileio/file.h>
#include <Helium/json/json.h>

#include <boost/iostreams/device/mapped_file.hpp>

namespace Helium {

  /**
   * @file fileio/moleculefile.h
   * @brief Classes for reading and writing molecule files.
   */

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
   * @brief Write a molecule to an STL output stream.
   *
   * @param os The STL output stream.
   * @param mol The molecule to write.
   */
  template<typename MoleculeType>
  void write_molecule(std::ostream &os, const MoleculeType &mol)
  {
    // hydrogen atoms are not written to the file
    // create a map to map atom indices & count the number of heavy atoms
    std::vector<unsigned short> indices(num_atoms(mol));
    unsigned short numAtoms = 0;
    FOREACH_ATOM (atom, mol) {
      indices[get_index(mol, *atom)] = numAtoms;
      if (!is_hydrogen(mol, *atom))
        ++numAtoms;
    }

    // count the number of bonds between heavy atoms
    unsigned short numBonds = 0;
    FOREACH_BOND (bond, mol)
      if (!is_hydrogen(mol, get_source(mol, *bond)) && !is_hydrogen(mol, get_target(mol, *bond)))
        ++numBonds;

    // write the number of atoms & bonds
    write16<unsigned short>(os, numAtoms);
    write16<unsigned short>(os, numBonds);

    // write atoms (6 byte / atom)
    FOREACH_ATOM (atom, mol) {
      if (is_hydrogen(mol, *atom))
        continue;

      // write the element
      write8<unsigned char>(os, get_element(mol, *atom));
      // write cyclic property
      write8<unsigned char>(os, 0);
      // write aromatic property
      write8<unsigned char>(os, is_aromatic(mol, *atom));
      // write isotope
      write8<unsigned char>(os, get_mass(mol, *atom));
      // write hydrogen count
      write8<unsigned char>(os, get_hydrogens(mol, *atom) + get_degree(mol, *atom) - get_heavy_degree(mol, *atom));
      // write formal charge
      write8<signed char>(os, get_charge(mol, *atom));
    }
  
    // write bonds (5 byte / bond
    FOREACH_BOND (bond, mol) {
      if (is_hydrogen(mol, get_source(mol, *bond)) || is_hydrogen(mol, get_target(mol, *bond)))
        continue;

      // write source & target indices
      write16<unsigned short>(os, indices[get_index(mol, get_source(mol, *bond))]);
      write16<unsigned short>(os, indices[get_index(mol, get_target(mol, *bond))]);
      // write bond order + aromatic & cyclic properties
      unsigned char props = get_order(mol, *bond);
      if (is_aromatic(mol, *bond))
        props |= 128;
      //if (bond->IsInRing())
      //  props |= 64;
      write8<unsigned char>(os, props);
    }
  }

  /**
   * @class MoleculeOutputFile fileio/moleculefile.h <Helium/fileio/moleculefile.h>
   * @brief Class for writing molecule files.
   */
  class MoleculeOutputFile
  {
    public:
      /**
       * @brief Constructor.
       *
       * @param filename The filename.
       */
      MoleculeOutputFile(const std::string &filename) : m_file(filename)
      {
      }

      ~MoleculeOutputFile()
      {
        // save the stream position where the molecule positions are stored
        Json::UInt64 positionsPos = m_file.stream().tellp();

        // write the molecule positions to the file
        m_file.write(&m_positions[0], m_positions.size() * sizeof(std::size_t));

        // create JSON header
        Json::Value data;
        data["filetype"] = "molecules";
        data["num_molecules"] = Json::UInt64(m_positions.size());
        data["molecule_indexes"] = positionsPos;

        // write JSON header
        Json::StyledWriter writer;
        m_file.writeHeader(writer.write(data));
      }

      /**
       * @brief Write a molecule to the file.
       *
       * @param mol The molecule to write.
       */
      template<typename MoleculeType>
      bool writeMolecule(const MoleculeType &mol)
      {
        m_positions.push_back(m_file.stream().tellp());
        write_molecule(m_file.stream(), mol);
        return true;
      }

      /**
       * @brief Get the last error.
       *
       * @return The last error.
       */
      const Error& error() const
      {
        return m_file.error();
      }

    private:
      BinaryOutputFile m_file; //!< The binary output file.
      std::vector<std::size_t> m_positions; //!< The molecule positions in the file.
      Error m_error; //!< The last error.
  };

  /**
   * @class MoleculeFile fileio/moleculefile.h <Helium/fileio/moleculefile.h>
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
      bool load(const std::string &filename)
      {
        TIMER("MoleculeFile::load():");

        // reset error
        m_error = Error();

        if (!m_file.open(filename)) {
          m_error = m_file.error();
          return false;
        }

        // read the josn header
        Json::Reader reader;
        Json::Value data;
        if (!reader.parse(m_file.header(), data)) {
          m_error = Error(reader.getFormattedErrorMessages());
          return false;
        }

        // make sure the required attributes are present
        if (!data.isMember("filetype") || data["filetype"].asString() != "molecules") {
          m_error = Error(make_string("JSON header for file ", filename, " does not contain 'filetype' attribute or is not 'molecules'"));
          return false;
        }
        if (!data.isMember("num_molecules")) {
          m_error = Error(make_string("JSON header for file ", filename, " does not contain 'num_molecules' attribute"));
          return false;
        }
        if (!data.isMember("molecule_indexes")) {
          m_error = Error(make_string("JSON header for file ", filename, " does not contain 'molecule_indexes' attribute"));
          return false;
        }

        // extract needed attributes
        m_numMolecules = data["num_molecules"].asUInt();
        m_positionsPos = data["molecule_indexes"].asUInt64();

        // read the molecule indexes
        m_positions.resize(m_numMolecules);
        m_file.stream().seekg(m_positionsPos);
        m_file.read(&m_positions[0], m_positions.size() * sizeof(uint64_t));

        // reset stream position to read first molecule
        m_file.seek(0);

        return true;
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

      /**
       * @brief Get the last error.
       *
       * @return The last error.
       */
      const Error& error() const
      {
        return m_error;
      }

    private:
      BinaryInputFile m_file; //!< The binary file.
      std::vector<uint64_t> m_positions; //!< The molecule positions.
      unsigned int m_numMolecules; //!< The number of molecules.
      Json::UInt64 m_positionsPos; //!< Position of molecule positions index.
      Error m_error; //!< The last error.
  };

  /**
   * @class MemoryMappedMoleculeFile fileio/moleculefile.h <Helium/fileio/moleculefile.h>
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
      bool load(const std::string &filename)
      {
        TIMER("MemoryMappedMoleculeFile::load():");

        // reset error
        m_error = Error();

        // use a temporary BinaryInputFile to easily read the header
        BinaryInputFile file(filename);
        if (file.error()) {
          m_error = file.error();
          return false;
        }

        // read the josn header
        Json::Reader reader;
        Json::Value data;
        if (!reader.parse(file.header(), data)) {
          m_error = Error(reader.getFormattedErrorMessages());
          return false;
        }

        // make sure the required attributes are present
        if (!data.isMember("filetype") || data["filetype"].asString() != "molecules") {
          m_error = Error(make_string("JSON header for file ", filename, " does not contain 'filetype' attribute or is not 'molecules'"));
          return false;
        }
        if (!data.isMember("num_molecules")) {
          m_error = Error(make_string("JSON header for file ", filename, " does not contain 'num_molecules' attribute"));
          return false;
        }
        if (!data.isMember("molecule_indexes")) {
          m_error = Error(make_string("JSON header for file ", filename, " does not contain 'molecule_indexes' attribute"));
          return false;
        }

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

        return true;
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

      /**
       * @brief Get the last error.
       *
       * @return The last error.
       */
      const Error& error() const
      {
        return m_error;
      }

    private:
      boost::iostreams::mapped_file_source m_mappedFile; //!< The memory mapped file.
      std::vector<uint64_t> m_positions; //!< The positions of the molecules in the file.
      unsigned int m_numMolecules; //!< The number of molecules in the file.
      Error m_error; //!< The last error.
  };

}

#endif
