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
#ifndef HELIUM_FILEIO_FINGERPRINTS_H
#define HELIUM_FILEIO_FINGERPRINTS_H

#include <Helium/bitvec.h>
#include <Helium/fileio/file.h>

#include <json/json.h>

#include <stdexcept>

namespace Helium {

  /**
   * @page fingerprints_page Fingerprint Storage and Indexes
   *
   * @section fingerprints_files Fingerprint File Classes
   *
   * To create Helium binary fingerprint files, the RowMajorFingerprintOutputFile
   * and ColumnMajorFingerprintOutputFile classes can be used. The difference
   * between row-major and column-major order files is explained below.
   *
   * For reading fingerprints, the InMemoryRowMajorFingerprintStorage and
   * InMemoryColumnMajorFingerprintStorage classes can be used. The classes
   * load and keep the fingerprint data in main memory. As an alternative,
   * the memory mapped equivalents MemoryMappedRowMajorFingerprintStorage and
   * MemoryMappedColumnMajorFingerprintStorage classes can be used. Additional
   * fingerprint storage classes can be added by following the fingerprint
   * storage concept described below.
   *
   * @section fingerprints_file_format Binary File Format
   *
   * Like all Helium binary files, the fingerprint file formats include a JSON
   * header containing all information about the binary data contained in the
   * file. Below is an example of such a header.
   @verbatim
   {
     'filetype': 'fingerprints',
     'order': 'row-major',
     'num_bits': 1024,
     'num_fingerprints': 1000000,
     'fingerprint': {
     'name': 'My custom fingerprint',
       'type': 'Helium::trees_fingerprint',
       'k': 7,
       'prime': 1021
     }
   }
   @endverbatim
   *
   * For fingerprint files the 'filetype' attribute will always be 'fingerprints'.
   * The 'order', 'num_bits', 'num_fingerprints' and 'fingerprint' attributes are
   * mandatory.
   * In the 'fingerprint' attribute, only the 'name' and 'type' attributes are
   * mandatory. Additional attributes may always be added such as the 'k' and
   * 'prime' attributes in the example above. (note: These last two attributes
   * may be required by an application that needs to compute the fingerprint for
   * a query molecule)
   *
   * While the 'name' attribute is mainly meant for end user purposes, the 'type'
   * attribute should be a string that corresponds to the method used for
   * generating the fingerprints (e.g. 'MyChemLib::some_fingerprint'). The
   * remaining attributes may then be used to specify parameters for this
   * fingerprint type.
   *
   * @section fingerprints_rowcol Row-Major vs. Column-Major Order
   *
   * Fingerprints can be stored in row-major order (i.e. each fingerprint is
   * stored in consecutive memory) or in column-major order (i.e. each
   * fingerprint bit is stored consecutively in memory). To illustrate the
   * difference, consider storing these 8-bit fingerprints for four molecules:
   *
   @verbatim
   molecule 1: 00101100
   molecule 2: 10001011
   molecule 3: 00000110
   molecule 4: 10101100
   @endverbatim
   *
   * In the row-major order storage, the fingerprints are stored by concatenating
   * the above fingerprints. In the figure below, the addresses are bit addresses
   * in hexadecimal. Using this storage method, the bits of an individual
   * molecule's fingerprint are stored consecutively.
   @verbatim
   0x00  00101100
   0x08  10001011
   0x10  00000110
   0x18  10101100
   @endverbatim
   *
   * In the column major order, the fingerprint are stored by concatenating the
   * columns. As a result, the individual fingerprint bits for all molecules are
   * stored consecutively. This can be seen in the figure below.
   *
   @verbatim
   0x00  0101
   0x04  0000
   0x08  1001
   0x0C  0000
   0x10  1101
   0x14  1011
   0x18  0110
   0x1C  0100
   @endverbatim
   *
   * The difference between these two methods of storing fingerprints can
   * dramatically influence the performance of an application. Due to various
   * hardware properties, sequential memory accesses are much faster than random
   * memory access patterns. Therefore it is vital to choose the correct storage
   * methods for each application.
   *
   * @section fingerprints_similarity Storing Fingerprints for Similarity Searches
   *
   * Similarity searches are performed by computing the similarity between the
   * query and queried molecules fingerprints. Although a variety of similarity
   * measures can be used the Tanimoto coefficient is very the most widely used.
   *
   * \f[
   *   T_{\mathrm{sim}} = \frac{| A \wedge B |}{| A \vee B |}
   * \f]
   *
   * Here \f$A\f$ and \f$B\f$ are the query and queried fingerprint bitstrings
   * and \f$|A|\f$ is the population count of fingerprint \f$A\f$.
   *
   * Since the whole fingerprint is needed to compute the similarity,
   * fingerprints for similarity searches are stored in row-major order.
   *
   * @section fingerprints_substructure Storing Fingerprints for Substructure Searches
   *
   * When performing a substructure search the fingerprints are used to filter out
   * molecules that will never match the query. All molecules that do not have
   * all bits set in their fingerprint that are set in the query's fingerprint
   * can never match the query. If the fingerprints (i.e. bitstrings) are
   * interpreted as sets where the value of the i-th bit would signify whether
   * the element i is in the set or not, this filtering can be seen as checking
   * if the query fingerprint is a subset of the molecule's fingerprint (which
   * would than be the superset).
   *
   * Let \f$M\f$ be set of all queried molecules\f$m\f$, \f$F\f$ be the set of all
   * possible fingerprints and \f$\mathrm{finger}: M \mapsto F\f$ be the function
   * assigning a fingerprint to a molecule. Then the set \f$I\f$ of molecules that
   * potentially contain the query is given by
   *
   * \f[
   *   I = \{ \forall m \in M \vert \mathrm{finger}(q) \subseteq \mathrm{finger}(m) \}
   * \f]
   *
   * where \f$q\f$ is the query.
   *
   * Constructing the set \f$I\f$ can easily be done by iterating over all
   * molecules (i.e. stored fingerprints) and adding them to \f$I\f$ if the
   * subset test passes. Although this would work, it is not very efficient.
   * A large portion of the data has to be processed that is not strictly
   * needed. A better way is to only check the bits that are set in the
   * query's fingerprint. While this is better in theory, care has to be taken
   * to implement this efficiently too. Checking only a single fingerprint bit
   * at a time using the row-major order storage would be very inefficient due
   * to random access memory patterns. To address this issue the column-major
   * order storage is used for substructure searching. It should also be noted
   * that constructing the final set \f$I\f$ can be done using the bitwise AND
   * operator on the bitstrings for the fingerprint bits that are set in the
   * query.
   *
   * A concrete example will illustrate how this works:
   *
   @verbatim
   query: 00001011 (bit 5, 7 & 8 are set)

   0x00  0101
   0x04  0000
   0x08  1001
   0x0C  0000
   0x10  1101 <- bit 5
   0x14  1011
   0x18  0110 <- bit 7
   0x1C  0100 <- bit 8

   1101 & 0110 & 0100 = 0100 -> only molecule 2 can contain the query
   @endverbatim
   *
   * In this example only 3 bitstrings have to be accessed. For this 8-bit
   * fingerprint this results in 62.5% less data that has to be processed.
   * Using a real fingerprint this may easily exceed 90% of the data (e.g.
   * a 1024-bit query fingerprint has 100 bits set).
   *
   * @section fingerprints_storage Fingerprint Storage Concepts
   *
   * In order for the indexing and search algorithms to work with various kind
   * of ways to actually store the fingerprints (e.g. read from disk on demand,
   * read from file and keep in main memory, memory mapped file, read from
   * database, ...), a concept for each major order is defined.
   *
   * @subsection fingerprints_storage_rowmajor Row-Major Order
   *
   * A class that is a model of the RowMajorFingerprintStorageConcept must support
   * the following operations:
   *
   @code
   std::string json = storage->header();
   unsigned int n = storage->numBits();
   unsigned int n = storage->numFingerprints();
   Helium::Word *fingerprint = storage->fingerprint(index);
   @endcode
   *
   * In the code above, storage is a pointer to an instance of a type that is a
   * model of the concept and index is an unsigned int. The index refers to the
   * molecule for which to get the fingerprint bitstring and should be in the
   * range [0, num_fingerprints). The returned pointer should point to a memory
   * location containing (for example) 1024 bits for a 1024-bit fingerprint.
   *
   * @subsection fingerprints_storage_colmajor Column-Major Order
   *
   * A class that is a model of the ColumnMajorFingerprintStorageConcept must support
   * the following operations:
   *
   @code
   std::string json = storage->header();
   unsigned int n = storage->numBits();
   unsigned int n = storage->numFingerprints();
   Helium::Word *bit = storage->bit(index);
   @endcode
   *
   * In the code above, storage is a pointer to an instance of a type that is a
   * model of the concept and index is an unsigned int. The index refers to the
   * bit in the fingerprint (e.g. in the range [0,1023] for a 1024-bit
   * fingerprint). The returned pointer should point to a memory location
   * containing num_fingerprints bits.
   */



  /**
   * @brief Output file for storing fingerprints in row-major order.
   */
  class RowMajorFingerprintOutputFile
  {
    public:
      /**
       * Constructor.
       *
       * @param filename The ouput filename.
       * @param numBits The number of bits in the fingerprints (e.g. 1024).
       */
      RowMajorFingerprintOutputFile(const std::string &filename, unsigned int numBits) : m_file(filename)
      {
        m_numBytes = bitvec_num_words_for_bits(numBits) * sizeof(Word);
      }

      /**
       * Write a single fingerprint to the file.
       *
       * @param fingerprint Pointer to the fingerprint.
       *
       * @return True if the fingerprint was successfully written to the file.
       */
      bool writeFingerprint(Word *fingerprint)
      {
        return m_file.write(fingerprint, m_numBytes);
      }

      /**
       * Write the JSON header to the file.
       *
       * @param header The JSON header.
       *
       * @return True if the header was successfully written to the file.
       */
      bool writeHeader(const std::string &header)
      {
        return m_file.writeHeader(header);
      }

    private:
      BinaryOutputFile m_file; //!< The output file.
      unsigned int m_numBytes; //!< The number of bytes in a fingerprint.
  };

  /**
   * @brief Output file for storing fingerprints in column-major order.
   */
  class ColumnMajorFingerprintOutputFile
  {
    public:
      /**
       * Constructor.
       *
       * @param filename The output filename.
       * @param numBits The number of bits in the fingerprint (e.g. 1024).
       * @param numFingerprints The number of fingerprints that will be written to the file.
       */
      ColumnMajorFingerprintOutputFile(const std::string &filename, unsigned int numBits,
          unsigned int numFingerprints) : m_file(filename), m_numBits(numBits),
          m_numFingerprints(numFingerprints), m_current(0)
      {
        // allocate data
        m_data = new Word[bitvec_num_words_for_bits(numFingerprints) * numBits];
        bitvec_zero(m_data, bitvec_num_words_for_bits(numFingerprints) * numBits);
      }

      /**
       * Destructor.
       */
      ~ColumnMajorFingerprintOutputFile()
      {
        // write the data
        m_file.write(m_data, bitvec_num_words_for_bits(m_numFingerprints) * m_numBits * sizeof(Word));
        delete [] m_data;
      }

      /**
       * Write a single fingerprint to the file.
       *
       * @param fingerprint Pointer to the fingerprint.
       *
       * @return True if the fingerprint was successfully written to the file.
       */
      bool writeFingerprint(Word *fingerprint)
      {
        // check each bit in the fingerprint
        for (int i = 0; i < m_numBits; ++i) {
          // skip this bit if it is not set
          if (!bitvec_get(i, fingerprint))
            continue;

          // set the correct bit
          bitvec_set(m_current, m_data + i * bitvec_num_words_for_bits(m_numFingerprints));
        }
        ++m_current;

        return true;
      }

      /**
       * Write the JSON header to the file.
       *
       * @param header The JSON header.
       *
       * @return True if the header was successfully written to the file.
       */
      bool writeHeader(const std::string &header)
      {
        return m_file.writeHeader(header);
      }

    private:
      BinaryOutputFile m_file; //!< The output file.
      unsigned int m_numBits; //!< The number of bits in the fingerprint.
      unsigned int m_numFingerprints; //!< The number of fingerprints that will be written.
      unsigned int m_current; //!< The current fingerprint being written.
      Word *m_data; //!< Pointer to all fingerprint bits.
  };

  /**
   * @brief Class for accessing row-major order fingerprints file from memory.
   */
  class InMemoryRowMajorFingerprintStorage
  {
    public:
      /**
       * @brief Constructor.
       */
      InMemoryRowMajorFingerprintStorage() : m_fingerprints(0), m_numBits(0),
          m_numFingerprints(0), m_init(false)
      {
      }

      /**
       * @brief Destructor.
       */
      ~InMemoryRowMajorFingerprintStorage()
      {
        if (m_fingerprints)
          delete [] m_fingerprints;
      }

      /**
       * @brief Get the JSON file header.
       *
       * @return The JSON file header.
       */
      std::string header() const
      {
        return m_json;
      }

      /**
       * @brief Get the number of fingerprint bits.
       *
       * @return The number of fingerprint bits.
       */
      unsigned int numBits() const
      {
        return m_numBits;
      }

      /**
       * @brief Get the number of fingerprints.
       *
       * @return The number of fingerprints.
       */
      unsigned int numFingerprints() const
      {
        return m_numFingerprints;
      }

      /**
       * @brief Get the fingerprint with the specified index.
       *
       * @param index The index of the molecule for which to get the fingerprints.
       *
       * @return The fingerprint with the specified index.
       */
      Word* fingerprint(unsigned int index) const
      {
        if (!m_init)
          return 0;
        return m_fingerprints + bitvec_num_words_for_bits(m_numBits) * index;
      }

      /**
       * @brief Load the fingerprints from a file in memory.
       *
       * An exception is thrown when an error occurs.
       *
       * @param filename The filename.
       */
      void load(const std::string &filename)
      {
        TIMER("InMemoryRowMajorFingerprintStorage::load():");

        // open the file
        BinaryInputFile file(filename);
        if (!file)
          throw std::runtime_error(make_string("Could not open fingerprint file \"", filename, "\""));

        // parse the JSON header
        m_json = file.header();
        Json::Reader reader;
        Json::Value data;
        if (!reader.parse(m_json, data))
          throw std::runtime_error(reader.getFormattedErrorMessages());

        // make sure the required attributes are present
        if (!data.isMember("filetype") || data["filetype"].asString() != "fingerprints")
          throw std::runtime_error(make_string("JSON header for file ", filename, " does not contain 'filetype' attribute or is not 'fingerprints'"));
        if (!data.isMember("order"))
          throw std::runtime_error(make_string("JSON header for file ", filename, " does not contain 'order' attribute"));
        if (!data.isMember("num_bits"))
          throw std::runtime_error(make_string("JSON header for file ", filename, " does not contain 'num_bits' attribute"));
        if (!data.isMember("num_fingerprints"))
          throw std::runtime_error(make_string("JSON header for file ", filename, " does not contain 'num_fingerprints' attribute"));

        // ensure this is a row-major order file
        if (data["order"].asString() != "row-major")
          throw std::runtime_error(make_string(filename, " is not a row-major order fingerprint file"));

        // get attributes from header
        m_numBits = data["num_bits"].asUInt();
        m_numFingerprints = data["num_fingerprints"].asUInt();

        // allocate memory
        m_fingerprints = new Word[bitvec_num_words_for_bits(m_numBits) * m_numFingerprints];
        file.read(m_fingerprints, bitvec_num_words_for_bits(m_numBits) * m_numFingerprints * sizeof(Word));

        m_init = true;
      }

    private:
      std::string m_json; //!< JSON header
      Word *m_fingerprints; //!< The fingerprint data.
      unsigned int m_numBits; //!< The number of fingerprint bits.
      unsigned int m_numFingerprints; //!< The number of fingerprints.
      bool m_init; //!< Flag to check for initialization.
  };

  /**
   * @brief Class for accessing column-major order fingerprints file from memory.
   */
  class InMemoryColumnMajorFingerprintStorage
  {
    public:
      /**
       * @brief Constructor.
       */
      InMemoryColumnMajorFingerprintStorage() : m_fingerprints(0), m_numBits(0),
          m_numFingerprints(0), m_init(false)
      {
      }

      /**
       * @brief Destructor.
       */
      ~InMemoryColumnMajorFingerprintStorage()
      {
        delete [] m_fingerprints;
      }

      /**
       * @brief Get the JSON file header.
       *
       * @return The JSON file header.
       */
      std::string header() const
      {
        return m_json;
      }

      /**
       * @brief Get the number of fingerprint bits.
       *
       * @return The number of fingerprint bits.
       */
      unsigned int numBits() const
      {
        return m_numBits;
      }

      /**
       * @brief Get the number of fingerprints.
       *
       * @return The number of fingerprints.
       */
      unsigned int numFingerprints() const
      {
        return m_numFingerprints;
      }

      /**
       * @brief Get the column with the specified index.
       *
       * @param index The index of the fingerprint bit.
       *
       * @return The column with the specified index.
       */
      Word* bit(unsigned int index) const
      {
        if (!m_init)
          return 0;
        return m_fingerprints + bitvec_num_words_for_bits(m_numFingerprints) * index;
      }

      /**
       * @brief Load the fingerprints from a file in memory.
       *
       * An exception is thrown when an error occurs.
       *
       * @param filename The filename.
       */
      void load(const std::string &filename)
      {
        TIMER("InMemoryColumnMajorFingerprintStorage::load():");

        // open the file
        BinaryInputFile file(filename);
        if (!file)
          throw std::runtime_error(make_string("Could not open fingerprint file \"", filename, "\""));

        // parse the JSON header
        m_json = file.header();
        Json::Reader reader;
        Json::Value data;
        if (!reader.parse(m_json, data))
          throw std::runtime_error(reader.getFormattedErrorMessages());

        // make sure the required attributes are present
        if (!data.isMember("filetype") || data["filetype"].asString() != "fingerprints")
          throw std::runtime_error(make_string("JSON header for file ", filename, " does not contain 'filetype' attribute or is not 'fingerprints'"));
        if (!data.isMember("order"))
          throw std::runtime_error(make_string("JSON header for file ", filename, " does not contain 'order' attribute"));
        if (!data.isMember("num_bits"))
          throw std::runtime_error(make_string("JSON header for file ", filename, " does not contain 'num_bits' attribute"));
        if (!data.isMember("num_fingerprints"))
          throw std::runtime_error(make_string("JSON header for file ", filename, " does not contain 'num_fingerprints' attribute"));

        // ensure this is a row-major order file
        if (data["order"].asString() != "column-major")
          throw std::runtime_error(make_string(filename, " is not a column-major order fingerprint file"));

        // get attributes from header
        m_numBits = data["num_bits"].asUInt();
        m_numFingerprints = data["num_fingerprints"].asUInt();

        // allocate memory
        m_fingerprints = new Word[bitvec_num_words_for_bits(m_numFingerprints) * m_numBits];
        file.read(m_fingerprints, bitvec_num_words_for_bits(m_numFingerprints) * m_numBits * sizeof(Word));

        m_init = true;
      }

    private:
      std::string m_json; //!< JSON header
      Word *m_fingerprints; //!< The fingerprint data.
      unsigned int m_numBits; //!< The number of fingerprint bits.
      unsigned int m_numFingerprints; //!< The number of fingerprints.
      bool m_init; //!< Flag to check for initialization.
  };

}

#endif
