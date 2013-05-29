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
#ifndef HELIUM_FILEIO_FPS_H
#define HELIUM_FILEIO_FPS_H

#include <Helium/bitvec.h>

namespace Helium {

  /**
   * @brief Class for reading FPS fingerprint files.
   *
   * FPS is a text format for storing fingerprints. This class can read FPS1
   * files and will store all read fingerprints and store them memory for this
   * class' lifetime. This class expects to find a num_bits field in the FPS
   * header (e.g. "#num_bits=1024").
   *
   * The FPS specification can be found here: https://code.google.com/p/chem-fingerprints/wiki/FPS
   */
  class FpsFile
  {
    public:
      /**
       * @brief Constructor.
       */
      FpsFile() : m_numFingerprints(0), m_numBits(0)
      {
      }

      /**
       * @brief Load the fingerprints form the specified filename.
       *
       * @param filename The FPS file.
       *
       * @return True if the file was parsed successfully.
       */
      bool load(const std::string &filename)
      {
        m_fingerprints.clear();

        std::ifstream ifs(filename.c_str());
        if (!ifs)
          return false;

        m_numFingerprints = 0;
        m_numBits = 0;
        unsigned int numWords = 0;

        m_type = "FPS";

        std::string line;
        // read fps header
        while (ifs && !ifs.eof() && ifs.peek() == '#') {
          std::getline(ifs, line);
          if (line.substr(0, 10) == "#num_bits=") {
            std::stringstream ss(line.substr(10));
            ss >> m_numBits;
            numWords = bitvec_num_words_for_bits(m_numBits);
          } else if (line.substr(0, 6) == "#type=") {
            m_type = line.substr(6);
          }
        }

        if (!m_numBits)
          return false;

        // allocate bit vector
        Word *fingerprint = new Word[numWords];

        // process fingerprints
        while (std::getline(ifs, line)) {
          if (line.empty())
            continue;

          std::size_t nibbles = line.find(' ');
          if (nibbles == std::string::npos)
            nibbles = line.find('\t');
          if (nibbles == std::string::npos)
            nibbles = line.size();

          // convert hex to fingerprint
          hex_to_bitvec(line.substr(0, nibbles), fingerprint, numWords);

          std::copy(fingerprint, fingerprint + numWords, std::back_inserter(m_fingerprints));
          ++m_numFingerprints;
        }

        return true;
      }

      /**
       * @brief Get the type of fingerprint.
       *
       * The type is extracted from the FPS file header (e.g. "#type=OpenBabel-FP2/1").
       * If there is no type field in the header, this function will return "FPS".
       *
       * @return The type of fingerprint.
       */
      const std::string& type() const
      {
        return m_type;
      }

      /**
       * @brief Get the number of bits.
       *
       * The number of bits are extracted from the FPS file header
       * (e.g. "#num_bits=1024").
       *
       * @return The number of bits in the fingerprints.
       */
      unsigned int numBits() const
      {
        return m_numBits;
      }

      /**
       * @brief Get the number of fingerprints loaded.
       *
       * @return The number of fingerprints loaded.
       */
      unsigned int numFingerprints() const
      {
        return m_numFingerprints;
      }

      /**
       * @brief Get the fingerprint with the specified index.
       *
       * @param index The index of the fingerprint to get.
       *
       * @return The fingerprint with the specified index.
       */
      Word* fingerprint(unsigned int index) const
      {
        if (index >= m_numFingerprints)
          return 0;
        return const_cast<Word*>(m_fingerprints.data()) + bitvec_num_words_for_bits(m_numBits) * index;
      }

    private:
      std::vector<Word> m_fingerprints; //!< The fingerprints.
      std::string m_type; //!< The fingerprint type.
      unsigned int m_numFingerprints; //!< The number of fingerprints.
      unsigned int m_numBits; //!< The number of bytes.
  };

}

#endif
