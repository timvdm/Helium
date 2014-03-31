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
#ifndef HELIUM_FILEIO_FILE_H
#define HELIUM_FILEIO_FILE_H

#include <Helium/util.h>

#include <fstream>
#include <stdexcept>
#include <string>

namespace Helium {

  /**
   * @page binary_files Helium Binary Files
   *
   * The exact format is described below but for using these file
   * formats it is enough to know that all files contain a JSON header
   * followed by the binary data itself.
   *
   * @section binary_files_format The Binary File Format
   *
   * The first 4 bytes are a magic number (0x48650001). The next 4 bytes contain
   * the length (in bytes) of the JSON header including zero termination
   * character. The JSON header starts at position 8 and ends at 8 plus the
   * number of bytes contained in the second 4 byte field. The first 128000
   * bytes are always reserved for the header in case it needs more space in
   * the future. What follows after the JSON header is binary data and the
   * format depends on the specific file format.
   */

  /**
   * @brief Helper class for all Helium binary input file formats.
   *
   * The HeliumInputFile class is a helper class used for all Helium binary
   * input file formats.
   */
  class BinaryInputFile
  {
    public:
      /**
       * @brief Constructor.
       */
      BinaryInputFile()
      {
      }

      /**
       * @brief Constructor.
       *
       * @param filename The filename.
       */
      BinaryInputFile(const std::string &filename)
      {
        open(filename);
      }

      /**
       * Open the file with the specified filename. This method also checks the
       * file's magic number and reads the JSON header. If this method returns
       * false the error() method can be used to get the type of error that
       * occured.
       *
       * @param filename The filename for the file to open.
       *
       * @return True if the file was opened successfully, the magic number is
       *         valid and the JSON header was read successfully.
       */
      bool open(const std::string &filename);

      /**
       * Close the file.
       */
      void close()
      {
        m_json.clear();
        m_ifs.close();
      }

      /**
       * Convert to bool, can be used to check if stream is valid.
       *
       * @return True if the file stream is valid.
       */
      operator bool() const
      {
        return m_ifs.is_open();
      }

      /**
       * Get the current file stream position. The binary data starts at
       * position 0.
       *
       * @return The file stream position.
       */
      std::ios_base::streampos tell()
      {
        return m_ifs.tellg() - m_offset;
      }

      /**
       * Set the file stream position.
       *
       * @param pos The new file stream position.
       */
      bool seek(std::ios_base::streampos pos)
      {
        return (bool)m_ifs.seekg(pos + m_offset);
      }

      /**
       * Get the JSON header.
       *
       * @return The JSON header.
       */
      const std::string& header() const
      {
        return m_json;
      }

      /**
       * Read data from the file.
       *
       * @param data The data to read from the file.
       * @param size The size of the data to read in bytes.
       *
       * @return True if the data was read successfully.
       */
      bool read(char *data, std::size_t size)
      {
        return (bool)m_ifs.read(data, size);
      }

      /**
       * @brief Read data from the file.
       *
       * The @p data pointer must be allocated to store the read data.
       *
       * @param data Pointer to memory to store the read data.
       * @param size The number of bytes to read.
       */
      template<typename T>
      bool read(T *data, std::size_t size)
      {
        return (bool)m_ifs.read(reinterpret_cast<char*>(data), size);
      }

      /**
       * Get the std::ifstream to manipulate it directly.
       *
       * @return The std::ifstream.
       */
      std::ifstream& stream()
      {
        return m_ifs;
      }

    private:
      std::ifstream m_ifs; //!< File handle
      std::string m_json; //!< JSON header
      std::ios_base::streampos m_offset; //!< Offset where binary data starts
  };

  /**
   * @brief Helper class for all Helium binary output file formats.
   *
   * The HeliumOutputFile class is a helper class used for all Helium binary
   * output file formats.
   */
  class BinaryOutputFile
  {
    public:
      /**
       * @brief Constructor.
       */
      BinaryOutputFile()
      {
      }

      /**
       * @brief Constructor.
       *
       * @param filename The filename.
       */
      BinaryOutputFile(const std::string &filename)
      {
        open(filename);
      }

      /**
       * Open the file with the specified filename.
       *
       * @param filename The filename for the file to open.
       *
       * @return True if the file was opened successfully.
       */
      bool open(const std::string &filename);

      /**
       * Close the file.
       */
      void close()
      {
        m_ofs.close();
      }

      /**
       * Convert to bool, can be used to check if stream is valid.
       *
       * @return True if the file stream is valid.
       */
      operator bool() const
      {
        return (bool)m_ofs;
      }

      /**
       * Get the current file stream position. The binary data starts at
       * position 0.
       *
       * @return The file stream position.
       */
      std::ios_base::streampos tell()
      {
        return m_ofs.tellp() - m_offset;
      }

      /**
       * Set the file stream position.
       *
       * @param pos The new file stream position.
       */
      bool seek(std::ios_base::streampos pos)
      {
        return (bool)m_ofs.seekp(pos + m_offset);
      }

      /**
       * Write data to the file.
       *
       * @param data The data to write to the file.
       * @param size The size of the data to write in bytes.
       *
       * @return True if the data was written successfully.
       */
      bool write(const char *data, std::size_t size)
      {
        return (bool)m_ofs.write(data, size);
      }

      /**
       * @overload
       */
      template<typename T>
      bool write(const T *data, std::size_t size)
      {
        return (bool)m_ofs.write(reinterpret_cast<const char*>(data), size);
      }

      /**
       * Write the JSON header to the file.
       *
       * @param header The JSON header.
       *
       * @return true if the header was written successfully.
       */
      bool writeHeader(const std::string &header);

      /**
       * Get the std::ofstream to manipulate it directly.
       */
      std::ofstream& stream()
      {
        return m_ofs;
      }

    private:
      std::ofstream m_ofs; //!< File handle
      std::ios_base::streampos m_offset; //!< Offset where binary data starts
  };

}

#endif
