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
#include "file.h"

#include <cassert>

namespace Helium {

  bool BinaryInputFile::open(const std::string &filename)
  {
    // reset JSON header
    m_json.clear();

    // try to open the file
    m_ifs.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!m_ifs)
      throw std::runtime_error(make_string("Could not open file \"", filename, "\""));

    // read the magic number
    unsigned int magic;
    if (!m_ifs.read(reinterpret_cast<char*>(&magic), sizeof(unsigned int)) || magic != 0x48650001) {
      m_ifs.close();
      throw std::runtime_error(make_string("Invalid magic number in file \"", filename, "\""));
    }

    // read length of JSON header
    unsigned int jsonSize;
    if (!m_ifs.read(reinterpret_cast<char*>(&jsonSize), sizeof(unsigned int))) {
      m_ifs.close();
      throw std::runtime_error(make_string("Invalid JSON header in file \"", filename, "\""));
    }

    // read json header
    char *buffer = new char[jsonSize];
    if (!m_ifs.read(buffer, jsonSize)) {
      m_ifs.close();
      throw std::runtime_error(make_string("Invalid JSON header in file \"", filename, "\""));
    }

    // copy header
    m_json = buffer;
    delete [] buffer;

    // set offset
    m_offset = m_ifs.tellg();

    return (bool)m_ifs;
  }

#define HEADER_SIZE 128000

  bool BinaryOutputFile::open(const std::string &filename)
  {
    // try to open the file
    m_ofs.open(filename.c_str(), std::ios_base::out | std::ios_base::binary);
    if (!m_ofs)
      throw std::runtime_error(make_string("Could not open file \"", filename, "\""));

    // write the magic number
    unsigned int magic = 0x48650001;
    m_ofs.write(reinterpret_cast<const char*>(&magic), sizeof(unsigned int));

    // set position to start of binary data
    m_ofs.seekp(2 * sizeof(unsigned int) + HEADER_SIZE);
    m_offset = m_ofs.tellp();

    return (bool)m_ofs;
  }

  bool BinaryOutputFile::writeHeader(const std::string &header)
  {
    assert(header.size() < HEADER_SIZE - 1);

    // store current file stream position
    std::ios_base::streampos pos = m_ofs.tellp();

    // set position to start of header
    m_ofs.seekp(sizeof(unsigned int));

    // write header size
    unsigned int jsonSize = HEADER_SIZE;
    m_ofs.write(reinterpret_cast<const char*>(&jsonSize), sizeof(unsigned int));

    // write header
    m_ofs.write(header.c_str(), header.size());

    // write padding
    const char padding = '\0';
    for (unsigned int i = 0; i < HEADER_SIZE - header.size(); ++i)
      m_ofs.write(&padding, 1);

    // restore file stream position
    m_ofs.seekp(pos);

    return (bool)m_ofs;
  }

}
