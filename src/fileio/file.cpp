#include "file.h"

#include <cassert>

namespace Helium {

  bool BinaryInputFile::open(const std::string &filename)
  {
    // reset error flag
    m_error = NoError;
    // reset JSON header
    m_json.clear();

    // try to open the file
    m_ifs.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
    if (!m_ifs) {
      m_error = CouldNotOpen;
      return false;
    }

    // read the magic number
    unsigned int magic;
    if (!m_ifs.read(reinterpret_cast<char*>(&magic), sizeof(unsigned int)) || magic != 0x48650001) {
      m_error = InvalidMagic;
      m_ifs.close();
      return false;
    }

    // read length of JSON header
    unsigned int jsonSize;
    if (!m_ifs.read(reinterpret_cast<char*>(&jsonSize), sizeof(unsigned int))) {
      m_error = InvalidHeader;
      m_ifs.close();
      return false;
    }

    // read json header
    char *buffer = new char[jsonSize];
    if (!m_ifs.read(buffer, jsonSize)) {
      m_error = InvalidHeader;
      m_ifs.close();
      return false;
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
    // reset error flag
    m_error = NoError;

    // try to open the file
    m_ofs.open(filename.c_str(), std::ios_base::out | std::ios_base::binary);
    if (!m_ofs) {
      m_error = CouldNotOpen;
      return false;
    }

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
