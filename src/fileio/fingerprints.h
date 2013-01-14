#ifndef HELIUM_FILEIO_FINGERPRINTS_H
#define HELIUM_FILEIO_FINGERPRINTS_H

#include "bitvec.h"
#include "enumeratepaths.h"
#include "enumeratesubgraphs.h"
#include "substructure.h"
#include "extendedconnectivities.h"
#include "canonical.h"

#include <stdexcept>

#include <boost/functional/hash.hpp>

namespace Helium {

  class FingerprintFile
  {
    public:
      FingerprintFile(const std::string &filename) : m_ifs(filename.c_str()), m_current(-1)
      {
        if (m_ifs)
          read32(m_ifs, m_numFingerprints);    
        m_numWords = 16;
      }

      unsigned int num_fingerprints() const
      {
        return m_numFingerprints;
      }

      unsigned int current() const
      {
        return m_current;
      }

      bool read_fingerprint(Word *fingerprint)
      {
        if (!m_ifs)
          return false;
        ++m_current;
        if (m_current == m_numFingerprints)
          return false;
        for (int i = 0; i < m_numWords; ++i)
          read64(m_ifs, fingerprint[i]);
        return m_ifs;
      }

    private:
      std::ifstream m_ifs;
      unsigned int m_numFingerprints;
      unsigned int m_current;
      int m_numWords;
  };

  struct InvertedFingerprintFileHeader
  {
    InvertedFingerprintFileHeader() : magic_number(get_magic_number())
    {
    }

    static unsigned int get_magic_number()
    {
      return 0x48650001;
    }

    unsigned int magic_number;
    unsigned int bits_per_word;
    unsigned int bits_per_fingerprint;
    unsigned int words_per_fingerprint;
    unsigned int words_per_fpbit;
    unsigned int num_fingerprints;
  };

  class InvertedFingerprintOutputFile
  {
    public:
      InvertedFingerprintOutputFile(unsigned int bitsPerFingerprint, unsigned int numFingerprints, const std::string &filename) 
          : m_ofs(filename.c_str(), std::ios_base::out | std::ios_base::binary), m_current(0)
      {
        // setup m_header
        m_header.bits_per_word = sizeof(Word) * 8; // e.g. 64
        m_header.bits_per_fingerprint = bitsPerFingerprint; // e.g. 1024
        m_header.words_per_fingerprint = m_header.bits_per_fingerprint / m_header.bits_per_word; // e.g. 16
        m_header.words_per_fpbit = (numFingerprints + numFingerprints % m_header.bits_per_word) / m_header.bits_per_word;
        m_header.num_fingerprints = numFingerprints;

        // write m_header
        if (!m_ofs)
          throw std::runtime_error(make_string("Could not open ", filename, " for writing."));
        
        // write the header
        m_ofs.write(reinterpret_cast<const char*>(&m_header), sizeof(InvertedFingerprintFileHeader));
        // allocate data
        m_data = new Word[m_header.words_per_fpbit * m_header.bits_per_fingerprint];
      }

      ~InvertedFingerprintOutputFile()
      {
        // write the data
        m_ofs.write(reinterpret_cast<const char*>(m_data), m_header.words_per_fpbit * m_header.bits_per_fingerprint * sizeof(Word));
        delete [] m_data;
      }

      void write(Word *fingerprint)
      {
        // check each bit in the fingerprint
        for (int i = 0; i < m_header.bits_per_fingerprint; ++i) {
          // skip this bit if it is not set
          if (!get(i, fingerprint))
            continue;

          // set the correct bit
          set(i * m_header.words_per_fpbit * 64 + m_current, m_data);
        }
        ++m_current;
      }

    private:
      InvertedFingerprintFileHeader m_header;
      std::ofstream m_ofs;
      unsigned int m_current;
      Word *m_data;
  };

  class InvertedFingerprintFile
  {
    public:
      InvertedFingerprintFile(const std::string &filename) : m_ifs(filename.c_str(), std::ios_base::in | std::ios_base::binary)
      {
        if (!m_ifs)
          throw std::runtime_error(make_string("Could not open ", filename, " for reading."));

        // read the header
        m_ifs.read(reinterpret_cast<char*>(&m_header), sizeof(InvertedFingerprintFileHeader));

        // check magic number
        if (m_header.magic_number != m_header.get_magic_number())
          throw std::runtime_error(make_string(filename, " is not an inverted fingerprint file."));

        // allocate memory
        m_data = new Word[m_header.words_per_fpbit];
      }

      ~InvertedFingerprintFile()
      {
        delete [] m_data;
      }

      unsigned int num_fingerprints() const
      {
        return m_header.num_fingerprints;
      }

      Word* allocate_result() const
      {
        return new Word[m_header.words_per_fpbit];
      }

      void search(Word *fingerprint, Word *result)
      {
        bool first = true;
        for (int i = 0; i < m_header.bits_per_fingerprint; ++i) { // foreach bit
          // skip this bit if it is not set
          if (!get(i, fingerprint))
            continue;
        

          // compute offset for this fingerprint bit
          unsigned int offset = i * m_header.words_per_fpbit;
          
          m_ifs.seekg(sizeof(InvertedFingerprintFileHeader) + offset * sizeof(Word));
          m_ifs.read(reinterpret_cast<char*>(m_data), m_header.words_per_fpbit * sizeof(Word));

          if (first) {
            // if this is the first bit, just set result
            for (unsigned int j = 0; j < m_header.words_per_fpbit; ++j)
              result[j] = m_data[j];
            first = false;
          } else {
            // do set intersection
            for (unsigned int j = 0; j < m_header.words_per_fpbit; ++j)
              result[j] &= m_data[j];
          }
        }
      }

    private:
      InvertedFingerprintFileHeader m_header;
      std::ifstream m_ifs;
      Word *m_data;
  };

  class InvertedFingerprintFileCached
  {
    public:
      InvertedFingerprintFileCached(const std::string &filename)
      {
        // open the file
        std::ifstream ifs(filename.c_str(), std::ios_base::in | std::ios_base::binary);
        if (!ifs)
          throw std::runtime_error(make_string("Could not open ", filename, " for reading."));

        // read the header
        ifs.read(reinterpret_cast<char*>(&m_header), sizeof(InvertedFingerprintFileHeader));

        // check magic number
        if (m_header.magic_number != m_header.get_magic_number())
          throw std::runtime_error(make_string(filename, " is not an inverted fingerprint file."));

        // allocate memory
        m_data = new Word[m_header.words_per_fpbit * m_header.bits_per_fingerprint];
        // read the data
        ifs.read(reinterpret_cast<char*>(m_data), m_header.words_per_fpbit * m_header.bits_per_fingerprint * sizeof(Word));
      }

      ~InvertedFingerprintFileCached()
      {
        delete [] m_data;
      }

      unsigned int num_fingerprints() const
      {
        return m_header.num_fingerprints;
      }

      Word* allocate_fingerprint() const
      {
        return new Word[m_header.words_per_fingerprint];
      }
      
      Word* allocate_result() const
      {
        return new Word[m_header.words_per_fpbit];
      }

      void search(Word *fingerprint, Word *result)
      {
        bool first = true;
        for (int i = 0; i < m_header.bits_per_fingerprint; ++i) { // foreach bit
          // skip this bit if it is not set
          if (!get(i, fingerprint))
            continue;

          // compute offset for this fingerprint bit
          unsigned int offset = i * m_header.words_per_fpbit;

          if (first) {
            // if this is the first bit, just set result
            for (unsigned int j = 0; j < m_header.words_per_fpbit; ++j)
              result[j] = m_data[offset + j];
            first = false;
          } else {
            // do set intersection
            for (unsigned int j = 0; j < m_header.words_per_fpbit; ++j)
              result[j] &= m_data[offset + j];
          }
        }
      }

    private:
      InvertedFingerprintFileHeader m_header;
      Word *m_data;
  };


  
}

#endif
