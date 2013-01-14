#ifndef HELIUM_BITVEC_H
#define HELIUM_BITVEC_H

#include <iostream>
#include <fstream>

namespace Helium {

  typedef unsigned long Word;

  const int BitsPerWord = sizeof(Word) * 8;

  inline unsigned int num_words_for_bits(unsigned int numBits)
  {
    return (numBits + numBits % BitsPerWord) / BitsPerWord;
  }

  inline void zero(Word *bitvec, int numWords)
  {
    for (int i = 0; i < numWords; ++i)
      bitvec[i] = 0;
  }

  inline Word* copy(const Word *bitvec, int numWords)
  {
    Word *result = new Word[numWords];
    std::copy(bitvec, bitvec + numWords, result);
    return result;
  }

  inline bool is_subset_superset(const Word *bitvec1, const Word *bitvec2, int numWords)
  {
    for (int i = 0; i < numWords; ++i)
      if (bitvec1[i] & ~bitvec2[i])
        return false;
    return true; 
  }

  inline int bit_count(Word word)
  {
#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
    return __builtin_popcountl(word);
#else
    int count = 0;
    Word bit = 1;
    for (int i = 0; i < BitsPerWord; ++i) {
      if (bit & word)
        ++count;
      bit <<= 1;
    }
    return count;
#endif
  }

  inline int bit_count(const Word *bitvec, int numWords)
  {
#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
    // use popcount
    int count = 0;
    for (int i = 0; i < numWords; ++i)
      count += __builtin_popcountl(*(bitvec + i));
    return count;
#else
    int count = 0;
    for (int i = 0; i < numWords; ++i) {
      Word bit = 1;
      for (int j = 0; j < BitsPerWord; ++j) {
        if (bit & word)
          ++count;
        bit <<= 1;
      }
    }
    return count;
#endif
  }

  inline void clear(Word *bitvec, int numWords)
  {
    for (int i = 0; i < numWords; ++i)
      *(bitvec + i) = 0;
  }

  inline bool get(int index, const Word *bitvec, int numWord = 0) // FIXME
  {
    int word = index / (sizeof(Word) * 8);
    int offset = index % (sizeof(Word) * 8);
    Word bit = static_cast<Word>(1) << offset;
    return *(bitvec + word) & bit;
  }

  inline void set(int index, Word *bitvec, int numWord = 0) // FIXME
  {
    int word = index / (sizeof(Word) * 8);
    int offset = index % (sizeof(Word) * 8);
    Word bit = static_cast<Word>(1) << offset;
    *(bitvec + word) |= bit;
  }

  inline void reset(int index, Word *bitvec, int numWord = 0) // FIXME
  {
    int word = index / (sizeof(Word) * 8);
    int offset = index % (sizeof(Word) * 8);
    Word bit = static_cast<Word>(1) << offset;
    *(bitvec + word) &= ~bit;
  }

  inline void print(Word word)
  {
    Word bit = 1;
    for (int j = 0; j < BitsPerWord; ++j) {
      //if (j != 0 && (j % 16) == 0)
      //  std::cout << " ";
      if (bit & word)
        std::cout << "1";
      else
        std::cout << "0";
      bit <<= 1;
    }
    std::cout << std::endl;
  }

  inline void print(const Word *bitvec, int numWords, bool spaces = true)
  {
    for (int i = 0; i < numWords; ++i) {
      Word bit = 1;
      for (int j = 0; j < BitsPerWord; ++j) {
        if (spaces && j != 0 && (j % 16) == 0)
          std::cout << " ";
        if (bit & *(bitvec + i))
          std::cout << "1";
        else
          std::cout << "0";
        bit <<= 1;
      }
      if (spaces && i + 1 < numWords)
        std::cout << " ";
    }
    std::cout << std::endl;
  }

  inline void write_size(std::ofstream &ofs, unsigned long size)
  {
    ofs.seekp(0);
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(unsigned long));
  }

  inline unsigned long read_size(std::ifstream &ifs)
  {
    unsigned long size;
    ifs.seekg(0);
    ifs.read(reinterpret_cast<char*>(&size), sizeof(unsigned long));
    return size;
  }

  inline void write(std::ofstream &ofs, const Word *bitvec, int numWords)
  {
    ofs.write(reinterpret_cast<const char*>(bitvec), sizeof(Word) * numWords);
  }

  inline void read(std::ifstream &ifs, Word *bitvec, int numWords)
  {
    ifs.read(reinterpret_cast<char*>(bitvec), sizeof(Word) * numWords);
  }


}

#endif
