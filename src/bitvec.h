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
#ifndef HELIUM_BITVEC_H
#define HELIUM_BITVEC_H

#include <iostream>
#include <fstream>

namespace Helium {

  /**
   * @file bitvec.h
   * @brief Functions for working with bit vectors.
   */

  /**
   * Type used for a bit vector word.
   */
  typedef unsigned long Word;

  /**
   * The number of bits per bit vector word.
   */
  const int BitsPerWord = sizeof(Word) * 8;

  /**
   * Get the number of words needed to store @p numBits words.
   *
   * @param numBits The number of bits.
   *
   * @return The number of words needed to store @p numBits bits.
   */
  inline unsigned int bitvec_num_words_for_bits(unsigned int numBits)
  {
    return (numBits + numBits % BitsPerWord) / BitsPerWord;
  }

  /**
   * Set all bits to 0.
   *
   * @param bitvec The bit vector to zero.
   * @param numWords The number of words for @p bitvec.
   */
  inline void bitvec_zero(Word *bitvec, int numWords)
  {
    for (int i = 0; i < numWords; ++i)
      bitvec[i] = 0;
  }

  /**
   * Copy a bit vector. The copied vector's memory should be freed using delete [].
   *
   * @param bitvec The bit vector to copy.
   * @param numWords The number of words for @p bitvec.
   *
   * @return A pointer to the copied bit vector.
   */
  inline Word* bitvec_copy(const Word *bitvec, int numWords)
  {
    Word *result = new Word[numWords];
    std::copy(bitvec, bitvec + numWords, result);
    return result;
  }

  /**
   * Get the bit value for the specified bit index.
   *
   * @param index The bit's index. The index should be less than the number of
   * bits in the @p bitvec.
   * @param bitvec The bit vector.
   *
   * @return The valuefor the bit at the specified @p index.
   */
  inline bool bitvec_get(int index, const Word *bitvec)
  {
    int word = index / (sizeof(Word) * 8);
    int offset = index % (sizeof(Word) * 8);
    Word bit = static_cast<Word>(1) << offset;
    return *(bitvec + word) & bit;
  }

  /**
   * Set the bit at the specified bit index (to 1).
   *
   * @param index The bit's index. The index should be less than the number of
   * bits in the @p bitvec.
   * @param bitvec The bit vector.
   */
  inline void bitvec_set(int index, Word *bitvec)
  {
    int word = index / (sizeof(Word) * 8);
    int offset = index % (sizeof(Word) * 8);
    Word bit = static_cast<Word>(1) << offset;
    *(bitvec + word) |= bit;
  }

  /**
   * Reset the bit at the specified bit index (to 0).
   *
   * @param index The bit's index. The index should be less than the number of
   * bits in the @p bitvec.
   * @param bitvec The bit vector.
   */
  inline void bitvec_reset(int index, Word *bitvec)
  {
    int word = index / (sizeof(Word) * 8);
    int offset = index % (sizeof(Word) * 8);
    Word bit = static_cast<Word>(1) << offset;
    *(bitvec + word) &= ~bit;
  }

  /**
   * Check if @p bitvec1 is a subset of @p bitvec2 (i.e. all bits set in
   * @p bitvec1 are also set in @p bitvec2). Both bitvectors should have the
   * same number of bits.
   *
   * @param bitvec1 The subset bit vector.
   * @param bitvec2 The superset bit vector.
   * @param numWords The number of words for @p bitvec1 and @p bitvec2.
   *
   * @return True if @p bitvec1 is a subset of @p bitvec2.
   */
  inline bool bitvec_is_subset_superset(const Word *bitvec1, const Word *bitvec2, int numWords)
  {
    for (int i = 0; i < numWords; ++i)
      if (bitvec1[i] & ~bitvec2[i])
        return false;
    return true;
  }

  /**
   * Get the bit count (i.e. number of bits set to 1 or the population count)
   * for a single bit vector word. When compiling with g++, the builtin function
   * wrapper for the POPCNT instruction is used.
   *
   * @param word The bit vector word.
   *
   * @return The bit count.
   */
  inline int bitvec_count(Word word)
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

  /**
   * Get the bit count (i.e. number of bits set to 1 or the population count)
   * for a single bit vector word. When compiling with g++, the builtin function
   * wrapper for the POPCNT instruction is used.
   *
   * @param bitvec The bit vector.
   * @param numWords The number of words for @p bitvec.
   *
   * @return The bit count.
   */
  inline int bitvec_count(const Word *bitvec, int numWords)
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


  /**
   * Get the bit count (i.e. number of bits set to 1 or the population count)
   * for a range of a bit vector.
   *
   * @param bitvec The bit vector.
   * @param begin The first bit index.
   * @param end The one-past-the-end bit index.
   *
   * @return The bit count for the range [begin,end).
   */
  inline int bitvec_count(const Word *bitvec, int begin, int end)
  {
    const int bits_per_word = 8 * sizeof(Word);
    const int lft = begin / bits_per_word;
    const int rgt = end / bits_per_word;

    int count = 0;
    // count bits in partial word on the left
    for (int i = lft * bits_per_word + (begin % bits_per_word); i < (lft + 1) * bits_per_word; ++i)
      count += bitvec_get(i, bitvec);
    // count bits in full words in the middle
    count += bitvec_count(bitvec + lft + 1, rgt - lft - 1);
    // count bits in partial word on the right
    for (int i = rgt * bits_per_word; i < rgt * bits_per_word + (end % bits_per_word); ++i)
      count += bitvec_get(i, bitvec);

    return count;
  }


  /**
   * Compute the Tanimoto coefficient of difference between two bit vectors.
   * This function is slower than the Tanimoto function below since both the
   * union and intersection bit counts have to be computed.
   *
   * \f[
   *   T_{\mathrm{sim}} = \frac{| A \wedge B |}{| A \vee B |}
   * \f]
   */
  inline double bitvec_tanimoto(const Word *bitvec1, const Word *bitvec2, int numWords)
  {
    int andCount = 0;
    int orCount = 0;

    for (int i = 0; i < numWords; ++i) {
      Word andbv = bitvec1[i] & bitvec2[i];
      Word orbv = bitvec1[i] | bitvec2[i];
#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
      //andCount += __builtin_popcountl(andfp);
      //orCount += __builtin_popcountl(orfp);
      andCount += __builtin_popcountll(andbv);
      orCount += __builtin_popcountll(orbv);
#else
      for (; andbv; andbv = andbv << 1)
        if(andbv < 0)
          ++andbits;
      for (; orbv; orbv = orbv << 1)
        if(orbv < 0)
          ++orbits;
#endif
    }

    return static_cast<double>(andCount) / orCount;
  }

  /**
   * Compute the Tanimoto coefficient of difference between two bit vectors.
   * When the bit counts of both bit vectors are known, the inclusion-exclusion
   * principle can be used to compute the Tanimoto coefficient faster.
   *
   * \f[
   *   T_{\mathrm{sim}} = \frac{| A \wedge B |}{|A| + |B| - | A \vee B |}
   * \f]
   */
  inline double bitvec_tanimoto(const Word *bitvec1, const Word *bitvec2, int bitCount1, int bitCount2, int numWords)
  {
    int andCount = 0;

    for (int i = 0; i < numWords; ++i) {
      Word andfp = bitvec1[i] & bitvec2[i];
#if __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
      andCount += __builtin_popcountl(andfp);
#else
      for(; andfp; andfp = andfp << 1)
        if(andfp < 0)
          ++andbits;
#endif
    }

    return static_cast<double>(andCount) / (bitCount1 + bitCount2 - andCount);
  }

  /**
   * Print a single bit vector word to std::cout.
   *
   * @param word The bit vector word to print.
   */
  inline void bitvec_print(Word word)
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

  /**
   * Print a bit vector to std::cout.
   *
   * @param bitvec The bit vector to print.
   * @param numWords The number of words for @p bitvec.
   * @param spaces If true, a space will be inserted between bit vector words.
   */
  inline void bitvec_print(const Word *bitvec, int numWords, bool spaces = true)
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

  /**
   * Write a bit vector's size to a STL file output stream. The size is written
   * as an unsigned long value.
   *
   * @param ofs The STL file output stream.
   * @param size The value to write to the stream.
   */
  inline void bitvec_write_size(std::ofstream &ofs, unsigned long size)
  {
    ofs.seekp(0);
    ofs.write(reinterpret_cast<const char*>(&size), sizeof(unsigned long));
  }

  /**
   * Read a bit vector's size from a STL file input stream.
   *
   * @param ifs The STL file input stream.
   *
   * @return The value read from the STL file input stream.
   */
  inline unsigned long bitvec_read_size(std::ifstream &ifs)
  {
    unsigned long size;
    ifs.seekg(0);
    ifs.read(reinterpret_cast<char*>(&size), sizeof(unsigned long));
    return size;
  }

  /**
   * Write a bit vector to a STL file output stream.
   *
   * @param ofs The STL file output stream.
   * @param bitvec The bit vector to write.
   * @param numWords The number of words for @p bitvec.
   */
  inline void bitvec_write(std::ofstream &ofs, const Word *bitvec, int numWords)
  {
    ofs.write(reinterpret_cast<const char*>(bitvec), sizeof(Word) * numWords);
  }

  /**
   * Read a bit vector from a STL file input stream.
   *
   * @param ifs The STL file input stream.
   * @param bitvec The bit vector to store the read data in.
   * @param numWords The number of words for @p bitvec.
   */
  inline void bitvec_read(std::ifstream &ifs, Word *bitvec, int numWords)
  {
    ifs.read(reinterpret_cast<char*>(bitvec), sizeof(Word) * numWords);
  }

}

#endif
