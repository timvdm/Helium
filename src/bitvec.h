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

#include <Helium/contract.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>

namespace Helium {

  /**
   * @file bitvec.h
   * @brief Functions for working with bit vectors.
   */

  /**
   * @brief Type used for a bit vector word.
   */
  typedef unsigned long Word;

  /**
   * @brief The number of bits per bit vector word.
   *
   * Since the Word type is unsigned long, on a 64-bit system this is usually 64 bits.
   */
  const int BitsPerWord = sizeof(Word) * 8;

  /**
   * @brief Get the number of words needed to store @p numBits words.
   *
   * @param numBits The number of bits.
   *
   * @return The number of words needed to store @p numBits bits.
   */
  inline unsigned int bitvec_num_words_for_bits(unsigned int numBits)
  {
    if ((numBits % BitsPerWord) == 0)
      return numBits / BitsPerWord;
    return (numBits + BitsPerWord - numBits % BitsPerWord) / BitsPerWord;
  }

  /**
   * @brief Set all bits to 0.
   *
   * @pre The bitvec pointer must be valid.
   *
   * @post All the bits in the bit vector will be set to 0.
   *
   * @param bitvec The bit vector to zero.
   * @param numWords The number of words for @p bitvec.
   */
  inline void bitvec_zero(Word *bitvec, int numWords)
  {
    PRE(bitvec);
    for (int i = 0; i < numWords; ++i)
      bitvec[i] = 0;
  }

  /**
   * @brief Copy a bit vector.
   *
   * The copied vector's memory should be freed using delete [].
   *
   * @pre The bitvec pointer must be valid.
   *
   * @param bitvec The bit vector to copy.
   * @param numWords The number of words for @p bitvec.
   *
   * @return A pointer to the copied bit vector.
   */
  inline Word* bitvec_copy(const Word *bitvec, int numWords)
  {
    PRE(bitvec);
    Word *result = new Word[numWords];
    std::copy(bitvec, bitvec + numWords, result);
    return result;
  }

  /**
   * @brief Get the bit value for the specified bit index.
   *
   * This function does not check if the index is within the valid range.
   *
   * @pre The bitvec pointer must be valid.
   *
   * @param index The bit's index. The index should be less than the number of
   *        bits in the @p bitvec.
   * @param bitvec The bit vector.
   *
   * @return The value for the bit at the specified @p index.
   */
  inline bool bitvec_get(int index, const Word *bitvec)
  {
    PRE(bitvec);
    int word = index / (sizeof(Word) * 8);
    int offset = index % (sizeof(Word) * 8);
    Word bit = static_cast<Word>(1) << offset;
    return *(bitvec + word) & bit;
  }

  /**
   * @brief Set the bit at the specified bit index (to 1).
   *
   * This function does not check if the index is within the valid range.
   *
   * @pre The bitvec pointer must be valid.
   *
   * @param index The bit's index. The index should be less than the number of
   *        bits in the @p bitvec.
   * @param bitvec The bit vector.
   */
  inline void bitvec_set(int index, Word *bitvec)
  {
    PRE(bitvec);
    int word = index / (sizeof(Word) * 8);
    int offset = index % (sizeof(Word) * 8);
    Word bit = static_cast<Word>(1) << offset;
    *(bitvec + word) |= bit;
  }

  /**
   * @brief Reset the bit at the specified bit index (to 0).
   *
   * @pre The bitvec pointer must be valid.
   *
   * @param index The bit's index. The index should be less than the number of
   *        bits in the @p bitvec.
   * @param bitvec The bit vector.
   */
  inline void bitvec_reset(int index, Word *bitvec)
  {
    PRE(bitvec);
    int word = index / (sizeof(Word) * 8);
    int offset = index % (sizeof(Word) * 8);
    Word bit = static_cast<Word>(1) << offset;
    *(bitvec + word) &= ~bit;
  }

  /**
   * @brief Check if @p bitvec1 is a subset of @p bitvec2.
   *
   * Check if @p bitvec1 is a subset of @p bitvec2 (i.e. all bits set in
   * @p bitvec1 are also set in @p bitvec2). Both bitvectors should have the
   * same number of bits.
   *
   * @pre Both bitvec1 and bitvec2 pointers must be valid.
   *
   * @param bitvec1 The subset bit vector.
   * @param bitvec2 The superset bit vector.
   * @param numWords The number of words for @p bitvec1 and @p bitvec2.
   *
   * @return True if @p bitvec1 is a subset of @p bitvec2.
   */
  inline bool bitvec_is_subset_superset(const Word *bitvec1, const Word *bitvec2, int numWords)
  {
    PRE(bitvec1);
    PRE(bitvec2);
    for (int i = 0; i < numWords; ++i)
      if (bitvec1[i] & ~bitvec2[i])
        return false;
    return true;
  }

  /**
   * @brief Get the population count for a single word.
   *
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
#ifdef HAVE_POPCNT
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
   * @brief Get the population count for a bit vector.
   *
   * Get the bit count (i.e. number of bits set to 1 or the population count)
   * for a single bit vector word. When compiling with g++, the builtin function
   * wrapper for the POPCNT instruction is used.
   *
   * @pre The bitvec pointer must be valid.
   *
   * @param bitvec The bit vector.
   * @param numWords The number of words for @p bitvec.
   *
   * @return The bit count.
   */
  inline int bitvec_count(const Word *bitvec, int numWords)
  {
    PRE(bitvec);
#ifdef HAVE_POPCNT
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
   * @brief Get the population count for a partial bitvec.
   *
   * Get the bit count (i.e. number of bits set to 1 or the population count)
   * for a range of a bit vector.
   *
   * @pre The bitvec pointer must be valid.
   *
   * @param bitvec The bit vector.
   * @param begin The first bit index.
   * @param end The one-past-the-end bit index.
   *
   * @return The bit count for the range [begin,end).
   */
  inline int bitvec_count(const Word *bitvec, int begin, int end)
  {
    PRE(bitvec);
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
   * @brief Get the population count for the union of two bit vectors.
   *
   * Get the bit count (i.e. number of bits set to 1 or the population count)
   * for a range of a bit vector.
   *
   * @pre Both bitvec1 and bitvec2 pointers must be valid.
   *
   * @param bitvec The bit vector.
   * @param begin The first bit index.
   * @param end The one-past-the-end bit index.
   *
   * @return The bit count for the range [begin,end).
   */
  inline int bitvec_union_count(const Word *bitvec1, const Word *bitvec2, int numWords)
  {
    PRE(bitvec1);
    PRE(bitvec2);
    int count = 0;

    for (int i = 0; i < numWords; ++i) {
      Word andfp = bitvec1[i] & bitvec2[i];
#ifdef HAVE_POPCNT
      count += __builtin_popcountl(andfp);
#else
      for (; andfp; andfp = andfp << 1)
        if (andfp < 0)
          ++count;
#endif
    }

    return count;
  }

  /**
   * @brief Compute the Tanimoto coefficient.
   *
   * Compute the Tanimoto coefficient of difference between two bit vectors.
   * This function is slower than the Tanimoto function below since both the
   * union and intersection bit counts have to be computed.
   *
   * \f[
   *   T_{\mathrm{sim}} = \frac{| A \wedge B |}{| A \vee B |}
   * \f]
   *
   * @pre Both bitvec1 and bitvec2 pointers must be valid.
   *
   * @param bitvec1 The first bit vector (\f$A\f$).
   * @param bitvec2 The first bit vector (\f$B\f$).
   * @param numWords The number of words in the bit vectors.
   *
   * @return The Tanimoto coefficient of difference.
   */
  inline double bitvec_tanimoto(const Word *bitvec1, const Word *bitvec2, int numWords)
  {
    PRE(bitvec1);
    PRE(bitvec2);
    int andCount = 0;
    int orCount = 0;

    for (int i = 0; i < numWords; ++i) {
      Word andbv = bitvec1[i] & bitvec2[i];
      Word orbv = bitvec1[i] | bitvec2[i];
#ifdef HAVE_POPCNT
      andCount += __builtin_popcountll(andbv);
      orCount += __builtin_popcountll(orbv);
#else
      for (; andbv; andbv = andbv << 1)
        if(andbv < 0)
          ++andCount;
      for (; orbv; orbv = orbv << 1)
        if(orbv < 0)
          ++orCount;
#endif
    }

    return static_cast<double>(andCount) / orCount;
  }

  /**
   * @brief Compute the Tanimoto coefficient.
   *
   * Compute the Tanimoto coefficient of difference between two bit vectors.
   * When the bit counts of both bit vectors are known, the inclusion-exclusion
   * principle can be used to compute the Tanimoto coefficient faster.
   *
   * \f[
   *   T_{\mathrm{sim}} = \frac{| A \wedge B |}{|A| + |B| - | A \vee B |}
   * \f]
   *
   * @pre Both bitvec1 and bitvec2 pointers must be valid.
   *
   * @param bitvec1 The first bit vector (\f$A\f$).
   * @param bitvec2 The first bit vector (\f$B\f$).
   * @param bitCount1 The bit count for the first bit vector (\f$|A|\f$).
   * @param bitCount2 The bit count for the second bit vector (\f$|B|\f$).
   * @param numWords The number of words in the bit vectors.
   *
   * @return The Tanimoto coefficient of difference.
   */
  inline double bitvec_tanimoto(const Word *bitvec1, const Word *bitvec2, int bitCount1, int bitCount2, int numWords)
  {
    PRE(bitvec1);
    PRE(bitvec2);
    int andCount = bitvec_union_count(bitvec1, bitvec2, numWords);
    return static_cast<double>(andCount) / (bitCount1 + bitCount2 - andCount);
  }

  /**
   * @brief Compute the Cosine coefficient of difference between two bit vectors.
   *
   * \f[
   *   C_{\mathrm{sim}} = \frac{| A \wedge B |}{\sqrt{|A| |B|}}
   * \f]
   *
   * @pre Both bitvec1 and bitvec2 pointers must be valid.
   *
   * @param bitvec1 The first bit vector (\f$A\f$).
   * @param bitvec2 The first bit vector (\f$B\f$).
   * @param bitCount1 The bit count for the first bit vector (\f$|A|\f$).
   * @param bitCount2 The bit count for the second bit vector (\f$|B|\f$).
   * @param numWords The number of words in the bit vectors.
   *
   * @return The Cosine coefficient of difference.
   */
  inline double bitvec_cosine(const Word *bitvec1, const Word *bitvec2, int bitCount1, int bitCount2, int numWords)
  {
    PRE(bitvec1);
    PRE(bitvec2);
    int andCount = bitvec_union_count(bitvec1, bitvec2, numWords);
    return static_cast<double>(andCount) / std::sqrt(bitCount1 * bitCount2);
  }

  /**
   * @brief Compute the Cosine coefficient of difference between two bit vectors.
   *
   * \f[
   *   C_{\mathrm{sim}} = \frac{| A \wedge B |}{\sqrt{|A| |B|}}
   * \f]
   *
   * @pre Both bitvec1 and bitvec2 pointers must be valid.
   *
   * @param bitvec1 The first bit vector (\f$A\f$).
   * @param bitvec2 The first bit vector (\f$B\f$).
   * @param numWords The number of words in the bit vectors.
   *
   * @return The Cosine coefficient of difference.
   */
  inline double bitvec_cosine(const Word *bitvec1, const Word *bitvec2, int numWords)
  {
    PRE(bitvec1);
    PRE(bitvec2);
    int andCount = bitvec_union_count(bitvec1, bitvec2, numWords);
    return static_cast<double>(andCount) / std::sqrt(bitvec_count(bitvec1, numWords) * bitvec_count(bitvec2, numWords));
  }

  /**
   * @brief Compute the Hamming coefficient of difference between two bit vectors.
   *
   * \f[
   *   H_{\mathrm{sim}} = |A| + |B| - 2 | A \wedge B |
   * \f]
   *
   * @pre Both bitvec1 and bitvec2 pointers must be valid.
   *
   * @param bitvec1 The first bit vector (\f$A\f$).
   * @param bitvec2 The first bit vector (\f$B\f$).
   * @param bitCount1 The bit count for the first bit vector (\f$|A|\f$).
   * @param bitCount2 The bit count for the second bit vector (\f$|B|\f$).
   * @param numWords The number of words in the bit vectors.
   *
   * @return The Hamming coefficient of difference.
   */
  inline double bitvec_hamming(const Word *bitvec1, const Word *bitvec2, int bitCount1, int bitCount2, int numWords)
  {
    PRE(bitvec1);
    PRE(bitvec2);
    int andCount = bitvec_union_count(bitvec1, bitvec2, numWords);
    return bitCount1 + bitCount2 - 2 * andCount;
  }

  /**
   * @brief Compute the Hamming coefficient of difference between two bit vectors.
   *
   * \f[
   *   H_{\mathrm{sim}} = |A| + |B| - 2 | A \wedge B |
   * \f]
   *
   * @pre Both bitvec1 and bitvec2 pointers must be valid.
   *
   * @param bitvec1 The first bit vector (\f$A\f$).
   * @param bitvec2 The first bit vector (\f$B\f$).
   * @param numWords The number of words in the bit vectors.
   *
   * @return The Hamming coefficient of difference.
   */
  inline double bitvec_hamming(const Word *bitvec1, const Word *bitvec2, int numWords)
  {
    PRE(bitvec1);
    PRE(bitvec2);
    int andCount = bitvec_union_count(bitvec1, bitvec2, numWords);
    return bitvec_count(bitvec1, numWords) + bitvec_count(bitvec2, numWords) - 2 * andCount;
  }

  /**
   * @brief Compute the Russell-Rao coefficient of difference between two bit vectors.
   *
   * \f[
   *   R_{\mathrm{sim}} = \frac{| A \wedge B |}{m}
   * \f]
   *
   * @pre Both bitvec1 and bitvec2 pointers must be valid.
   *
   * @param bitvec1 The first bit vector (\f$A\f$).
   * @param bitvec2 The first bit vector (\f$B\f$).
   * @param numWords The number of words in the bit vectors (\f$m\f$ / (8 * sizeof(Word))).
   *
   * @return The Russell-Rao coefficient of difference.
   */
  inline double bitvec_russell_rao(const Word *bitvec1, const Word *bitvec2, int numWords)
  {
    PRE(bitvec1);
    PRE(bitvec2);
    int andCount = bitvec_union_count(bitvec1, bitvec2, numWords);
    return static_cast<double>(andCount) / (numWords * 8 * sizeof(Word));
  }

  /**
   * @brief Compute the Forbes coefficient of difference between two bit vectors.
   *
   * \f[
   *   F_{\mathrm{sim}} = \frac{| A \wedge B | m}{|A| |B|}
   * \f]
   *
   * @pre Both bitvec1 and bitvec2 pointers must be valid.
   *
   * @param bitvec1 The first bit vector (\f$A\f$).
   * @param bitvec2 The first bit vector (\f$B\f$).
   * @param bitCount1 The bit count for the first bit vector (\f$|A|\f$).
   * @param bitCount2 The bit count for the second bit vector (\f$|B|\f$).
   * @param numWords The number of words in the bit vectors (\f$m\f$ / (8 * sizeof(Word))).
   *
   * @return The Forbes coefficient of difference.
   */
  inline double bitvec_forbes(const Word *bitvec1, const Word *bitvec2, int bitCount1, int bitCount2, int numWords)
  {
    PRE(bitvec1);
    PRE(bitvec2);
    int andCount = bitvec_union_count(bitvec1, bitvec2, numWords);
    return static_cast<double>(andCount * numWords * 8 * sizeof(Word)) / (bitCount1 * bitCount2);
  }

  /**
   * @brief Compute the Forbes coefficient of difference between two bit vectors.
   *
   * \f[
   *   F_{\mathrm{sim}} = \frac{| A \wedge B | m}{|A| |B|}
   * \f]
   *
   * @pre Both bitvec1 and bitvec2 pointers must be valid.
   *
   * @param bitvec1 The first bit vector (\f$A\f$).
   * @param bitvec2 The first bit vector (\f$B\f$).
   * @param numWords The number of words in the bit vectors (\f$m\f$ / (8 * sizeof(Word))).
   *
   * @return The Forbes coefficient of difference.
   */
  inline double bitvec_forbes(const Word *bitvec1, const Word *bitvec2, int numWords)
  {
    PRE(bitvec1);
    PRE(bitvec2);
    int andCount = bitvec_union_count(bitvec1, bitvec2, numWords);
    return static_cast<double>(andCount * numWords * 8 * sizeof(Word)) / (bitvec_count(bitvec1, numWords) * bitvec_count(bitvec2, numWords));
  }

  /**
   * @brief Print a single bit vector word to std::cout.
   *
   * @post The word's bits will be printed to stdout.
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
   * @brief Print a bit vector to std::cout.
   *
   * @post The bit vector's bits will be printed to stdout.
   *
   * @param bitvec The bit vector to print.
   * @param numWords The number of words for @p bitvec.
   * @param spaces If true, a space will be inserted between bit vector words.
   */
  inline void bitvec_print(const Word *bitvec, int numWords, bool spaces = true)
  {
    PRE(bitvec);
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
   * @brief Write a bit vector's size to a STL file output stream.
   *
   * The size is written as an unsigned long value.
   *
   * @post The specified size will be written to output stream.
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
   * @brief Read a bit vector's size from a STL file input stream.
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
   * @brief Write a bit vector to a STL file output stream.
   *
   * @post The bit vector will be written to the output stream.
   *
   * @param ofs The STL file output stream.
   * @param bitvec The bit vector to write.
   * @param numWords The number of words for @p bitvec.
   */
  inline void bitvec_write(std::ofstream &ofs, const Word *bitvec, int numWords)
  {
    PRE(bitvec);
    ofs.write(reinterpret_cast<const char*>(bitvec), sizeof(Word) * numWords);
  }

  /**
   * @brief Read a bit vector from a STL file input stream.
   *
   * @param ifs The STL file input stream.
   * @param bitvec The bit vector to store the read data in.
   * @param numWords The number of words for @p bitvec.
   */
  inline void bitvec_read(std::ifstream &ifs, Word *bitvec, int numWords)
  {
    PRE(bitvec);
    ifs.read(reinterpret_cast<char*>(bitvec), sizeof(Word) * numWords);
  }

  namespace impl {

    inline int hex_to_dec(char digit)
    {
      switch (digit) {
        case '0':
          return 0;
        case '1':
          return 1;
        case '2':
          return 2;
        case '3':
          return 3;
        case '4':
          return 4;
        case '5':
          return 5;
        case '6':
          return 6;
        case '7':
          return 7;
        case '8':
          return 8;
        case '9':
          return 9;
        case 'a':
        case 'A':
          return 10;
        case 'b':
        case 'B':
          return 11;
        case 'c':
        case 'C':
          return 12;
        case 'd':
        case 'D':
          return 13;
        case 'e':
        case 'E':
          return 14;
        case 'f':
        case 'F':
          return 15;
      }

      assert(0);
      return 0;
    }

    inline char dec_to_hex(int value)
    {
      assert(value >= 0 && value < 16);

      switch (value) {
        case 0:
          return '0';
        case 1:
          return '1';
        case 2:
          return '2';
        case 3:
          return '3';
        case 4:
          return '4';
        case 5:
          return '5';
        case 6:
          return '6';
        case 7:
          return '7';
        case 8:
          return '8';
        case 9:
          return '9';
        case 10:
          return 'a';
        case 11:
          return 'b';
        case 12:
          return 'c';
        case 13:
          return 'd';
        case 14:
          return 'e';
        case 15:
          return 'f';
      }

      assert(0);
      return '0';
    }

  }

  /**
   * @brief Convert a hexadecimal string to a bit vector.
   *
   * Both upper case and lower case letters can be used for the hexadecimal
   * bit vector.
   *
   * @pre There must be an even number of characters in the @p hex string
   *      (i.e. half-byte or nibbles are not allowed). There may not be more
   *      bytes in the hexadecimal string than the size of the bitvec. The
   *      @p bitvec pointer must be valid.
   * @code
   * PRE(hex.size() / 16 <= numWords);
   * @endcode
   *
   * @post @p bitvec will contain the bits specified by the hexadecimal string.
   *
   * @param hex The hexadecimal bit vector.
   * @param bitvec The bit vector to store the result.
   * @param numWords The number of words for @p bitvec.
   */
  inline void hex_to_bitvec(const std::string &hex, Word *bitvec, int numWords)
  {
    PRE(bitvec);
    PRE((hex.size() % 2) == 0);
    PRE(hex.size() / 16 <= numWords);

    bitvec_zero(bitvec, numWords);

    unsigned char *fp = reinterpret_cast<unsigned char*>(bitvec);

    for (std::size_t i = 0; i < hex.size() / 2; ++i) {
      fp[i] = impl::hex_to_dec(hex[2 * i]) << 4;
      fp[i] += impl::hex_to_dec(hex[2 * i + 1]);
    }
  }

  /**
   * @brief Convert a bit vector to a hexadecimal string.
   *
   * Both upper case and lower case letters can be used for the hexadecimal
   * bit vector.
   *
   * @pre The @p bitvec pointer must be valid.
   *
   * @param bitvec The bit vector to store the result.
   * @param numWords The number of words for @p bitvec.
   *
   * @return The hexadecimal string containing the bits specified by the bit vector.
   */
  inline std::string bitvec_to_hex(const Word *bitvec, int numWords)
  {
    PRE(bitvec);
    std::size_t bytes = numWords * sizeof(Word);
    const unsigned char *fp = reinterpret_cast<const unsigned char*>(bitvec);

    std::stringstream ss;
    for (std::size_t i = 0; i < bytes; ++i)
      ss << impl::dec_to_hex((fp[i] & 240) >> 4) << impl::dec_to_hex(fp[i] & 15);

    return ss.str();
  }

}

#endif
