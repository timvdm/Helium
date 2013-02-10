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
#ifndef HELIUM_SIMILARITY_H
#define HELIUM_SIMILARITY_H

#include "fingerprints.h"

#include <numeric>

#ifdef HAVE_CPP11
#include <future>
#endif

namespace Helium {

  /**
   * @brief Brute force similarity search.
   *
   * This function computes the tanimoto coefficient between the query and
   * all fingerprints in the fingerprint sotrage. If the tanimoto is above the
   * specified threshold, it is added to the list of results.
   *
   * @param query The query fingerprint.
   * @param storage The fingerprints to search.
   * @param Tmin The minimum tanimoto score, must be in the range [0,1].
   * @param k When non-zero, this limits the number of requested results.
   */
  template<typename RowMajorFingerprintStorageType>
  std::vector<std::pair<unsigned int, double> > brute_force_similarity_search(const Word *query,
      RowMajorFingerprintStorageType &storage, double Tmin)
  {
    TIMER("brute_force_fimilarity_search():");
    int numWords = bitvec_num_words_for_bits(storage.numBits());

    std::vector<std::pair<unsigned int, double> > result;
    for (unsigned int i = 0; i < storage.numFingerprints(); ++i) {
      const Word *fingerprint = storage.fingerprint(i);
      double T = bitvec_tanimoto(query, fingerprint, numWords);
      if (T >= Tmin)
        result.push_back(std::make_pair(i, T));
    }

    return result;
  }

#ifdef HAVE_CPP11

  namespace impl {

    /**
     * Functor for threaded brute force similarity search.
     */
    template<typename RowMajorFingerprintStorageType>
    struct BruteForceSimilaritySearch
    {
      BruteForceSimilaritySearch(const Word *query_, RowMajorFingerprintStorageType &storage_,
          unsigned int begin_, unsigned int end_, double Tmin_) : query(query_),
          storage(storage_), begin(begin_), end(end_), Tmin(Tmin_)
      {
      }

      std::vector<std::pair<unsigned int, double> > operator()()
      {
        int numWords = bitvec_num_words_for_bits(storage.numBits());

        std::vector<std::pair<unsigned int, double> > result;
        for (unsigned int i = begin; i < end; ++i) {
          const Word *fingerprint = storage.fingerprint(i);
          double T = bitvec_tanimoto(query, fingerprint, numWords);
          if (T >= Tmin)
            result.push_back(std::make_pair(i, T));
        }

        return result;
      }

      const Word *query;
      RowMajorFingerprintStorageType &storage;
      unsigned int begin;
      unsigned int end;
      double Tmin;
    };

  }

  /**
   * @brief Brute force similarity search.
   *
   * This function computes the tanimoto coefficient between the query and
   * all fingerprints in the fingerprint sotrage. If the tanimoto is above the
   * specified threshold, it is added to the list of results.
   *
   * @param query The query fingerprint.
   * @param storage The fingerprints to search.
   * @param Tmin The minimum tanimoto score, must be in the range [0,1].
   */
  template<typename RowMajorFingerprintStorageType>
  std::vector<std::pair<unsigned int, double> > brute_force_similarity_search_threaded(const Word *query,
      RowMajorFingerprintStorageType &storage, double Tmin)
  {
    TIMER("brute_force_fimilarity_search_threaded():");

    unsigned int numFingerprints = storage.numFingerprints();
    unsigned int taskSize = numFingerprints / 4;

    auto future1 = std::async(impl::BruteForceSimilaritySearch<RowMajorFingerprintStorageType>(query, storage, 0, taskSize, Tmin));
    auto future2 = std::async(impl::BruteForceSimilaritySearch<RowMajorFingerprintStorageType>(query, storage, taskSize, 2 * taskSize, Tmin));
    auto future3 = std::async(impl::BruteForceSimilaritySearch<RowMajorFingerprintStorageType>(query, storage, 2 * taskSize, 3 * taskSize, Tmin));
    auto future4 = std::async(impl::BruteForceSimilaritySearch<RowMajorFingerprintStorageType>(query, storage, 3 * taskSize, numFingerprints, Tmin));

    std::vector<std::pair<unsigned int, double> > tmp;
    std::vector<std::pair<unsigned int, double> > result(future1.get());

    tmp = future2.get();
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(result));
    tmp = future3.get();
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(result));
    tmp = future4.get();
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(result));

    return result;
  }

#endif

  /**
   * @brief Similarity search fingerpint index.
   *
   */
  template<typename FingerprintStorageType>
  class SimilaritySearchIndex
  {
      /**
       * Base class for tree nodes in the kD-grid.
       */
      struct Node
      {
      };

      /**
       * An intermediate tree node in the kD-grid.
       */
      struct TreeNode : public Node
      {
        /**
         * Constructor.
         *
         * @param size The number of children.
         */
        TreeNode(int size)
        {
          children.resize(size, 0);
        }

        std::vector<Node*> children; //!< Child nodes
      };

      /**
       * Leaf node for tree in the kD-grid.
       */
      struct LeafNode : public Node
      {
        std::vector<unsigned int> fingerprints; //!< The fingerprint indices
      };

      /**
       * Get the bit count for the portion of the fingerprint at a certain
       * depth.
       */
      int bitCount(const Word *fingerprint, int depth) const
      {
        int beginBit = depth * (m_numBits / m_k);
        int endBit = (depth + 1) * (m_numBits / m_k);

        if (depth == m_k - 1)
          endBit = (depth + 1) * (m_numBits / m_k) + (m_numBits % ((depth + 1) * (m_numBits / m_k)));

        return bitvec_count(fingerprint, beginBit, endBit);
      }

      int childSize() const
      {
        return m_numBits / m_k + 1;
      }

      std::vector<unsigned int>& findLeaf(const Word *fingerprint)
      {
        int depth = 0;
        TreeNode *node = m_tree;

        while (depth < m_k) {
          int count = bitCount(fingerprint, depth);

          if (depth == m_k - 1) {
            LeafNode *leaf = static_cast<LeafNode*>(node->children[count]);
            if (!leaf) {
              leaf = new LeafNode;
              node->children[count] = leaf;
            }
            return leaf->fingerprints;
          }

          TreeNode *child = static_cast<TreeNode*>(node->children[count]);
          if (!child) {
            child = new TreeNode(childSize());
            node->children[count] = child;
          }

          node = child;
          ++depth;
        }

        UNREACHABLE_RETURN_REF(std::vector<unsigned int>);
      }

      void tanimotoDFS(const Word *fingerprint, double threshold, unsigned int maxResults,
          TreeNode *node, int depth, int *n_j, int *bitCounts, int bitCount,
          std::vector<std::pair<unsigned int, double> > &hits) const
      {
        // bound on the number of 1-bits in logical and (first part of bitstring)
        int A_i_min = 0;
        for (int i = 0; i < depth; ++i)
          A_i_min += std::min(bitCounts[i], n_j[i]);

        // bound on the number of 1-bits in logical or (first part of bitstring)
        int A_i_max = 0;
        for (int i = 0; i < depth; ++i)
          A_i_max += std::max(bitCounts[i], n_j[i]);

        // bound on last part of bitstring
        int A_i_last = 0;
        for (int i = depth + 1; i < m_k; ++i)
          A_i_last += bitCounts[i];

        // bitcount of current part
        int A_i = bitCounts[depth];
        // number of bins in the current part
        int N_i = childSize() - 1;

        int n_lower = std::max(static_cast<int>(std::ceil(threshold * (A_i_max + A_i + A_i_last) - (A_i_min + A_i_last))), 0);
        int n_upper = threshold == 0.0 ? N_i : std::min(static_cast<int>(std::floor((A_i_min + A_i + A_i_last - threshold * (A_i_max + A_i_last)) / threshold)), N_i);
        assert(n_lower >= 0);
        assert(n_upper <= N_i);

        ++depth;

        if (depth == m_k) {
          int bitCountB = std::accumulate(n_j, n_j + m_k, 0);

          for (int i = n_lower; i < n_upper + 1; ++i) {
            if (!node->children[i])
              continue;
            LeafNode *leaf = static_cast<LeafNode*>(node->children[i]);

            assert(leaf);
            for (std::size_t j = 0; j < leaf->fingerprints.size(); ++j) {
              double S = bitvec_tanimoto(fingerprint, m_storage.fingerprint(leaf->fingerprints[j]), bitCount, bitCountB + i, bitvec_num_words_for_bits(m_numBits));
              if (S >= threshold)
                hits.push_back(std::make_pair(leaf->fingerprints[j], S));
              if (hits.size() == maxResults)
                return;
            }
          }
          return;
        }

        for (int i = n_lower; i < n_upper + 1; ++i) {
          if (!node->children[i])
            continue;
          n_j[depth - 1] = i;
          tanimotoDFS(fingerprint, threshold, maxResults, static_cast<TreeNode*>(node->children[i]), depth, n_j, bitCounts, bitCount, hits);
          if (hits.size() == maxResults)
            return;
        }
      }

      void clearDFS(TreeNode *node, int depth)
      {
        ++depth;

        for (int i = 0; i < node->children.size(); ++i) {
          if (!node->children[i])
            continue;
          if (depth < m_k)
            clearDFS(static_cast<TreeNode*>(node->children[i]), depth); // recurse
          else
            delete node->children[i]; // delete leaf node
        }

        delete node;
      }

    public:
      SimilaritySearchIndex(const FingerprintStorageType &storage, int k)
          : m_storage(storage), m_k(k), m_numBits(storage.numBits())
      {
        TIMER("Loading SimilaritySearchIndex:");

        m_tree = new TreeNode(m_numBits / m_k + 1);

        for (unsigned int i = 0; i < m_storage.numFingerprints(); ++i)
          findLeaf(m_storage.fingerprint(i)).push_back(i);

      }

      ~SimilaritySearchIndex()
      {
        clearDFS(m_tree, 0);
        //delete m_tree;
      }

      std::vector<std::pair<unsigned int, double> > search(const Word *fingerprint, double threshold, unsigned int maxResults = 0) const
      {
        TIMER("SimilaritySearchIndex::search():");

        if (!maxResults)
          maxResults = std::numeric_limits<unsigned>::max();

        int count = 0;
        std::vector<int> n_j(m_k);
        std::vector<int> bitCounts(m_k);
        for (int i = 0; i < m_k; ++i) {
          int part = bitCount(fingerprint, i);
          bitCounts[i] = part;
          count += part;
        }

        std::vector<std::pair<unsigned int, double> > hits;
        tanimotoDFS(fingerprint, threshold, maxResults, m_tree, 0, &n_j[0], &bitCounts[0], count, hits);

        return hits;
      }

    private:
      const FingerprintStorageType &m_storage; //!< Fingerprint storage
      TreeNode *m_tree; //!< Tree root node
      int m_k; //!< Dimentionality of the kD-grid
      unsigned int m_numBits; //!< Number of bits in the fingerprint

  };

}

#endif
