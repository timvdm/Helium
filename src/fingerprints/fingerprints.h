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
#ifndef HELIUM_FINGERPRINTS_H
#define HELIUM_FINGERPRINTS_H

#include <Helium/bitvec.h>
#include <Helium/algorithms/enumeratepaths.h>
#include <Helium/algorithms/enumeratesubgraphs.h>
#include <Helium/substructure.h>
#include <Helium/algorithms/extendedconnectivities.h>
#include <Helium/algorithms/canonical.h>

#include <stdexcept>

#include <boost/functional/hash.hpp>

namespace Helium {

  /**
   * Compute the path-based fingerprint for the specified molecule. All paths
   * in the molecular graph will be enumerated upto the specified size. For each
   * path, a canonical code is generated which is hashed using the @p hashPrime
   * number to set a bit in the @p fingerprint corresponding to that path.
   *
   * @param mol The molecule which models the MoleculeConcept.
   * @param fingerprint Pointer to the fingerprint memory. This memory must be
   *        the correct size (see @p numWords).
   * @param size Maximum number of atoms in the paths.
   * @param numWords The number of words the fingerprint has. Each word has
   *        8 * sizeof(Word) bits which is usually 64 bit. This means the
   *        default value 16 results in fingerprints of 1024 bits.
   * @param hashPrime A prime number to hash the paths so they will fit in the
   *        fingerprint. The largest prime, less than or equal to the number of
   *        bits in the fingerprint is ideal.
   */
  template<typename MoleculeType>
  void path_fingerprint(MoleculeType &mol, Word *fingerprint, int size = 7, int numWords = 16, int hashPrime = 1021)
  {
    assert(hashPrime <= numWords * sizeof(Word) * 8);
    boost::hash<std::vector<unsigned long> > hash;
    // set all bits to 0
    bitvec_zero(fingerprint, numWords);
    // enumerate the paths
    std::vector<std::vector<unsigned int> > paths = enumerate_paths(mol, size);
    // set the bits
    for (std::size_t i = 0; i < paths.size(); ++i) {
      std::vector<bool> atoms(num_atoms(mol));
      std::vector<bool> bonds(num_bonds(mol));

      // set bits for atoms/bonds in the path
      for (std::size_t j = 0; j < paths[i].size(); ++j) {
        atoms[paths[i][j]] = true;
        if (j + 1 < paths[i].size())
          bonds[get_index(mol, get_bond(mol, get_atom(mol, paths[i][j]), get_atom(mol, paths[i][j + 1])))] = true;
      }

      // create path molecule
      Substructure<MoleculeType> substruct(mol, atoms, bonds);
      // compute symmetry classes
      std::vector<unsigned long> symmetry = extended_connectivities(substruct);
      // canonicalize the path
      std::vector<unsigned long> code = canonicalize(substruct, symmetry, AtomElementAttribute(), BondOrderAttribute()).second;
      // set the bit for the hashed canonical code modulo the hash prime.
      bitvec_set(hash(code) % hashPrime, fingerprint);
    }
  }

  namespace impl {

    template<typename MoleculeType>
    struct SubgraphsFingerprint
    {
      /**
       * The subgraph enumeration callback. This functor computes the
       * fingerprint by canonicalizing each subgraph and setting the bit for
       * the hashed value of the subgraph's canonical code modulo the hash
       * prime.
       */
      struct EnumerateSubgraphsCallback
      {
        /**
         * Constructor.
         *
         * @param mol_ The molecule,
         * @param fp The fingerprint bit vector.
         * @param words The number of words for the @p fp bit vector.
         * @param prime The prime to use for taking the modulo (e.g. the largest
         *              prime smaller than the number of bits in the fingerprint.
         */
        EnumerateSubgraphsCallback(MoleculeType &mol_, Word *fp, int words, int prime)
            : mol(mol_), fingerprint(fp), numWords(words), hashPrime(prime)
        {
          bitvec_zero(fingerprint, numWords);
        }

        /**
         * Callback function, gets called for every subgraph found in the molecule.
         */
        void operator()(const Subgraph &subgraph)
        {
          // create the subgraph molecule
          Substructure<MoleculeType> substruct(mol, subgraph.atoms, subgraph.bonds);
          // compute symmetry classes
          std::vector<unsigned long> symmetry = extended_connectivities(substruct);
          // canonicalize the subgraph
          std::vector<unsigned long> code = canonicalize(substruct, symmetry, AtomElementAttribute(), BondOrderAttribute()).second;
          // set the bit for the hashed canonical code modulo the hash prime.
          bitvec_set(m_hash(code) % hashPrime, fingerprint);
        }

        MoleculeType &mol; //!< The molecule
        Word *fingerprint; //!< The fingerprint bit vector
        boost::hash<std::vector<unsigned long> > m_hash; //!< The hash function
        int numWords; //!< The number of words for the fingerprint bit vector
        int hashPrime; //!< The modulo prime number
      };

      SubgraphsFingerprint(MoleculeType &mol, Word *fp, int size, bool trees, int numWords, int hashPrime)
          : callback(mol, fp, numWords, hashPrime)
      {
        // enumerate subgraphs
        enumerate_subgraphs(mol, callback, size, trees);
      }

      EnumerateSubgraphsCallback callback; //!< Subgraph enumerator callback
    };

  }

  /**
   * Compute the tree-based fingerprint for the specified molecule. All trees
   * in the molecular graph will be enumerated upto the specified size. For each
   * tree, a canonical code is generated which is hashed using the @p hashPrime
   * number to set a bit in the @p fingerprint corresponding to that tree.
   *
   * @param mol The molecule which models the MoleculeConcept.
   * @param fingerprint Pointer to the fingerprint memory. This memory must be
   *        the correct size (see @p numWords).
   * @param size Maximum number of atoms in the trees.
   * @param numWords The number of words the fingerprint has. Each word has
   *        8 * sizeof(Word) bits which is usually 64 bit. This means the
   *        default value 16 results in fingerprints of 1024 bits.
   * @param hashPrime A prime number to hash the trees so they will fit in the
   *        fingerprint. The largest prime, less than or equal to the number of
   *        bits in the fingerprint is ideal.
   */
  template<typename MoleculeType>
  void tree_fingerprint(MoleculeType &mol, Word *fingerprint, int size = 7, int numWords = 16, int hashPrime = 1021)
  {
    assert(hashPrime <= numWords * sizeof(Word) * 8);
    impl::SubgraphsFingerprint<MoleculeType>(mol, fingerprint, size, true, numWords, hashPrime);
  }

  /**
   * Compute the subgraph-based fingerprint for the specified molecule. All subgraphs
   * in the molecular graph will be enumerated upto the specified size. For each
   * subgraph, a canonical code is generated which is hashed using the @p hashPrime
   * number to set a bit in the @p fingerprint corresponding to that subgraph.
   *
   * @param mol The molecule which models the MoleculeConcept.
   * @param fingerprint Pointer to the fingerprint memory. This memory must be
   *        the correct size (see @p numWords).
   * @param size Maximum number of atoms in the subgraphs.
   * @param numWords The number of words the fingerprint has. Each word has
   *        8 * sizeof(Word) bits which is usually 64 bit. This means the
   *        default value 16 results in fingerprints of 1024 bits.
   * @param hashPrime A prime number to hash the subgraphs so they will fit in the
   *        fingerprint. The largest prime, less than or equal to the number of
   *        bits in the fingerprint is ideal.
   */
  template<typename MoleculeType>
  void subgraph_fingerprint(MoleculeType &mol, Word *fingerprint, int size = 7, int numWords = 16, int hashPrime = 1021)
  {
    assert(hashPrime <= numWords * sizeof(Word) * 8);
    impl::SubgraphsFingerprint<MoleculeType>(mol, fingerprint, size, false, numWords, hashPrime);
  }

}

#endif
