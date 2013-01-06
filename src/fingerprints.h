#ifndef HELIUM_FINGERPRINTS_H
#define HELIUM_FINGERPRINTS_H

#include "../src/bitvec.h"
#include "../src/enumeratepaths.h"
#include "../src/enumeratesubgraphs.h"
#include "../src/substructure.h"
#include "../src/extendedconnectivities.h"
#include "../src/canonical.h"

#include <boost/functional/hash.hpp>

namespace Helium {

  /**
   * Calculate the path-based fingerprint for the specified molecule. All paths
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
  void path_fingerprint(MoleculeType *mol, Word *fingerprint, int size = 7, int numWords = 16, int hashPrime = 1021)
  {
    boost::hash<std::vector<unsigned long> > hash;
    zero(fingerprint, numWords);
    std::vector<std::vector<unsigned int> > paths = enumerate_paths(mol, size);
    for (std::size_t i = 0; i < paths.size(); ++i) {
      std::vector<unsigned long> canonicalCode = canonicalize_path(mol, paths[i]).second;
      set(hash(canonicalCode) % hashPrime, fingerprint, numWords);
    }
  }

  namespace impl {

    template<typename MoleculeType>
    struct SubgraphsFingerprint
    {
      struct EnumerateSubgraphsCallback
      {
        EnumerateSubgraphsCallback(MoleculeType *mol_, Word *fp, int words, int prime)
            : mol(mol_), fingerprint(fp), numWords(words), hashPrime(prime)
        {
          zero(fingerprint, numWords);
        }

        void operator()(const Subgraph &subgraph)
        {
          Substructure substruct(mol, subgraph.atoms, subgraph.bonds);

          std::vector<unsigned long> symmetry = extended_connectivities(&substruct);
          std::vector<unsigned long> canonicalCode = canonicalize(&substruct, symmetry).second;

          set(m_hash(canonicalCode) % hashPrime, fingerprint, numWords);
        }

        MoleculeType *mol;
        Word *fingerprint;
        boost::hash<std::vector<unsigned long> > m_hash;
        int numWords;
        int hashPrime;
      };

      SubgraphsFingerprint(MoleculeType *mol, Word *fp, int size, bool trees, int numWords, int hashPrime)
          : callback(mol, fp, numWords, hashPrime)
      {
        // enumerate subgraphs
        enumerate_subgraphs(mol, callback, size, trees);
      }

      EnumerateSubgraphsCallback callback;
    };

  }

  template<typename MoleculeType>
  void tree_fingerprint(MoleculeType *mol, Word *fingerprint, int size = 7, int numWords = 16, int hashPrime = 1021)
  {
    impl::SubgraphsFingerprint<MoleculeType>(mol, fingerprint, size, true, numWords, hashPrime);
  }

  template<typename MoleculeType>
  void subgraph_fingerprint(MoleculeType *mol, Word *fingerprint, int size = 7, int numWords = 16, int hashPrime = 1021)
  {
    impl::SubgraphsFingerprint<MoleculeType>(mol, fingerprint, size, false, numWords, hashPrime);
  }
  
}

#endif
