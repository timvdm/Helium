#include "../src/molecule.h"
#include "../src/fileio.h"

#include "../src/fingerprints.h"
#include "../src/fileio/fingerprints.h"

#include <numeric>
#include <boost/functional/hash.hpp>

#include "args.h"

using namespace Helium;

const int NUM_BITS = 1024;
const int PRIME = 1021;
const int NUM_WORDS = NUM_BITS / (8 * sizeof(Word));


int main(int argc, char**argv)
{
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " <method> <query> <molecule_file>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Methods:" << std::endl;
    std::cerr << "    -paths        Create hashed fingerprints from paths" << std::endl;
    std::cerr << "    -trees        Create hashed fingerprints from trees" << std::endl;
    std::cerr << "    -subgraphs    Create hashed fingerprints from subgraphs" << std::endl;
    return 0;
  }

  ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("method", "query", "molecule_file"));
  std::string methodString = args.GetArgString("method");
  std::string query = args.GetArgString("query");
  std::string filename = args.GetArgString("molecule_file");

  enum Method {
    PathsMethod,
    TreesMethod,
    SubgraphsMethod
  };

  Method method = PathsMethod;
  if (methodString == "-paths")
    method = PathsMethod;
  else if (methodString == "-trees")
    method = TreesMethod;
  else if (methodString == "-subgraphs")
    method = SubgraphsMethod;
  else {
    std::cerr << "Method \"" << methodString << "\" not recognised" << std::endl;
    return -1;
  }

  // convert query to fingerprint
  HeMol mol;
  read_smiles(query, mol);
  Word queryFingerprint[NUM_WORDS];
  switch (method) {
    case PathsMethod:
      path_fingerprint(mol, queryFingerprint, 7, NUM_WORDS, PRIME);
      break;
    case TreesMethod:
      tree_fingerprint(mol, queryFingerprint, 7, NUM_WORDS, PRIME);
      break;
    case SubgraphsMethod:
      subgraph_fingerprint(mol, queryFingerprint, 7, NUM_WORDS, PRIME);
      break;
  }

//  print(queryFingerprint, NUM_WORDS);

  // open fingerprint file
  //InvertedFingerprintFile file(filename);
  InvertedFingerprintFileCached file(filename);

  // process fingerprints
  unsigned int hits = 0;
  Word *screen = file.allocate_result();
  file.search(queryFingerprint, screen);

  for (unsigned int i = 0; i < file.num_fingerprints(); ++i)
    if (bitvec_get(i, screen))
      ++hits;
  /*
  */
 
  /*
  // open fingerprint file
  FingerprintFile file(filename);
  // process fingerprints
  unsigned int hits = 0;
  Word fingerprint[NUM_WORDS];
  while (file.read_fingerprint(fingerprint)) {
    //print(fingerprint, NUM_WORDS);
    if (is_subset_superset(queryFingerprint, fingerprint, NUM_WORDS))
      ++hits;
  }
  */

  std::cout << "Hits: " << hits << "/" << file.num_fingerprints() << std::endl;
}
