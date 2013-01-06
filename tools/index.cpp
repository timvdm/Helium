#include "../src/molecule.h"
#include "../src/fileio.h"

#include "../src/fingerprints.h"

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
    std::cerr << "Usage: " << argv[0] << " <method> <in_file> <out_file>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Methods:" << std::endl;
    std::cerr << "    -paths        Create hashed fingerprints from paths" << std::endl;
    std::cerr << "    -trees        Create hashed fingerprints from trees" << std::endl;
    std::cerr << "    -subgraphs    Create hashed fingerprints from subgraphs" << std::endl;
    return 0;
  }

  ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("method", "in_file", "out_file"));
  std::string methodString = args.GetArgString("method");
  std::string inFile = args.GetArgString("in_file");
  std::string outFile = args.GetArgString("out_file");

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


  // open index file
  std::ofstream ofs(outFile.c_str(), std::ios_base::out | std::ios_base::binary);
  write32(ofs, 0);

  // open molecule file
  MoleculeFile file(inFile);
  Molecule mol;

  Word fingerprint[NUM_WORDS];
  std::vector<int> bitCounts;

  // process molecules
  while (file.read_molecule(mol)) {
    if ((file.current() % 100) == 0)
      std::cout << file.current() << std::endl;

    switch (method) {
      case PathsMethod:
        path_fingerprint(&mol, fingerprint, 7, NUM_WORDS, PRIME);
        break;
      case TreesMethod:
        tree_fingerprint(&mol, fingerprint, 7, NUM_WORDS, PRIME);
        break;
      case SubgraphsMethod:
        subgraph_fingerprint(&mol, fingerprint, 7, NUM_WORDS, PRIME);
        break;
    }

    int bitCount = bit_count(fingerprint, NUM_WORDS);
    bitCounts.push_back(bitCount);

    std::cout << "Molecule " << file.current() << ": " <<  bitCount << " bits " << std::endl;

    for (int i = 0; i < NUM_WORDS; ++i)
      write64(ofs, fingerprint[i]);
  }

  unsigned int sum = std::accumulate(bitCounts.begin(), bitCounts.end(), 0);

  std::cout << "Average bits: " << sum / static_cast<double>(bitCounts.size()) << std::endl;


  // write number of molecules
  ofs.seekp(0);
  write32(ofs, file.numMolecules());

}
