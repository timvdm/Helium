#include "../src/molecule.h"
#include "../src/fileio.h"

#include "../src/fingerprints.h"

#include <numeric>
#include <boost/functional/hash.hpp>

#include "args.h"

using namespace Helium;

int main(int argc, char**argv)
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [options] <method> <in_file> <out_file>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Methods:" << std::endl;
    std::cerr << "    -paths        Create hashed fingerprints from paths" << std::endl;
    std::cerr << "    -trees        Create hashed fingerprints from trees" << std::endl;
    std::cerr << "    -subgraphs    Create hashed fingerprints from subgraphs" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "    -k <number>   The maximum size of the path/tree/subgraph (default is 7)" << std::endl;

    std::cerr << std::endl;
    return 0;
  }

  ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("in_file", "out_file"));
  std::string inFile = args.GetArgString("in_file");
  std::string outFile = args.GetArgString("out_file");

  // open molecule file
  MoleculeFile file(inFile);
  HeMol mol;

  std::map<int, std::vector<unsigned int> > sizes;

  // process molecules
  while (file.read_molecule(mol)) {
    sizes[num_atoms(mol)].push_back(file.current());
  }

  //std::cout << "num_sizes: " << sizes.size() << std::endl;

  // open filter file
  std::ofstream ofs(outFile.c_str(), std::ios_base::out | std::ios_base::binary);

  // write number of molecules
  write32(ofs, file.numMolecules());
  // write number of sizes
  write32(ofs, sizes.size());
  // write sizes
  for (std::map<int, std::vector<unsigned int> >::iterator i = sizes.begin(); i != sizes.end(); ++i)
    write32(ofs, i->first);
  // write filters
  Word *filter = new Word[bitvec_num_words_for_bits(file.numMolecules())];
  for (std::map<int, std::vector<unsigned int> >::iterator i = sizes.begin(); i != sizes.end(); ++i) {
    //std::cout << i->first << " atoms: " << i->second.size() << " molecules" << std::endl;
    bitvec_zero(filter, bitvec_num_words_for_bits(file.numMolecules()));
    std::cout << i->first << ": ";
    for (std::size_t j = 0; j < i->second.size(); ++j) {
      bitvec_set(i->second[j], filter);
      std::cout << i->second[j] << " ";
    }
    std::cout << std::endl;
    //std::cout << i->first << " -> " << (static_cast<unsigned int>(ofs.tellp()) - 4 * (2 + sizes.size())) / sizeof(Word) << std::endl;
    ofs.write(reinterpret_cast<const char*>(filter), sizeof(Word) * bitvec_num_words_for_bits(file.numMolecules()));
  }

  delete [] filter;
}
