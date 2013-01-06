#include "../src/molecule.h"
#include "../src/fileio.h"
#include "../src/bitvec.h"

#include "../src/enumeratesubgraphs.h"
#include "../src/substructure.h"
#include "../src/extendedconnectivities.h"
#include "../src/canonical.h"

#include <numeric>

#include <boost/timer/timer.hpp>

#include "args.h"

const boost::timer::nanosecond_type one_milisecond(1000000L);

using namespace Helium;

std::string get_filename(const std::string &filename)
{
  return DATADIR + std::string("/") + filename;
}

/**
 *
 * Benchmark reading *.hem files.
 *
 */

void benchmark_read(const std::string &filename)
{
  boost::timer::cpu_timer timer;

  // open molecule file
  MoleculeFile file(get_filename(filename));
  Molecule mol;

  // read molecules
  while (file.read_molecule(mol));

  boost::timer::cpu_times elapsed = timer.elapsed();
  boost::timer::nanosecond_type ns = elapsed.system + elapsed.user;
  std::cout << "read_file (" << file.current() << " molecules): " << ns / one_milisecond << " ms" << std::endl;
}

/**
 *
 * Benchmark subgraph enumeration.
 *
 */

template<typename MoleculeType>
struct EnumerateSubgraphsCallback
{
  EnumerateSubgraphsCallback() : count(0) {}
  void operator()(const Subgraph &subgraph) { ++count; }
  int count;
};

void benchmark_enumerate_subgraphs(const std::string &filename, int size = 7)
{
  boost::timer::cpu_timer timer;

  // open molecule file
  MoleculeFile file(get_filename(filename));
  Molecule mol;

  // read molecules
  while (file.read_molecule(mol)) {
    EnumerateSubgraphsCallback<Molecule> callback;
    enumerate_subgraphs(&mol, callback, size);
  }

  boost::timer::cpu_times elapsed = timer.elapsed();
  boost::timer::nanosecond_type ns = elapsed.system + elapsed.user;
  std::cout << "subgraph_enumeration (" << file.current() << " molecules): " << ns / one_milisecond << " ms" << std::endl;
}









int main(int argc, char**argv)
{
  /*
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <molecule_file>" << std::endl;
    return 0;
  }
  */

  benchmark_read("1M.hem");
  benchmark_enumerate_subgraphs("10K.hem");

  /*
  ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("in_file", "out_file"));
  std::string inFile = args.GetArgString("in_file");
  std::string outFile = args.GetArgString("out_file");


  // open index file
  std::ofstream ofs(outFile.c_str(), std::ios_base::out | std::ios_base::binary);
  write32(ofs, 0);

  // open molecule file
  MoleculeFile file(inFile);
  Molecule mol;

  std::vector<int> bitCounts;

  // process molecules
  while (file.read_molecule(mol)) {
    if ((file.current() % 100) == 0)
      std::cout << file.current() << std::endl;

    // enumerate subgraphs
    EnumerateSubgraphsCallback<Molecule> callback(&mol);
    enumerate_subgraphs(&mol, callback, 7);

    int bitCount = bit_count(callback.fingerprint, NUM_WORDS);
    bitCounts.push_back(bitCount);

    std::cout << "Molecule " << file.current() << ": " <<  bitCount << " bits " << std::endl;

    for (int i = 0; i < NUM_WORDS; ++i)
      write64(ofs, callback.fingerprint[i]);
  }

  unsigned int sum = std::accumulate(bitCounts.begin(), bitCounts.end(), 0);

  std::cout << "Average bits: " << sum / static_cast<double>(bitCounts.size()) << std::endl;


  // write number of molecules
  ofs.seekp(0);
  write32(ofs, file.numMolecules());
  */
}
