#include "../src/molecule.h"
#include "../src/fileio.h"

#include "../src/fingerprints.h"

#include <numeric>
#include <boost/functional/hash.hpp>

#include "args.h"

using namespace Helium;

const int NUM_BITS = 1024;
const int PRIME = 1021;
const int BITS_PER_WORD = 8 * sizeof(Word);
const int NUM_WORDS = NUM_BITS / BITS_PER_WORD;


int main(int argc, char**argv)
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " [options] <method> <in_file> <out_file>" << std::endl;
    return 0;
  }

  ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("in_file", "out_file"));
  std::string inFile = args.GetArgString("in_file");
  std::string outFile = args.GetArgString("out_file");

  // get number of fingerprints
  unsigned int numFingerprints = FingerprintFile(inFile).num_fingerprints();


  // open input fingerprint file
  FingerprintFile file(inFile);
  // open output fingerprint file
  InvertedFingerprintOutputFile invfile(NUM_BITS, numFingerprints, outFile);

  // process fingerprints
  Word fingerprint[NUM_WORDS];
  while (file.read_fingerprint(fingerprint)) {
    invfile.write(fingerprint);
  }








  /*
  unsigned int numBits = numFingerprints + numFingerprints % BITS_PER_WORD;
  assert(numBits % BITS_PER_WORD == 0);

  // open inverted index file
  std::ofstream ofs(outFile.c_str(), std::ios_base::out | std::ios_base::binary);
  write32(ofs, numFingerprints);

  Word **buffer = new Word*[BITS_PER_WORD];
  for (int i = 0; i < BITS_PER_WORD; ++i)
    buffer[i] = new Word[numBits];

  for (unsigned int i = 0; i < NUM_WORDS;  ++i) {
    std::cout << "bit " << i * BITS_PER_WORD << "-" << (i + 1) * BITS_PER_WORD << std::endl;
    // open molecule file
    FingerprintFile file(inFile);

    // process fingerprints
    Word fingerprint[NUM_WORDS];
    while (file.read_fingerprint(fingerprint)) {
      int start_index = i * BITS_PER_WORD;

      for (int j = 0; j < BITS_PER_WORD; ++j) {
        if (get(start_index + j, fingerprint, NUM_WORDS))
          set(file.current() - 1, buffer[j], NUM_WORDS);
      }
    }

    for (unsigned int j = 0; j < BITS_PER_WORD;  ++j) // for each bit
      for (int k = 0; k < numBits / BITS_PER_WORD; ++k) // for each word
        write64(ofs, buffer[j][k]);
  }
  */
}
