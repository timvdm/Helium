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
#include "tool.h"

#include <Helium/fileio/fingerprints.h>
#include <Helium/fingerprints/fingerprints.h>

#include <numeric> // std::accumulate

#include "args.h"

using namespace Helium;

namespace Helium {

  /**
   * Tool for folding fingerprint indexes.
   */
  class FoldTool : public HeliumTool
  {
    public:
      /**
       * Perform tool action.
       */
      int run(int argc, char **argv)
      {
        ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("bits", "in_file", "out_file"));
        // required arguments
        const int bits = args.GetArgInt("bits");
        const int words = bits / (8 * sizeof(Word));
        const int prime = previous_prime(bits);
        std::string inFile = args.GetArgString("in_file");
        std::string outFile = args.GetArgString("out_file");

        // open input file
        InMemoryRowMajorFingerprintStorage inputFile;
        inputFile.load(inFile);

        // open output file
        RowMajorFingerprintOutputFile outputFile(outFile, bits);

        // allocate bit vector
        Word *folded = new Word[words];
        // keep track of bit counts
        std::vector<int> bitCounts;

        // process molecules
        for (unsigned int i = 0; i < inputFile.numFingerprints(); ++i) {

          const Word *unfolded = inputFile.fingerprint(i);

          bitvec_zero(folded, words);
          for (int j = 0; j < inputFile.numBits(); ++j)
            if (bitvec_get(j, unfolded))
              bitvec_set(j % prime, folded);

          // record bit count
          int bitCount = bitvec_count(folded, words);
          bitCounts.push_back(bitCount);

          outputFile.writeFingerprint(folded);
        }

        // free bit vector
        delete [] folded;

        unsigned int average_count = std::accumulate(bitCounts.begin(), bitCounts.end(), 0) / inputFile.numFingerprints();
        unsigned int min_count = *std::min_element(bitCounts.begin(), bitCounts.end());
        unsigned int max_count = *std::max_element(bitCounts.begin(), bitCounts.end());

        // create JSON header
        std::string json = inputFile.header();
        Json::Reader reader;
        Json::Value data;
        reader.parse(json, data);
        data["num_bits"] = bits;
        data["fingerprint"]["prime"] = prime;
        data["statistics"]["average_count"] = average_count;
        data["statistics"]["min_count"] = min_count;
        data["statistics"]["max_count"] = max_count;

        // write JSON header
        Json::StyledWriter writer;
        outputFile.writeHeader(writer.write(data));

        return 0;
      }

  };

  class FoldToolFactory : public HeliumToolFactory
  {
    public:
      HELIUM_TOOL("fold", "Fold fingerprint index files", 3, FoldTool);

      /**
       * Get usage information.
       */
      std::string usage(const std::string &command) const
      {
        std::stringstream ss;
        ss << "Usage: " << command << " <bits> <in_file> <out_file>" << std::endl;
        ss << std::endl;
        ss << "The fold tool can be used to fold fingerprint index files. Any bits argument specifies" << std::endl;
        ss << "the new number of bits, this must be less than the number of bits in the input file." << std::endl;
        ss << std::endl;
        return ss.str();
      }
  };

  FoldToolFactory theFoldToolFactory;

}
