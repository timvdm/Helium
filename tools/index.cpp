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

#include <Helium/molecule.h>
#include <Helium/fileio/molecules.h>
#include <Helium/fileio/fingerprints.h>
#include <Helium/fingerprints/fingerprints.h>

#include <numeric> // std::accumulate

#include "args.h"

namespace Helium {

  /**
   * Tool for  creating fingerprint indexes.
   */
  class IndexTool : public HeliumTool
  {
    public:
      /**
       * Fingerprint generation method.
       */
      enum Method {
        PathsMethod,
        TreesMethod,
        SubgraphsMethod
      };

      /**
       * Perform tool action.
       */
      int run(int argc, char **argv)
      {
        ParseArgs args(argc, argv, ParseArgs::Args("-k(number)", "-bits(number)"), ParseArgs::Args("method", "in_file", "out_file"));
        // optional arguments
        const int k = args.IsArg("-k") ? args.GetArgInt("-k", 0) : 7;
        const int bits = args.IsArg("-bits") ? args.GetArgInt("-bits", 0) : 1024;
        const int words = bits / (8 * sizeof(Word));
        const int prime = previous_prime(bits);
        // required arguments
        std::string methodString = args.GetArgString("method");
        std::string inFile = args.GetArgString("in_file");
        std::string outFile = args.GetArgString("out_file");

        // parse method argument
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

        // print fingerprint settings
        std::cerr << "Fingerprint settings:" << std::endl;
        std::cerr << "    method: " << methodString.substr(1) << std::endl;
        std::cerr << "    k: " << k << std::endl;
        std::cerr << "    bits: " << bits << std::endl;
        std::cerr << "    prime: " << prime << std::endl;

        std::cerr << "Indexing " << inFile << "..." << std::endl;

        // open index file
        RowMajorFingerprintOutputFile indexFile(outFile, bits);

        // open molecule file
        MoleculeFile file(inFile);
        HeMol mol;

        // allocate bit vector
        Word *fingerprint = new Word[words];
        // keep track of bit counts
        std::vector<int> bitCounts;

        // process molecules
        for (unsigned int i = 0; i < file.numMolecules(); ++i) {
          file.read_molecule(mol);
          if ((i % 100) == 0)
            std::cout << i << std::endl;

          // compute the fingerprint
          switch (method) {
            case PathsMethod:
              path_fingerprint(mol, fingerprint, k, words, prime);
              break;
            case TreesMethod:
              tree_fingerprint(mol, fingerprint, k, words, prime);
              break;
            case SubgraphsMethod:
              subgraph_fingerprint(mol, fingerprint, k, words, prime);
              break;
          }

          // record bit count
          int bitCount = bitvec_count(fingerprint, words);
          bitCounts.push_back(bitCount);

          indexFile.writeFingerprint(fingerprint);
        }

        // free bit vector
        delete [] fingerprint;

        unsigned int average_count = std::accumulate(bitCounts.begin(), bitCounts.end(), 0) / file.numMolecules();
        unsigned int min_count = *std::min_element(bitCounts.begin(), bitCounts.end());
        unsigned int max_count = *std::max_element(bitCounts.begin(), bitCounts.end());

        // create JSON header
        Json::Value data;
        data["filetype"] = "fingerprints";
        data["order"] = "row-major";
        data["num_bits"] = bits;
        data["num_fingerprints"] = file.numMolecules();
        data["fingerprint"] = Json::Value(Json::objectValue);
        std::stringstream ss;
        ss << "Helium::" << methodString.substr(1) << "_fingerprint (k = " << k << ", bits = " << bits << ")";
        data["fingerprint"]["name"] = ss.str();
        ss.str("");
        ss << "Helium::" << methodString.substr(1) << "_fingerprint";
        data["fingerprint"]["type"] = ss.str();
        data["fingerprint"]["k"] = k;
        data["fingerprint"]["prime"] = prime;
        data["statistics"] = Json::Value(Json::objectValue);
        data["statistics"]["average_count"] = average_count;
        data["statistics"]["min_count"] = min_count;
        data["statistics"]["max_count"] = max_count;

        // write JSON header
        Json::StyledWriter writer;
        indexFile.writeHeader(writer.write(data));


        return 0;
      }

  };

  class IndexToolFactory : public HeliumToolFactory
  {
    public:
      HELIUM_TOOL("index", "Create fingerprint index files", 3, IndexTool);

      /**
       * Get usage information.
       */
      std::string usage(const std::string &command) const
      {
        std::stringstream ss;
        ss << "Usage: " << command << " [options] <method> <in_file> <out_file>" << std::endl;
        ss << std::endl;
        ss << "The index tool can be used to create fingerprint index files. Any contents in the output" << std::endl;
        ss << "file will be overwritten." << std::endl;
        ss << std::endl;
        ss << "Methods:" << std::endl;
        ss << "    -paths        Create hashed fingerprints from paths" << std::endl;
        ss << "    -trees        Create hashed fingerprints from trees" << std::endl;
        ss << "    -subgraphs    Create hashed fingerprints from subgraphs" << std::endl;
        ss << std::endl;
        ss << "Options:" << std::endl;
        ss << "    -k <n>        The maximum size of the path/tree/subgraph (default is 7)" << std::endl;
        ss << "    -bits <n>     The number of bits in the fingerprint (default is 1024)" << std::endl;
        ss << std::endl;
        return ss.str();
      }
  };

  IndexToolFactory theIndexToolFactory;

}
