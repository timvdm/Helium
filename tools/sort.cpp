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
#include "tool.h"

#include <Helium/molecule.h>
#include <Helium/fileio/moleculefile.h>
#include <Helium/fileio/fingerprints.h>
#include <Helium/fingerprints/fingerprints.h>

#include <numeric> // std::accumulate

#include "args.h"

namespace Helium {

  /**
   * Tool for sroting fingerprint indexes based on population count.
   */
  class SortTool : public HeliumTool
  {
    public:
      /**
       * Perform tool action.
       */
      int run(int argc, char **argv)
      {
        ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("in_file", "out_file"));
        // required arguments
        std::string inFile = args.GetArgString("in_file");
        std::string outFile = args.GetArgString("out_file");



        //
        // open input fingerprint file
        //
        InMemoryRowMajorFingerprintStorage storage;
        if (!storage.load(inFile)) {
          std::cerr << storage.error().what() << std::endl;
          return -1;
        }

        //
        // open output fingerprint file
        //
        RowMajorFingerprintOutputFile indexFile(outFile, storage.numBits());

        //
        // sort the fingerprints
        //
        int numWords = bitvec_num_words_for_bits(storage.numBits());
        for (int c = 0; c < storage.numBits(); ++c) {
          for (unsigned int i = 0; i < storage.numFingerprints(); ++i) {
            int popcnt = bitvec_count(storage.fingerprint(i), numWords);
            if (c == popcnt)
              indexFile.writeFingerprint(storage.fingerprint(i));
          }
        }

        // write JSON header
        Json::StyledWriter writer;
        indexFile.writeHeader(storage.header());

        return 0;
      }

  };

  class SortToolFactory : public HeliumToolFactory
  {
    public:
      HELIUM_TOOL("sort", "Sort fingerprint index files by population count", 2, SortTool);

      /**
       * Get usage information.
       */
      std::string usage(const std::string &command) const
      {
        std::stringstream ss;
        ss << "Usage: " << command << " <in_file> <out_file>" << std::endl;
        ss << std::endl;
        return ss.str();
      }
  };

  SortToolFactory theSortToolFactory;

}
