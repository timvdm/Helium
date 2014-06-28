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

#include "args.h"

using namespace Helium;

namespace Helium {

  /**
   * Tool for converting row-major to column-major order fingerprint files.
   */
  class TransposeTool : public HeliumTool
  {
    public:
      /**
       * Perform tool action.
       */
      int run(int argc, char**argv)
      {
        ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("in_file", "out_file"));
        // required arguments
        std::string inFile = args.GetArgString("in_file");
        std::string outFile = args.GetArgString("out_file");

        // open input file
        InMemoryRowMajorFingerprintStorage inputFile;
        if (!inputFile.load(inFile)) {
          std::cerr << inputFile.error().what() << std::endl;
          return -1;
        }

        // open output file
        ColumnMajorFingerprintOutputFile outputFile(outFile, inputFile.numBits(), inputFile.numFingerprints());

        // process fingerprints
        for (unsigned int i = 0; i < inputFile.numFingerprints(); ++i) {
          Word *fingerprint = inputFile.fingerprint(i);
          outputFile.writeFingerprint(fingerprint);
        }

        // create JSON header
        std::string json = inputFile.header();
        Json::Reader reader;
        Json::Value data;
        reader.parse(json, data);
        data["order"] = "column-major";

        // write JSON header
        Json::StyledWriter writer;
        outputFile.writeHeader(writer.write(data));

        return 0;
      }

  };

  class TransposeToolFactory : public HeliumToolFactory
  {
    public:
      HELIUM_TOOL("transpose", "Convert fingerprint files from row-major to column-major order", 2, TransposeTool);

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

  TransposeToolFactory theTransposeToolFactory;


}
