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

#include <Helium/smarts.h>
#include <Helium/algorithms/cycles.h>
#include <Helium/smiles.h>
#include <Helium/fileio/moleculefile.h>
#include <Helium/json/json.h>

#include "args.h"

namespace Helium {

  class SmartsTool : public HeliumTool
  {
    public:
      /**
       * Perform tool action.
       */
      int run(int argc, char **argv)
      {
        ParseArgs args(argc, argv, ParseArgs::Args("-styled", "-smiles"),
            ParseArgs::Args("smarts", "molecule_file"));
        // optional arguments
        const bool styled = args.IsArg("-styled");
        const bool writeSmiles = args.IsArg("-smiles");
        // required arguments
        std::string smarts = args.GetArgString("smarts");
        std::string moleculeFilename = args.GetArgString("molecule_file");

        // creacte SMARTS query
        Smarts s;
        if (!s.init(smarts)) {
          std::cerr << "Error: " << s.error().what() << std::endl;
          return -1;
        }

        // perform search

        MemoryMappedMoleculeFile moleculeFile;
        try {
          moleculeFile.load(moleculeFilename);
        } catch (const std::exception &e) {
          std::cerr << e.what() << std::endl;
          return -1;
        }

        HeMol mol;
        std::vector<Index> result;
        std::vector<std::string> smiles;
        for (unsigned int i = 0; i < moleculeFile.numMolecules(); ++i) {
          moleculeFile.readMolecule(i, mol);

          bool found = false;
          if (s.requiresCycles())
            found = s.find(mol, relevant_cycles(mol));
          else
            found = s.find(mol, RingSet<HeMol>(mol));

          if (found) {
            result.push_back(i);
            if (writeSmiles) {
              Smiles SMILES;
              smiles.push_back(SMILES.write(mol));
            }
          }
        }

        // print results
        Json::Value data;
        data["num_hits"] = Json::Int(result.size());
        data["hits"] = Json::Value(Json::arrayValue);
        for (std::size_t i = 0; i < result.size(); ++i) {
          data["hits"][Json::ArrayIndex(i)] = Json::Value(Json::objectValue);
          Json::Value &obj = data["hits"][Json::ArrayIndex(i)];
          obj["index"] = result[i];
          if (writeSmiles)
            obj["smiles"] = smiles[i];
        }

        if (styled) {
          Json::StyledWriter writer;
          std::cout << writer.write(data);
        } else {
          Json::FastWriter writer;
          std::cout << writer.write(data);
        }

        return 0;
      }

  };

  class SmartsToolFactory : public HeliumToolFactory
  {
    public:
      HELIUM_TOOL("smarts", "Perform a SMARTS search", 2, SmartsTool);

      /**
       * Get usage information.
       */
      std::string usage(const std::string &command) const
      {
        std::stringstream ss;
        ss << "Usage: " << command << " [options] <SMARTS> <molecule_file>" << std::endl;
        ss << std::endl;
        ss << "Perform a SMARTS search." << std::endl;
        ss << std::endl;
        ss << "Options:" << std::endl;
        ss << "    -styled       Output nicely formatted JSON (default is fast non-human friendly JSON)" << std::endl;
        ss << "    -smiles       Output SMILES" << std::endl;
        ss << std::endl;
        return ss.str();
      }
  };

  SmartsToolFactory theSmartsToolFactory;

}
