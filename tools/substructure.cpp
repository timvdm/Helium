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

#include "../src/similarity.h"
#include "../src/smiles.h"
#include "../src/fileio/fingerprints.h"
#include "../src/fileio/molecules.h"
#include "../src/isomorphism.h"

#include <json/json.h>

#include "args.h"

namespace Helium {

  template<typename MoleculeType>
  Word* compute_fingerprint(const std::string &settings, MoleculeType &mol)
  {
    Json::Reader reader;
    Json::Value data;

    if (!reader.parse(settings, data)) {
      std::cerr << reader.getFormattedErrorMessages() << std::endl;
      return 0;
    }

    int bits = data["num_bits"].asInt();
    int words = bitvec_num_words_for_bits(bits);
    int k = data["fingerprint"]["k"].asInt();
    int prime = data["fingerprint"]["prime"].asInt();
    std::string type = data["fingerprint"]["type"].asString();

    Word *fingerprint = new Word[words];

    if (type == "Helium::paths_fingerprint") {
      path_fingerprint(mol, fingerprint, k, words, prime);
      return fingerprint;
    }

    if (type == "Helium::trees_fingerprint") {
      tree_fingerprint(mol, fingerprint, k, words, prime);
      return fingerprint;
    }

    if (type == "Helium::subgraph_fingerprint") {
      subgraph_fingerprint(mol, fingerprint, k, words, prime);
      return fingerprint;
    }

    std::cerr << "Fingerprint type \"" << type << "\" not recognised" << std::endl;

    delete [] fingerprint;
    return 0;
  }


  void screen(InMemoryColumnMajorFingerprintStorage &storage, Word *fingerprint, Word *result)
  {
    bool first = true;
    for (int i = 0; i < storage.numBits(); ++i) { // foreach bit
      // skip this bit if it is not set
      if (!bitvec_get(i, fingerprint))
        continue;

      if (first) {
        // if this is the first bit, just set result
        for (unsigned int j = 0; j < bitvec_num_words_for_bits(storage.numFingerprints()); ++j)
          result[j] = storage.bit(i)[j];
        first = false;
      } else {
        // do set intersection
        for (unsigned int j = 0; j < bitvec_num_words_for_bits(storage.numFingerprints()); ++j)
          result[j] &= storage.bit(i)[j];
      }
    }
  }


  class SubstructureTool : public HeliumTool
  {
    public:
      /**
       * Perform tool action.
       */
      int run(int argc, char **argv)
      {
        ParseArgs args(argc, argv, ParseArgs::Args("-brute"), 
            ParseArgs::Args("query", "molecule_file", "fingerprint_file"));
        // optional arguments
        const bool brute = args.IsArg("-brute");
        // required arguments
        std::string smiles = args.GetArgString("query");
        std::string moleculeFilename = args.GetArgString("molecule_file");
        std::string fingerprintFilename = args.GetArgString("fingerprint_file");

        // open fingerprint file
        InMemoryColumnMajorFingerprintStorage storage;
        try {
          storage.load(fingerprintFilename);
        } catch (const std::exception &e) {
          std::cerr << e.what() << std::endl;
          return -1;
        }

        // compute query fingerprint
        HeMol query;
        parse_smiles(smiles, query);
        Word *queryFingerprint = compute_fingerprint(storage.header(), query);

        // perform search
        Word *candidates = new Word[bitvec_num_words_for_bits(storage.numFingerprints())];

        screen(storage, queryFingerprint, candidates);

        std::ifstream ifs(moleculeFilename.c_str(), std::ios_base::binary | std::ios_base::in);
        unsigned int numMolecules;
        read32(ifs, numMolecules);

        HeMol mol;
        std::vector<Index> result;
        for (Index i = 0; i < numMolecules; ++i) {
          read_molecule(ifs, mol);
          if (bitvec_get(i, candidates)) {
            if (isomorphism_search<DefaultAtomMatcher, DefaultBondMatcher>(mol, query))
              result.push_back(i);
          }
        }

        // print results
        Json::Value data;
        data["hits"] = Json::Value(Json::arrayValue);
        data["screened"] = Json::Value(bitvec_count(candidates, bitvec_num_words_for_bits(storage.numFingerprints())));
        data["confirmed"] = Json::Value(static_cast<unsigned int>(result.size()));
        data["false_positives"] = Json::Value(1.0 - static_cast<double>(result.size()) / bitvec_count(candidates, bitvec_num_words_for_bits(storage.numFingerprints())));
        for (std::size_t i = 0; i < result.size(); ++i) {
          data["hits"][Json::ArrayIndex(i)] = Json::Value(Json::objectValue);
          Json::Value &obj = data["hits"][Json::ArrayIndex(i)];
          obj["index"] = result[i];
        }

        Json::StyledWriter writer;
        std::cout << writer.write(data);

        // deallocate fingerprint
        if (queryFingerprint)
          delete [] queryFingerprint;
        delete [] candidates;

        return 0;
      }

  };
  
  class SubstructureToolFactory : public HeliumToolFactory
  {
    public:
      HELIUM_TOOL("substructure", "Perform a substructure search", 3, SubstructureTool);

      /**
       * Get usage information.
       */
      std::string usage(const std::string &command) const
      {
        std::stringstream ss;
        ss << "Usage: " << command << " [options] <query> <molecule_file> <fingerprint_file>" << std::endl;
        ss << std::endl;
        ss << "Perform a substructure search. The fingerprint file must store the" << std::endl;
        ss << "fingerprints in column-major order. The query has to be a SMILES string." << std::endl;
        ss << std::endl;
        ss << "Optionally, the <query> can be replaced with 'interactive' to start an interactive session" << std::endl;
        ss << std::endl;
        ss << "Options:" << std::endl;
        ss << "    -brute        Do brute force search (default is to use index)" << std::endl;
        ss << std::endl;
        return ss.str();
      }
  };

  SubstructureToolFactory theSubstructureToolFactory;

}
