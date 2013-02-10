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
#include "similarity.h"

#include "../src/similarity.h"
#include "../src/fileio.h"
#include "../src/fileio/fingerprints.h"

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

  std::string SimilarityTool::usage(const std::string &command) const
  {
    std::stringstream ss;
    ss << "Usage: " << command << " [options] <query> <fingerprint_file>" << std::endl;
    ss << std::endl;
    ss << "Perform a similarity search on a fingerprint file. The fingerprint file must store the" << std::endl;
    ss << "fingerprints in row-major order. The query has to be a SMILES string." << std::endl;
    ss << std::endl;
    ss << "Options:" << std::endl;
    ss << "    -Tmin <n>     The minimum tanimoto score (default is 0.7)" << std::endl;
    ss << "    -brute        Do brute force search (default is to use index)" << std::endl;
#ifdef HAVE_CPP11
    ss << "    -brute        Do brute force search (default is to use index)" << std::endl;
#endif
    ss << "    -k <n>        When using an index (i.e. no -brute), specify the dimension for the kD-grid (default is 3)" << std::endl;
    ss << std::endl;
    return ss.str();
  }

  int SimilarityTool::run(int argc, char**argv)
  {
    ParseArgs args(argc, argv, ParseArgs::Args("-Tmin(number)", "-brute", "-brute-mt", "-k(number)"), ParseArgs::Args("query", "fingerprint_file"));
    // optional arguments
    const double Tmin = args.IsArg("-Tmin") ? args.GetArgDouble("-Tmin", 0) : 0.7;
    const bool brute = args.IsArg("-brute");
#ifdef HAVE_CPP11
    const bool brute_mt = args.IsArg("-brute-mt");
#endif
    const int k = args.IsArg("-k") ? args.GetArgInt("-k", 0) : 3;
    // required arguments
    std::string smiles = args.GetArgString("query");
    std::string filename = args.GetArgString("fingerprint_file");

    // open fingerprint file
    InMemoryRowMajorFingerprintStorage storage(filename);

    // compute query fingerprint
    HeMol query;
    read_smiles(smiles, query);
    Word *queryFingerprint = compute_fingerprint(storage.header(), query);

    // perform search
    std::vector<std::pair<unsigned int, double> > result;
#ifdef HAVE_CPP11
    if (brute_mt) {
      result = brute_force_similarity_search_threaded(queryFingerprint, storage, Tmin);
    } else
#endif
    if (brute) {
      result = brute_force_similarity_search(queryFingerprint, storage, Tmin);
    } else {
      SimilaritySearchIndex<InMemoryRowMajorFingerprintStorage> index(storage, k);
      result = index.search(queryFingerprint, Tmin);
    }


    std::sort(result.begin(), result.end(), compare_first<unsigned int, double>());

    // print results
    for (std::size_t i = 0; i < result.size(); ++i)
      std::cout << result[i].first << "\t" << result[i].second << std::endl;

    if (!queryFingerprint)
      return 0;

    delete [] queryFingerprint;
    return 0;
  }

}
