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

#include <Helium/fingerprints/similarity.h>
#include <Helium/fileio/fingerprints.h>
#include <Helium/fileio/fps.h>
#include <Helium/smiles.h>

#ifdef HAVE_CPP11
#include <Helium/concurrent.h>
#endif

#include <json/json.h>

#include "args.h"

namespace Helium {

  Word* compute_fingerprint(const std::string &settings, const std::string &smiles)
  {
    Json::Reader reader;
    Json::Value data;

    if (!reader.parse(settings, data)) {
      std::cerr << reader.getFormattedErrorMessages() << std::endl;
      return 0;
    }

    HeMol mol;
    parse_smiles(smiles, mol);

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


  bool is_fps_file(const std::string &filename)
  {
    std::ifstream ifs(filename.c_str());
    if (!ifs)
      return false;
    std::string line;
    std::getline(ifs, line);
    return line == "#FPS1";
  }

  bool load_queries(const std::string &query, const std::string &storage_header, std::vector<Word*> &queries)
  {
    if (is_fps_file(query)) {
      // load fps file
      FpsFile fpsFile;
      fpsFile.load(query);

      // copy fingerprints
      int numWords = bitvec_num_words_for_bits(fpsFile.numBits());
      for (unsigned int i = 0; i < fpsFile.numFingerprints(); ++i)
        queries.push_back(bitvec_copy(fpsFile.fingerprint(i), numWords));
    } else {
      // query is SMILES
      Word *fingerprint = compute_fingerprint(storage_header, query);
      if (!fingerprint)
        return false;
      queries.push_back(fingerprint);
    }

    return queries.size();
  }

  void free_queries(std::vector<Word*> &queries)
  {
    for (std::size_t i = 0; i < queries.size(); ++i)
      delete [] queries[i];
  }

  // functor used by Concurrent
  template<typename FingerprintStorageType>
  struct RunSimilaritySearch
  {
    RunSimilaritySearch(const SimilaritySearchIndex<FingerprintStorageType> &index_, double Tmin_)
      : index(index_), Tmin(Tmin_)
    {
    }

    void operator()(const Word *query, std::vector<std::pair<unsigned int, double> > &result) const
    {
      result = index.search(query, Tmin);
    }

    const SimilaritySearchIndex<FingerprintStorageType> &index;
    const double Tmin;
  };

  // alternative method for making similarity search threaded
  template<typename FingerprintStorageType>
  void run_similarity_search(const SimilaritySearchIndex<FingerprintStorageType> &index, double Tmin,
      const std::vector<Word*> &queries, std::vector<std::vector<std::pair<unsigned int, double> > > &result,
      unsigned int begin, unsigned int end)
  {
    for (unsigned int i = begin; i < end; ++i)
      result[i] = index.search(queries[i], Tmin);
  }

  class SimilarityTool : public HeliumTool
  {
    public:
      /**
       * Perform tool action.
       */
      int run(int argc, char **argv)
      {
        //
        // Argument handling
        //
        ParseArgs args(argc, argv, ParseArgs::Args("-Tmin(number)", "-brute",
#ifdef HAVE_CPP11
              "-brute-mt", "-mt",
#endif
              "-k(number)"), ParseArgs::Args("query", "fingerprint_file"));
        // optional arguments
        const double Tmin = args.IsArg("-Tmin") ? args.GetArgDouble("-Tmin", 0) - 10e-5 : 0.7 - 10e-5;
        bool brute = args.IsArg("-brute");
#ifdef HAVE_CPP11
        const bool brute_mt = args.IsArg("-brute-mt");
        const bool mt = args.IsArg("-mt");
#endif
        const int k = args.IsArg("-k") ? args.GetArgInt("-k", 0) : 3;
        // required arguments
        std::string query = args.GetArgString("query");
        std::string filename = args.GetArgString("fingerprint_file");

        //
        // check for incompatible arguments
        //
#ifdef HAVE_CPP11
        if (brute && brute_mt) {
          std::cerr << "Options -brute and -brute-mt can not be used simultaneously, -brute will be ignored." << std::endl;
          brute = false;
        }
        if (brute_mt && args.IsArg("-k"))
          std::cerr << "Option -k <n> has no effect when using option -brute-mt, -k will be ignored." << std::endl;
        if (brute && mt)
          std::cerr << "Option -mt has no effect when using option -brute, -mt will be ignored." << std::endl;
        if (brute_mt && mt)
          std::cerr << "Option -mt has no effect when using option -brute-mt, -mt will be ignored." << std::endl;
#endif
        if (brute && args.IsArg("-k"))
          std::cerr << "Option -k <n> has no effect when using option -brute, -k will be ignored." << std::endl;


        //
        // open fingerprint file
        //
        InMemoryRowMajorFingerprintStorage storage;
        try {
          storage.load(filename);
        } catch (const std::exception &e) {
          std::cerr << e.what() << std::endl;
          return -1;
        }

        //
        // load queries
        //
        std::vector<Word*> queries;
        if (!load_queries(query, storage.header(), queries)) {
          std::cerr << "Could not load queries" << std::endl;
          return -1;
        }

        //
        // perform search
        //
        std::vector<std::vector<std::pair<unsigned int, double> > > result(queries.size());
#ifdef HAVE_CPP11
        if (brute_mt) {
          for (std::size_t i = 0; i < queries.size(); ++i)
            result[i] = brute_force_similarity_search_threaded(queries[i], storage, Tmin);
        } else
#endif
        if (brute) {
          for (std::size_t i = 0; i < queries.size(); ++i)
            result[i] = brute_force_similarity_search(queries[i], storage, Tmin);
        } else {
          SimilaritySearchIndex<InMemoryRowMajorFingerprintStorage> index(storage, k);
#ifdef HAVE_CPP11
          if (mt) {
            /*
            // run searches in concurrently using multiple threads
            unsigned numThreads = std::thread::hardware_concurrency();
            // c++ implementations may return 0
            if (!numThreads)
              numThreads = 2;

            unsigned int taskSize = queries.size() / numThreads;

            std::vector<std::thread> threads;
            for (int i = 0; i < numThreads; ++i) {
              unsigned int begin = i * taskSize;
              unsigned int end = std::min(static_cast<unsigned int>(queries.size()), (i + 1) * taskSize);
              std::cout << "(" << begin << ", " << end << ")" << std::endl;
              threads.push_back(std::thread(run_similarity_search<InMemoryRowMajorFingerprintStorage>,
                    std::ref(index), Tmin, std::ref(queries), std::ref(result), begin, end));
            }

            for (auto &thread : threads)
              thread.join();
            */

            // run searches in concurrently using multiple threads
            typedef RunSimilaritySearch<InMemoryRowMajorFingerprintStorage> CallableType;
            typedef Word* TaskType;
            typedef std::vector<std::pair<unsigned int, double> > ResultType;

            Concurrent<const CallableType&, TaskType, ResultType> concurrent;
            concurrent.run(CallableType(index, Tmin), queries, result);
          } else {
            // run all searches sequentially using a single thread
            for (std::size_t i = 0; i < queries.size(); ++i)
              result[i] = index.search(queries[i], Tmin);
          }
#else
          // run all searches sequentially using a single thread
          for (std::size_t i = 0; i < queries.size(); ++i)
            result[i] = index.search(queries[i], Tmin);
#endif
        }

        // deallocate fingerprint
        free_queries(queries);

        // sort the results
        for (std::size_t i = 0; i < queries.size(); ++i)
          std::sort(result[i].begin(), result[i].end(), compare_first<unsigned int, double>());


        //
        // print results
        //
        Json::Value data;
        data["hits"] = Json::Value(Json::arrayValue);
        for (std::size_t i = 0; i < result.size(); ++i) {
          data["hits"][Json::ArrayIndex(i)] = Json::Value(Json::arrayValue);
          for (std::size_t j = 0; j < result[i].size(); ++j) {
            data["hits"][Json::ArrayIndex(i)][Json::ArrayIndex(j)] = Json::Value(Json::objectValue);
            Json::Value &obj = data["hits"][Json::ArrayIndex(i)][Json::ArrayIndex(j)];
            obj["index"] = result[i][j].first;
            obj["tanimoto"] = result[i][j].second;
          }
        }

        Json::StyledWriter writer;
        std::cout << writer.write(data);

        return 0;
      }

  };

  class SimilarityToolFactory : public HeliumToolFactory
  {
    public:
      HELIUM_TOOL("similarity", "Perform a similarity search on a fingerprint index file", 2, SimilarityTool);

      /**
       * Get usage information.
       */
      std::string usage(const std::string &command) const
      {
        std::stringstream ss;
        ss << "Usage: " << command << " [options] <query> <fingerprint_file>" << std::endl;
        ss << std::endl;
        ss << "Perform a similarity search on a fingerprint file. The fingerprint file must store the" << std::endl;
        ss << "fingerprints in row-major order. The query has to be a SMILES string." << std::endl;
        ss << std::endl;
        ss << "Optionally, the <query> can be replaced with 'interactive' to start an interactive session" << std::endl;
        ss << std::endl;
        ss << "Options:" << std::endl;
        ss << "    -Tmin <n>     The minimum tanimoto score (default is 0.7)" << std::endl;
        ss << "    -brute        Do brute force search (default is to use index)" << std::endl;
#ifdef HAVE_CPP11
        ss << "    -brute-mt     Do threaded brute force search (default is to use index)" << std::endl;
        ss << "    -mt           Do threaded index search (default is not to use threads)" << std::endl;
#endif
        ss << "    -k <n>        When using an index (i.e. no -brute), specify the dimension for the kD-grid (default is 3)" << std::endl;
        ss << std::endl;
        return ss.str();
      }
  };

  SimilarityToolFactory theSimilarityToolFactory;

}
