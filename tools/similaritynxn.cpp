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

#include <numeric>
#include <algorithm>
#include <functional>

#include <Helium/fingerprints/similarity.h>
#include <Helium/fileio/fingerprints.h>
#include <Helium/fileio/fps.h>
#include <Helium/smiles.h>
#include <Helium/json/json.h>

#ifdef HAVE_CPP11
#include <Helium/concurrent.h>
#endif

#ifdef HAVE_OPENCL
#include <CL/cl.hpp>
#endif


#include "args.h"

namespace Helium {

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
 
#ifdef HAVE_OPENCL 
  // defined in tools/opencl.cpp
  void checkError(cl_int ret, const std::string &msg);

  void contextCallback(const char *error_info, const void *private_info, size_t cb, void *user_data)
  {
    std::cerr << "Error: " << error_info << std::endl;
    std::exit(-1);
  }

  /**
   * @brief Read a files contents into a std::string.
   *
   * This function is used to get the source for an OpenCL program.
   *
   * @param filename The file containing the OpenCL source code.
   *
   * @return std::string containing the files contents (source code).
   */
  std::string read_source(const std::string &filename)
  {
    std::ifstream ifs(filename.c_str());
    if (!ifs) {
      std::cerr << "Error: Could not open file " << filename << " for reading" << std::endl;
      return std::string();
    }

    std::string line, source;
    while (std::getline(ifs, line))
      source += line + "\n";

    return source;
  }
#endif


  // alternative method for making similarity search threaded
  template<typename FingerprintStorageType>
  void run_similarity_search(const SimilaritySearchIndex<FingerprintStorageType> &index, double Tmin,
      const std::vector<Word*> &queries, std::vector<std::vector<std::pair<unsigned int, double> > > &result,
      unsigned int begin, unsigned int end)
  {
    for (unsigned int i = begin; i < end; ++i)
      result[i] = index.search(queries[i], Tmin);
  }

  template<typename T>
  T padded(T value, std::size_t unit)
  {
    return value + unit - (value % unit);
  }

#ifdef HAVE_OPENCL
  struct OpenCLResult
  {
    cl_uint index;
    cl_float T;
  };
#endif

  class SimilarityNxNTool : public HeliumTool
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
        ParseArgs args(argc, argv, ParseArgs::Args("-Tmin(number)", "-brute", "-benchmark",
#ifdef HAVE_CPP11
              "-brute-mt", "-mt",
#endif
#ifdef HAVE_OPENCL
              "-opencl", "-platform(number)", "-device(number)",
#endif
              "-k(number)", "-N(number)"), ParseArgs::Args("fingerprint_file"));
        // optional arguments
        const double Tmin = args.IsArg("-Tmin") ? args.GetArgDouble("-Tmin", 0) - 10e-5 : 0.7 - 10e-5;
        bool brute = args.IsArg("-brute");
#ifdef HAVE_CPP11
        const bool brute_mt = args.IsArg("-brute-mt");
        const bool mt = args.IsArg("-mt");
#endif
#ifdef HAVE_OPENCL
        const bool opencl = args.IsArg("-opencl");
        const int platform_id = args.IsArg("-platform") ? args.GetArgInt("-platform", 0) : 1;
        const int device_id = args.IsArg("-device") ? args.GetArgInt("-device", 0) : 1;
#endif
        const bool benchmark = args.IsArg("-benchmark");
        const int k = args.IsArg("-k") ? args.GetArgInt("-k", 0) : 3;
        const int N = args.IsArg("-N") ? args.GetArgInt("-N", 0) : 10;
        // required arguments
        std::string filename = args.GetArgString("fingerprint_file");

        //
        // check for incompatible arguments
        //
        if (brute && args.IsArg("-k"))
          std::cerr << "Option -k <n> has no effect when using option -brute, -k will be ignored." << std::endl;
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
#ifdef HAVE_OPENCL
        if (opencl && args.IsArg("-k"))
          std::cerr << "Option -k <n> has no effect when using option -opencl, -k will be ignored." << std::endl;
        if (opencl && brute)
          std::cerr << "Option -brute has no effect when using option -opencl, -brute will be ignored." << std::endl;
#endif
#if defined(HAVE_CPP11) && defined(HAVE_OPENCL)
        if (opencl && brute_mt)
          std::cerr << "Option -brute-mt has no effect when using option -opencl, -brute-mt will be ignored." << std::endl;
        if (opencl && mt)
          std::cerr << "Option -mt has no effect when using option -opencl, -mt will be ignored." << std::endl;
#endif


        //
        // open fingerprint file
        //
        InMemoryRowMajorFingerprintStorage storage;
        if (!storage.load(filename)) {
          std::cerr << storage.error().what() << std::endl;
          return -1;
        }

        //
        // perform search
        //
        std::vector<std::vector<std::pair<unsigned int, double> > > result(storage.numFingerprints());
#ifdef HAVE_CPP11
        if (brute_mt) {
          std::vector<std::pair<unsigned int, double> > tmp;
          for (std::size_t i = 0; i < storage.numFingerprints(); ++i) {
            tmp = brute_force_similarity_search_threaded(storage.fingerprint(i), storage, Tmin);
            std::sort(tmp.begin(), tmp.end(), compare_second<unsigned int, double, std::greater>());
            tmp.resize(std::min(N, static_cast<int>(tmp.size())));
            result[i] = tmp;
          }
        } else
#endif
        if (brute) {
          std::vector<std::pair<unsigned int, double> > tmp;
          for (std::size_t i = 0; i < storage.numFingerprints(); ++i) {
            tmp = brute_force_similarity_search(storage.fingerprint(i), storage, Tmin);
            std::sort(tmp.begin(), tmp.end(), compare_second<unsigned int, double, std::greater>());
            tmp.resize(std::min<int>(N, static_cast<int>(tmp.size())));
            result[i] = tmp;
          }
        } else
#ifdef HAVE_OPENCL
        if (opencl) {

          // get a list of OpenCL platforms
          std::vector<cl::Platform> platforms;
          checkError(cl::Platform::get(&platforms), "Could not get OpenCL platforms");

          // check if platform is valid
          if (platform_id < 1 || platform_id > platforms.size()) {
            std::cerr << "Error: invalid OpenCL platform id (" << platform_id << ")" << std::endl;
            return -1;
          }

          cl::Platform &platform = platforms[platform_id - 1];
          // get a list of devices for this platform
          std::vector<cl::Device> devices;
          checkError(platform.getDevices(CL_DEVICE_TYPE_ALL, &devices), "Could not get OpenCL platform devices");
        
          // check if device is valid
          if (device_id < 1 || device_id > devices.size()) {
            std::cerr << "Error: invalid OpenCL device id (" << device_id << ")" << std::endl;
            return -1;
          }

          std::vector<cl::Device> device(1, devices[device_id - 1]);
          
          // create OpenCL context
          cl_int err;
          cl::Context context(device, // devices
                              NULL, // properties
                              contextCallback, // callback
                              NULL, // uer data
                              &err);
          checkError(err, "Could not create OpenCL context");

          // create OpenCL program
          cl::Program program(context, read_source("similarity.cl").c_str(), err);
          checkError(err, "Could not create OpenCL program");

          // build the program
          err = program.build(device, //devices
                              //"-cl-std=CL1.2", // options
                              NULL, // options
                              NULL, // callback
                              NULL); // user data

          // get build log
          if (err != CL_SUCCESS) {
            for (std::size_t i = 0; i < device.size(); ++i) {
              std::string build_log;
              checkError(program.getBuildInfo(device[i], CL_PROGRAM_BUILD_LOG, &build_log), "Could not get OpenCL build log");
              if (build_log.size()) {
                std::cerr << "Device " << i + 1 << " build log:" << std::endl;
                std::cerr << build_log << std::endl;
              }
            }
          }

          checkError(err, "Could not build OpenCL program");

          // create OpenCL kernels
          std::vector<cl::Kernel> kernels;
          checkError(program.createKernels(&kernels), "Could not create OpenCL kernels");

          // create OpenCL command queue
          cl::CommandQueue queue(context, device[0], 0, &err);
          checkError(err, "Could not create OpenCL command queue");

          // create OpenCL buffer to hold fingerprints
          std::size_t fingerprints_buffer_size = storage.numFingerprints() * bitvec_num_words_for_bits(storage.numBits()) * sizeof(Word);
          cl::Buffer fingerprints_buffer(context, // context
                                         CL_MEM_READ_ONLY, // flags
                                         fingerprints_buffer_size, // size
                                         NULL, // host pointer
                                         &err);
          checkError(err, "Could not create OpenCL buffer to hold fingerprints");

          // create OpenCL buffer to hold fingerprints
          cl::Buffer results_buffer(context, // context
                                    //CL_MEM_WRITE_ONLY, // flags
                                    CL_MEM_READ_WRITE, // flags
                                    N * storage.numFingerprints() * sizeof(OpenCLResult), // size
                                    NULL, // host pointer
                                    &err);
          checkError(err, "Could not create OpenCL buffer to hold results");

          // copy fingerprints to OpenCL device
          checkError(queue.enqueueWriteBuffer(fingerprints_buffer, // buffer
                                              CL_TRUE, // blocking write
                                              0, // offset
                                              fingerprints_buffer_size, // size
                                              storage.fingerprint(0), // host pointer
                                              NULL, // events
                                              NULL), // event
                     "Could not copy fingerprints to OpenCL device");
          
          // set OpenCL kernel arguments
          cl::Kernel &kernel = kernels[0];
          checkError(kernel.setArg(0, (cl_uint)storage.numFingerprints()), "Could not set num_fingerprints kernel argument (arg 0)");
          checkError(kernel.setArg(1, (cl_uint)bitvec_num_words_for_bits(storage.numBits())), "Could not set num_fingerprints kernel argument (arg 1)");
          checkError(kernel.setArg(2, (cl_float)Tmin), "Could not set threshold kernel argument (arg 2)");
          checkError(kernel.setArg(3, (cl_uint)N), "Could not set N kernel argument (arg 3)");
          checkError(kernel.setArg(4, fingerprints_buffer), "Could not set fingerprints kernel argument (arg 4)");
          checkError(kernel.setArg(5, results_buffer), "Could not set results kernel argument (arg 5)");

          size_t work_group_size;
          checkError(kernel.getWorkGroupInfo(device[0], CL_KERNEL_WORK_GROUP_SIZE, &work_group_size), "Could not get OpenCL kernel work group size");
          work_group_size = std::min(work_group_size, static_cast<size_t>(64));

          unsigned int num_fingerprints = padded(storage.numFingerprints(), work_group_size);

          // enueue the OpenCL kernel to run
          checkError(queue.enqueueNDRangeKernel(kernel, // kernel
                                                cl::NDRange(0), // offset
                                                cl::NDRange(num_fingerprints), // global
                                                cl::NDRange(work_group_size), // local
                                                NULL, // events
                                                NULL), // event
                     "Could not enqueue ND-range kernel to run");

          checkError(queue.finish(), "Could not finish queue");

          // copy results back to host
          std::vector<OpenCLResult> results(N * storage.numFingerprints());
          checkError(queue.enqueueReadBuffer(results_buffer, // buffer
                                             CL_TRUE, // blocking read,
                                             0, // offset
                                             results.size() * sizeof(OpenCLResult), // size
                                             &results[0], // host pointer,
                                             NULL, // events
                                             NULL), // event
                     "Could not copy results from OpenCL device");
    
          // convert results to std::vector result
          for (std::size_t i = 0; i < storage.numFingerprints(); ++i) {
            for (std::size_t j = 0; j < N; ++j) {
              const OpenCLResult &hit = results[i * N + j];
              if (hit.T < 0.0)
                break;
              if (hit.T >= Tmin)
                result[i].push_back(std::make_pair(hit.index, hit.T));
            }
            std::sort(result[i].begin(), result[i].end(), compare_second<unsigned int, double, std::greater>());
          }











        } else
#endif
        {
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

            std::vector<Word*> queries(storage.numFingerprints());
            for (std::size_t i = 0; i < storage.numFingerprints(); ++i)
              queries[i] = storage.fingerprint(i);

            Concurrent<const CallableType&, TaskType, ResultType> concurrent;
            concurrent.run(CallableType(index, Tmin), queries, result);
          } else {
            // run all searches sequentially using a single thread
            std::vector<std::pair<unsigned int, double> > tmp;
            for (std::size_t i = 0; i < storage.numFingerprints(); ++i) {
              tmp = index.search(storage.fingerprint(i), Tmin);
              std::sort(tmp.begin(), tmp.end(), compare_second<unsigned int, double, std::greater>());
              tmp.resize(std::min(N, static_cast<int>(tmp.size())));
              result[i] = tmp;
            }
          }
#else
          // run all searches sequentially using a single thread
          std::vector<std::pair<unsigned int, double> > tmp;
          for (std::size_t i = 0; i < storage.numFingerprints(); ++i) {
            tmp = index.search(storage.fingerprint(i), Tmin);
            std::sort(tmp.begin(), tmp.end(), compare_second<unsigned int, double, std::greater>());
            tmp.resize(std::min<int>(N, static_cast<int>(tmp.size())));
            result[i] = tmp;
          }
#endif
        }

        if (benchmark)
          return 0;

        // sort the results
        for (std::size_t i = 0; i < storage.numFingerprints(); ++i)
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

  class SimilarityNxNToolFactory : public HeliumToolFactory
  {
    public:
      HELIUM_TOOL("similarityNxN", "Perform an NxN similarity search on a fingerprint index file", 1, SimilarityNxNTool);

      /**
       * Get usage information.
       */
      std::string usage(const std::string &command) const
      {
        std::stringstream ss;
        ss << "Usage: " << command << " [options] <fingerprint_file>" << std::endl;
        ss << std::endl;
        ss << "Perform a similarity search on a fingerprint file. The fingerprint file must store the" << std::endl;
        ss << "fingerprints in row-major order. The query has to be a SMILES string." << std::endl;
        ss << std::endl;
        ss << "Options:" << std::endl;
        ss << "    -Tmin <n>     The minimum tanimoto score (default is 0.7)" << std::endl;
        ss << "    -N <n>        The number of nearest matches to find for each fingerint (default is 10)" << std::endl;
        ss << "    -brute        Do brute force search (default is to use index)" << std::endl;
#ifdef HAVE_CPP11
        ss << "    -brute-mt     Do threaded brute force search (default is to use index)" << std::endl;
        ss << "    -mt           Do threaded index search (default is not to use threads)" << std::endl;
#endif
#ifdef HAVE_OPENCL
        ss << "    -opencl       Use OpenCL (default is not to use OpenCL)" << std::endl;
        ss << "    -platform <n> The OpenCL platform to use (default is to use platform 1)" << std::endl;
#endif
        ss << "    -k <n>        When using an index (i.e. no -brute), specify the dimension for the kD-grid (default is 3)" << std::endl;
        ss << std::endl;
        return ss.str();
      }
  };

  SimilarityNxNToolFactory theSimilarityNxNToolFactory;

}
