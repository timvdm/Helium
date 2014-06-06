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

#include <Helium/diagram.h>
#include <Helium/hemol.h>
#include <Helium/smiles.h>
#include <Helium/fileio/moleculefile.h>
#include <Helium/algorithms/gtd.h>
#include <Helium/algorithms/cycles.h>
#include <Helium/depict/svgpainter.h>
#include <Helium/depict/depict.h>
#include <Helium/algorithms/kekulize.h>

#include <boost/timer/timer.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "args.h"

const boost::timer::nanosecond_type one_milisecond(1000000L);

namespace Helium {

  template<typename MoleculeType, template<typename> class Algorithm>
  int measure_ms(MoleculeType &mol, const Algorithm<MoleculeType> &algorithm)
  {
    boost::timer::cpu_timer timer;

    algorithm(mol);

    boost::timer::cpu_times elapsed = timer.elapsed();
    int ms = (elapsed.system + elapsed.user) / one_milisecond;
    return ms;
  }

  template<typename MoleculeType>
  struct RelevantCycles
  {
    void operator()(const MoleculeType &mol) const
    {
      relevant_cycles(mol);
    }
  };

  template<typename MoleculeType>
  struct Kekulize
  {
    void operator()(MoleculeType &mol) const
    {
      kekulize(mol);
    }
  };

  template<typename MoleculeType>
  struct GenerateDiagram
  {
    void operator()(const MoleculeType &mol) const
    {
      generate_diagram(mol);
    }
  };

  class BenchmarkTool : public HeliumTool
  {
    public:
      /**
       * Perform tool action.
       */
      int run(int argc, char **argv)
      {
        ParseArgs args(argc, argv, ParseArgs::Args(),
            ParseArgs::Args("algorithm", "filename"));
        // required arguments
        std::string algorithm = args.GetArgString("algorithm");
        std::string filename = args.GetArgString("filename");

        HeMol mol;
        MoleculeFile file;
        try {
          file.load(filename);
        } catch (const std::exception &e) {
          std::cerr << e.what() << std::endl;
          return -1;
        }

        int algo = 0;
        if (algorithm == "relevant_cycles")
          algo = 1;
        else if (algorithm == "kekulize")
          algo = 2;
        else if (algorithm == "generate_diagram")
          algo = 3;

        if (!algo) {
          std::cerr << "Unkown algorithm: " << algorithm << std::endl;
          return -1;
        }

        int ms;
        boost::timer::cpu_timer timer;

        for (std::size_t i = 0; i < file.numMolecules(); ++i) {
          file.readMolecule(mol);

          switch (algo) {
            case 1:
              ms = measure_ms(mol, RelevantCycles<HeMol>());
              break;
            case 2:
              ms = measure_ms(mol, Kekulize<HeMol>());
              break;
            case 3:
              ms = measure_ms(mol, GenerateDiagram<HeMol>());
              break;
          }

          std::cout << i << ":" << ms << " ms" << std::endl;
        }

        boost::timer::cpu_times elapsed = timer.elapsed();
        ms = (elapsed.system + elapsed.user) / one_milisecond;
        std::cout << "total: " << ms / 1000 << " s" << std::endl;

        return 0;
      }

  };

  class BenchmarkToolFactory : public HeliumToolFactory
  {
    public:
      HELIUM_TOOL("benchmark", "Benchmark algorithms", 2, BenchmarkTool);

      /**
       * Get usage information.
       */
      std::string usage(const std::string &command) const
      {
        std::stringstream ss;
        ss << "Usage: " << command << " <algorithm> <molecule_file>" << std::endl;
        ss << std::endl;
        ss << "Algorithms:" << std::endl;
        ss << "    relevant_cycles" << std::endl;
        ss << "    kekulize" << std::endl;
        ss << "    generate_diagram" << std::endl;
        ss << std::endl;
        return ss.str();
      }
  };

  BenchmarkToolFactory theBenchmarkToolFactory;

}
