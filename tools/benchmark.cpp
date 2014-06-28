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
#include <Helium/algorithms/extendedconnectivities.h>
#include <Helium/algorithms/canonical.h>
#include <Helium/algorithms/components.h>
#include <Helium/algorithms/cycles.h>
#include <Helium/algorithms/cycles.h>
#include <Helium/depict/svgpainter.h>
#include <Helium/depict/depict.h>
#include <Helium/algorithms/kekulize.h>

#include <boost/timer/timer.hpp>
//#include <boost/date_time/gregorian/gregorian.hpp>
//#include <boost/date_time/posix_time/posix_time.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>


#include "args.h"

const boost::timer::nanosecond_type one_second(1000000000L);
const boost::timer::nanosecond_type one_millisecond(1000000L);
//const boost::timer::nanosecond_type one_microsecond(1000L);

namespace Helium {

  template<typename MoleculeType>
  bool is_planar(const MoleculeType &mol)
  {
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::property<boost::vertex_index_t, int> > graph;

    graph g(num_atoms(mol));
    FOREACH_BOND (bond, mol)
      add_edge(get_index(mol, get_source(mol, *bond)), get_index(mol, get_target(mol, *bond)), g);

    return boost::boyer_myrvold_planarity_test(g);
  }

  template<typename MoleculeType, template<typename> class Algorithm>
  unsigned long measure_ns(MoleculeType &mol, const Algorithm<MoleculeType> &algorithm)
  {
    boost::timer::cpu_timer timer;

    algorithm(mol);

    boost::timer::cpu_times elapsed = timer.elapsed();
    return elapsed.wall;
  }

  template<typename MoleculeType>
  struct CycleMembership
  {
    void operator()(const MoleculeType &mol) const
    {
      std::vector<bool> atoms, bonds;
      cycle_membership(mol, atoms, bonds);
    }
  };

  template<typename MoleculeType>
  struct AtomComponents
  {
    void operator()(const MoleculeType &mol) const
    {
      std::vector<unsigned int> atoms = connected_atom_components(mol);
    }
  };

  template<typename MoleculeType>
  struct BondComponents
  {
    void operator()(const MoleculeType &mol) const
    {
      std::vector<unsigned int> bonds = connected_bond_components(mol);
    }
  };

  template<typename MoleculeType>
  struct ExtendedConnectivities
  {
    void operator()(const MoleculeType &mol) const
    {
      std::vector<unsigned long> ec = extended_connectivities(mol, DefaultAtomInvariant());
    }
  };

  template<typename MoleculeType>
  struct Canonicalize
  {
    void operator()(const MoleculeType &mol) const
    {
      std::vector<unsigned int> atoms = connected_atom_components(mol);
      std::vector<unsigned int> bonds = connected_bond_components(mol);
      std::vector<unsigned long> ec = extended_connectivities(mol, DefaultAtomInvariant());
      canonicalize(mol, ec, DefaultAtomInvariant(), DefaultBondInvariant(), atoms, bonds);
    }
  };

  template<typename MoleculeType>
  struct RelevantCycles
  {
    void operator()(const MoleculeType &mol) const
    {
      relevant_cycles_vismara(mol);
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

  enum Mode {
    Individual = 1, // measure individual molecules
    Graph = 2, // create a graph
    LastMode
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
        if (!file.load(filename)) {
          std::cerr << file.error().what() << std::endl;
          return -1;
        }

        int algo = 0;
        if (algorithm == "read")
          algo = 1;
        else if (algorithm == "cycle_membership")
          algo = 2;
        else if (algorithm == "relevant_cycles")
          algo = 3;
        else if (algorithm == "kekulize")
          algo = 4;
        else if (algorithm == "generate_diagram")
          algo = 5;
        else if (algorithm == "connected_atom_components")
          algo = 6;
        else if (algorithm == "connected_bond_components")
          algo = 7;
        else if (algorithm == "extended_connectivities")
          algo = 8;
        else if (algorithm == "canonicalize")
          algo = 9;

        if (!algo) {
          std::cerr << "Unkown algorithm: " << algorithm << std::endl;
          return -1;
        }

        int mode = Graph;

        unsigned long ns;
        boost::timer::cpu_timer timer;

        for (std::size_t i = 0; i < file.numMolecules(); ++i) {
          file.readMolecule(mol);

          //std::cout << i << std::endl;

          switch (algo) {
            case 1:
              break;
            case 2:
              ns = measure_ns(mol, CycleMembership<HeMol>());
              break;
            case 3:
              ns = measure_ns(mol, RelevantCycles<HeMol>());
              break;
            case 4:
              ns = measure_ns(mol, Kekulize<HeMol>());
              break;
            case 5:
              ns = measure_ns(mol, GenerateDiagram<HeMol>());
              break;
            case 6:
              ns = measure_ns(mol, AtomComponents<HeMol>());
              break;
            case 7:
              ns = measure_ns(mol, BondComponents<HeMol>());
              break;
            case 8:
              ns = measure_ns(mol, ExtendedConnectivities<HeMol>());
              break;
            case 9:
              ns = measure_ns(mol, Canonicalize<HeMol>());
              break;
          }

          switch (mode) {
            case Individual:
              std::cout << i << ":" << ns << " ns" << std::endl;
              break;
            case Graph:
              if ((i % 1000) == 0) {
                boost::timer::cpu_times elapsed = timer.elapsed();
                ns = elapsed.wall;
                std::cout << i << "," << ns / one_millisecond << std::endl;
              }
              break;
          }
        }

        boost::timer::cpu_times elapsed = timer.elapsed();
        ns = elapsed.wall;
        switch (mode) {
          case Individual:
            std::cout << "total: " << ns / one_second << " s" << std::endl;
            break;
          case Graph:
            std::cout << file.numMolecules() << "," << ns / one_millisecond << std::endl;
            break;
        }

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
        ss << "    read" << std::endl;
        ss << "    cycle_membership" << std::endl;
        ss << "    relevant_cycles" << std::endl;
        ss << "    kekulize" << std::endl;
        ss << "    generate_diagram" << std::endl;
        ss << "    atom_components" << std::endl;
        ss << "    bond_components" << std::endl;
        ss << "    extended_connectivities" << std::endl;
        ss << "    canonicalize" << std::endl;
        ss << std::endl;
        return ss.str();
      }
  };

  BenchmarkToolFactory theBenchmarkToolFactory;

}
