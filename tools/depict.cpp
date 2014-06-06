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
#include <Helium/algorithms/cycles.h>
#include <Helium/algorithms/kekulize.h>
#include <Helium/depict/svgpainter.h>
#include <Helium/depict/depict.h>
#include <boost/filesystem.hpp>
#include <fstream>

#include "args.h"

namespace Helium {

  bool depict(HeMol &mol, const std::string &filename)
  {
    std::ofstream ofs(filename.c_str());
    if (!ofs)
      return false;

    kekulize(mol);
    std::vector<std::pair<double, double> > coords = generate_diagram(mol);

    SVGPainter painter(ofs);
    Depict depict(&painter);
    depict.setFontFamily("Arial");
    depict.setPenWidth(3);
    depict.setOption(Depict::asymmetricDoubleBond);
    depict.drawMolecule(mol, relevant_cycles(mol), coords);

    return true;
  }

  std::string separator()
  {
    boost::filesystem::path slash("/");
    return slash.make_preferred().native();
  }

  class DepictTool : public HeliumTool
  {
    public:
      /**
       * Perform tool action.
       */
      int run(int argc, char **argv)
      {
        ParseArgs args(argc, argv, ParseArgs::Args(),
            ParseArgs::Args("input", "output"));
        // required arguments
        std::string input = args.GetArgString("input");
        std::string output = args.GetArgString("output");

        if (!boost::filesystem::is_regular_file(input)) {
          // try to parse SMILES
          HeMol mol;
          Smiles SMILES;
          if (!SMILES.read(input, mol)) {
            std::cerr << SMILES.error().what() << std::endl;
            return -1;
          }

          if (!depict(mol, output)) {
            std::cerr << "Could not open " << output << " for writing." << std::endl;
            return -1;
          }

          return 0;
        }

        HeMol mol;
        MoleculeFile file;
        try {
          file.load(input);
        } catch (const std::exception &e) {
          std::cerr << e.what() << std::endl;
          return -1;
        }

        if (file.numMolecules() == 1) {
          file.readMolecule(mol);

          if (!depict(mol, output)) {
            std::cerr << "Could not open " << output << " for writing." << std::endl;
            return -1;
          }

          return 0;
        }

        if (!boost::filesystem::is_directory(output)) {
          std::cerr << "<output> must be a directory for multi-molecule files." << std::endl;
          return -1;
        }

        for (std::size_t i = 0; i < file.numMolecules(); ++i) {
          std::cerr << i << std::endl;
          file.readMolecule(mol);

          std::stringstream filename;
          filename << output << separator() << i + 1 << ".svg";
          depict(mol, filename.str());
        }


        return 0;
      }

  };

  class DepictToolFactory : public HeliumToolFactory
  {
    public:
      HELIUM_TOOL("depict", "Depict molecule(s)", 2, DepictTool);

      /**
       * Get usage information.
       */
      std::string usage(const std::string &command) const
      {
        std::stringstream ss;
        ss << "Usage: " << command << " [options] <input> <output>" << std::endl;
        ss << std::endl;
        ss << "<input> can be a SMILES or file and depict will determine what to do based on this:" << std::endl;
        ss << std::endl;
        ss << "- If <input> is a file containing a single molecule, <output> should be the .svg filename." << std::endl;
        ss << "- If <input> is a multi-molecule file, <output> should be a directory and files will" << std::endl;
        ss << "  be named 1.svg 2.svg etc." << std::endl;
        ss << "- If <input> is a SMILES string, <output> should be the .svg filename" << std::endl;
        ss << std::endl;
        ss << "Options:" << std::endl;
        ss << "    -subdir <n>   Create subdirectories containing n depictions at most" << std::endl;
        ss << std::endl;
        return ss.str();
      }
  };

  DepictToolFactory theDepictToolFactory;

}
