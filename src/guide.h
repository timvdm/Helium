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

namespace Helium {

/**
 * @page user_guide Helium User Guide
 *
 * @section compiling Building Helium
 *
 * @subsection compiling_deps Dependencies
 *
 * Helium depends on the following libraries:
 *
 * @li Boost C++ Libraries (1.48 or later): http://www.boost.org/
 * @li @b optional: OpenBabel 2.x: http://www.openbabel.org/
 *
 * Although OpenBabel is optional, this is highly recommended since it is used
 * for file conversion.
 *
 * CMake (http://www.cmake.org/) is required for building.
 *
 * @subsection compiling_make Compiling
 *
 * Once the Helium source code has been obtained, navigate to the top-level
 * directory. There, the following sequence of commands will configure, build,
 * test and install Helium:
 *
 @verbatim
 ../Helium $ mkdir build
 ../Helium $ cd build
 ../Helium/build $ cmake ..
 ../Helium/build $ make
 ../Helium/build $ make test
 ../Helium/build $ sudo make install
 @endverbatim
 *
 * @subsection compiling_use Using Helium in Your Own Projects
 *
 * In the examples/external directory, there is an example showing how external
 * programs can be build that depend on Helium. The 2 top lines in the
 * CMakeLists.txt file are all that need to be changed for simple applications.
 *
 * @section user_guide_files Creating Helium Molecule Files
 *
 * Helium uses binary molecule files for better performance. These can be created
 * with the obhelium tool that comes with Helium if the OpenBabel libraries were
 * found. The usage is very simple:
 *
 @verbatim
 $ obhelium
 $ Usage: obhelium <in_file> <out_file>
 @endverbatim
 *
 * For example, the SMILES file molecules.smi is converted to molecules.hel using
 * this command:
 *
 @verbatim
 $ obhelium molecules.smi molecules.hel
 @endverbatim
 *
 * @section user_guide_helium The Helium Tool
 *
 * The helium tool performs various tasks depending on it's first argument. Running
 * helium without any parameters lists the options:
 *
 @verbatim
 $ helium
 Usage: helium <tool> [options]

 Tools:
     header         Extract the JSON header from binary Helium files
     index          Create fingerprint index files
     fold           Fold fingerprint index files
     transpose      Convert fingerprint files from row-major to column-major order
     similarity     Perform a similarity search on a fingerprint index file
     similarityNxN  Perform an NxN similarity search on a fingerprint index file
     sort           Sort fingerprint index files by population count
     substructure   Perform a substructure search
     smarts         Perform a SMARTS search
 @endverbatim
 *
 * More information on a specific task can is given when only the first
 * argument is given:
 *
 @verbatim
 $ helium smarts
 Usage: helium smarts [options] <SMARTS> <molecule_file>

 Perform a SMARTS search.

 Options:
     -styled       Output nicely formatted JSON (default is fast non-human friendly JSON)
     -smiles       Output SMILES
 @endverbatim
 *
 * The @b header (sub)tool can be used to print out the JSON header that all
 * Helium binary files (e.g. molecule file, fingerprints file, ...) contain:
 *
 @verbatim
 $ helium header molecules.hel
 {
      "filetype" : "molecules",
      "molecule_indexes" : 289823869,
      "num_molecules" : 1000000
 }
 @endverbatim
 *
 * The @b smarts (sub)tool can be used to perform direct (i.e. without any index)
 * SMARTS queries on molecule files:
 *
 @verbatim
 $ smarts -smiles -styled "c1ccccc1[N+](=O)[O-]" ../data/1K.hel
 {
      "hits" : [
              {
                    "index" : 29,
                    "smiles" : "[O-][N+](=O)c1c(N)ccc(Br)c1F"
              },
              {
                    "index" : 126,
                    "smiles" : "Brc1ccc(CC(=O)C(=O)O)c(c1)[N+](=O)[O-]"
              },
              {
                "index" : 170,
                "smiles" : "COc1ccc(cc1)Nc2cc(Cl)c(cc2C(=O)O)[N+](=O)[O-]"
              },
              ...
        ],
   "num_hits" : 18
 }
 @endverbatim
 *
 *
 * The ramining tools currently deal with fingerprints.
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */

}
