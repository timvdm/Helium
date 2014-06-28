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
#include <Helium/molecule.h>
#include <Helium/fileio/file.h>
#include <Helium/bitvec.h>
#include <Helium/fileio/fingerprints.h>
#include <Helium/json/json.h>

#include "args.h"
#include "progress.h"

#include <numeric>

using namespace Helium;

bool is_hex(char digit)
{
  switch (digit) {
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
    case 'a':
    case 'A':
    case 'b':
    case 'B':
    case 'c':
    case 'C':
    case 'd':
    case 'D':
    case 'e':
    case 'E':
    case 'f':
    case 'F':
      return true;
    default:
      return false;
  }
}

int main(int argc, char**argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <in_file> <out_file>" << std::endl;
    return 0;
  }

  ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("in_file", "out_file"));
  std::string inFile = args.GetArgString("in_file");
  std::string outFile = args.GetArgString("out_file");


  // open the input file
  std::ifstream ifs(inFile.c_str());

  int bits = -1;
  std::string type = "fps";

  std::string line;
  // read fps header
  while (ifs && !ifs.eof() && ifs.peek() == '#') {
    std::getline(ifs, line);
    if (line.substr(0, 10) == "#num_bits=") {
      std::stringstream ss(line.substr(10));
      ss >> bits;
    } else if (line.substr(0, 6) == "#type=") {
      type = line.substr(6);
    }
  }

  if (bits == -1) {
    std::cerr << "Could not read num_bits header field." << std::endl;
    return -1;
  }

  // open index file
  RowMajorFingerprintOutputFile indexFile(outFile, bits);

  // keep track of bit counts
  std::vector<int> bitCounts;

  // process fingerprints
  int nibbles = bits / 4 + ((bits % 4) ? 1 : 0);
  while (std::getline(ifs, line)) {

    // convert hex to fingerprint
    std::pair<Word*, int> fingerprint = bitvec_from_hex(line.substr(0, nibbles));

    // record bit count
    int bitCount = bitvec_count(fingerprint.first, fingerprint.second);
    bitCounts.push_back(bitCount);

    indexFile.writeFingerprint(fingerprint.first);

    // free bit vector
    delete [] fingerprint.first;
  }

  unsigned int average_count = std::accumulate(bitCounts.begin(), bitCounts.end(), 0) / bitCounts.size();
  unsigned int min_count = *std::min_element(bitCounts.begin(), bitCounts.end());
  unsigned int max_count = *std::max_element(bitCounts.begin(), bitCounts.end());

  // create JSON header
  Json::Value data;
  data["filetype"] = "fingerprints";
  data["order"] = "row-major";
  data["num_bits"] = bits;
  data["num_fingerprints"] = static_cast<int>(bitCounts.size());
  data["fingerprint"] = Json::Value(Json::objectValue);
  data["fingerprint"]["name"] = type;
  data["fingerprint"]["type"] = type;
  data["statistics"] = Json::Value(Json::objectValue);
  data["statistics"]["average_count"] = average_count;
  data["statistics"]["min_count"] = min_count;
  data["statistics"]["max_count"] = max_count;

  // write JSON header
  Json::StyledWriter writer;
  indexFile.writeHeader(writer.write(data));
}
