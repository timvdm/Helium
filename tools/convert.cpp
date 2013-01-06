#include "../src/molecule.h"
#include "../src/fileio.h"

#include "args.h"

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

using namespace Helium;

int main(int argc, char**argv)
{
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <in_file> <out_file>" << std::endl;
    return 0;
  }

  ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("in_file", "out_file"));
  std::string inFile = args.GetArgString("in_file");
  std::string outFile = args.GetArgString("out_file");


  std::ifstream ifs(inFile.c_str());
  std::ofstream ofs(outFile.c_str(), std::ios_base::out | std::ios_base::binary);

  unsigned int numMolecules = 0;
  write32(ofs, numMolecules);

  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv(&ifs);
  conv.SetInFormat(conv.FormatFromExt(inFile));

  while (conv.Read(&mol)) {
    ++numMolecules;
    if ((numMolecules % 1000) == 0)
      std::cout << numMolecules << std::endl;
    write_molecule(ofs, &mol);
  }

  ofs.seekp(0);
  write32(ofs, numMolecules);

}
